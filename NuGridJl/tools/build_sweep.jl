# build_sweep.jl — build an Iliadis (2002)-style one-reaction-at-a-time rate
# sweep from a template ppn/ run directory.
#
# Ports the reaction-resolution logic from NovaSensitivityStudy's
# single-zone/tools/NovaRunTools.jl (`resolve_reaction_index`,
# `matching_rows_for_reaction`, `validate_reverse_index`, `copy_ppn`) onto
# NuGridJl's typed `Network`/`Reaction`/`Isotope` instead of that module's
# own hand-rolled `Row`, and drops the nova_cases/-specific path assumptions
# so it works on any template directory.
#
# Also ports NovaRunTools.jl's `source_priority` (see SOURCE_PRESETS below)
# as an optional, named `preset` -- the reaction-type/target-mass-dependent
# priority ladder Iliadis et al. 2002 imply (NACRE below A=20, ILI01 from
# 20-40, its own ladders for weak and (n,g) reactions), for when the plain
# `--priorities` flat list isn't reaction-aware enough to pick a unique
# active row. Deliberately NOT ported: NovaRunTools.jl's
# `enforce_activation_policy`/`network_edits.json` machinery -- that's
# network *curation* (making sure each physical reaction has exactly one
# active row before a sweep is ever built), a separate concern from
# sweep-building itself.
#
# Not part of the NuGridJl package: this is the one place allowed to write
# ppn_physics.input rate overrides and launch ppn.exe. NuGridJl itself only
# ever reads what this script (or anything following the same directory
# convention) produces, via `PPNSweep`.
#
# Usage:
#   julia --project=<path to NuGridJl> tools/build_sweep.jl \
#       <template_dir> <reaction_plan.json> <out_dir> \
#       [--jobs N] [--dry-run] [--priorities SRC,...] [--preset NAME]

using NuGridJl
using JSON

include("run_parallel.jl")

# ---------------------------------------------------------------------------
# reaction plan
# ---------------------------------------------------------------------------
# A JSON array (or `{"reactions": [...]}`) of entries:
#   {"name": "13N_pg_14O", "factors": [0.5, 2.0]}
#   {"name": "23Na_pa_20Ne", "factors": [...], "reverse_index": 31}   (optional)
#   {"name": "...", "factors": [...], "index": 154}                   (optional, forces the row)
# `name` is `<target><channel><product>` (isomer products end in g/m),
# matching the convention already used throughout NovaSensitivityStudy's
# reaction_plan.json files and sweep directory names.

struct ReactionPlanEntry
    name::String
    factors::Vector{Float64}
    index::Union{Nothing,Int}          # force this network index instead of resolving by name
    reverse_index::Union{Nothing,Int}  # also factor this reverse reaction, by the same amount
end

function read_reaction_plan(path::AbstractString; default_factors::Vector{Float64} = Float64[])
    data = JSON.parsefile(path)
    entries = data isa AbstractDict ? data["reactions"] : data
    plan_defaults = data isa AbstractDict ? Float64.(get(data, "default_factors", default_factors)) : default_factors
    return [
        ReactionPlanEntry(
            e["name"],
            Float64.(get(e, "factors", plan_defaults)),
            haskey(e, "index") ? Int(e["index"]) : nothing,
            haskey(e, "reverse_index") ? Int(e["reverse_index"]) : nothing,
        )
        for e in entries
    ]
end

# ---------------------------------------------------------------------------
# reaction name -> Network Reaction, with NovaRunTools.jl's ambiguity rules
# ---------------------------------------------------------------------------

const _CHANNEL_RTYPE = Dict(
    "pg" => "(p,g)", "pa" => "(p,a)", "pn" => "(p,n)", "pp" => "(p,p)",
    "ag" => "(a,g)", "an" => "(a,n)", "ap" => "(a,p)",
    "ng" => "(n,g)", "na" => "(n,a)", "np" => "(n,p)",
)
const _ALPHA_REVERSE_RTYPES = Set(["(p,a)", "(a,p)", "(n,a)", "(a,n)"])
const _REVERSE_RTYPE = Dict("(p,a)" => "(a,p)", "(a,p)" => "(p,a)", "(n,a)" => "(a,n)", "(a,n)" => "(n,a)")

# The light species a channel letter refers to -- used to pick the target out
# of the reactant list by *role*, not by mass. Mass ordering breaks for the
# one case a light projectile outweighs its target (3He(a,g)7Be: the alpha
# projectile, A=4, is heavier than the 3He target, A=3) -- 'g' (gamma) has no
# tracked species, meaning a single-reactant/single-product reaction.
const _CHANNEL_SPECIES = Dict('p' => Isotope(1, 1, 0), 'a' => Isotope(2, 4, 0), 'n' => Isotope(0, 1, 0))

function _parse_reaction_name(name::AbstractString)
    m = match(r"^(\d+[A-Za-z]+[gm]?)_([a-z]{2})_(\d+[A-Za-z]+[gm]?)$", name)
    m === nothing && throw(ArgumentError("reaction name must look like \"13N_pg_14O\" (got \"$name\")"))
    channel = m.captures[2]
    rtype = get(_CHANNEL_RTYPE, channel, nothing)
    rtype === nothing && throw(ArgumentError("unknown channel code \"$channel\" in \"$name\""))
    return (target = Isotope(String(m.captures[1])), rtype = rtype, product = Isotope(String(m.captures[3])),
            projectile = get(_CHANNEL_SPECIES, channel[1], nothing), ejectile = get(_CHANNEL_SPECIES, channel[2], nothing))
end

function _reactants_match(r::Reaction, target::Isotope, projectile::Union{Nothing,Isotope})
    projectile === nothing && return length(r.reactants) == 1 && only(r.reactants) == target
    return length(r.reactants) == 2 && target in r.reactants && projectile in r.reactants
end

function _products_match(r::Reaction, product::Isotope, ejectile::Union{Nothing,Isotope})
    ejectile === nothing && return length(r.products) == 1 && only(r.products) == product
    return length(r.products) == 2 && product in r.products && ejectile in r.products
end

"""
    matching_reactions(net::Network, name::AbstractString) -> Vector{Reaction}

Every reaction in `net` matching `name` (e.g. `"13N_pg_14O"`, or
`"3He_ag_7Be"` where the alpha projectile outweighs the 3He target):
reactants are exactly `{target, projectile}` (or just `{target}` for a
`g`-initiated channel), reaction-type code matches, and products are exactly
`{product, ejectile}` (or just `{product}`) — or, if nothing matches on
product too (the network's own boundary remapping can redirect an
out-of-range product to a different species), every reaction matching just
target+projectile+rtype.
"""
function matching_reactions(net::Network, name::AbstractString)
    t = _parse_reaction_name(name)
    same_target_rtype = filter(net.reactions) do r
        r.rtype == t.rtype && _reactants_match(r, t.target, t.projectile)
    end
    exact = filter(r -> _products_match(r, t.product, t.ejectile), same_target_rtype)
    return isempty(exact) ? same_target_rtype : exact
end

const _WEAK_RTYPES = Set(["(-,g)", "(+,g)"])
const _NOVA_MAX_A = 40

"Reaction-aware source priority ladder for the `\"iliadis2002\"` preset — see [`SOURCE_PRESETS`](@ref)."
function _iliadis2002_source_priority(rtype::AbstractString, target_a::Integer)
    if rtype in _WEAK_RTYPES
        return ["ODA94", "NETB1", "FFW85", "LMP00", "JINAC", "BASEL", "JINAR", "JINAV"]
    elseif rtype == "(n,g)"
        return target_a < 20 ?
            ["NACRR", "NACRL", "NACRU", "KADON", "JINAC", "BASEL", "JINAR", "JINAV"] :
            ["ILI01", "NACRR", "NACRL", "NACRU", "KADON", "JINAC", "BASEL", "JINAR", "JINAV"]
    elseif 1 <= target_a < 20
        return ["NACRR", "NACRL", "NACRU", "JINAC", "BASEL", "JINAR", "JINAV", "VITAL", "RVRSE"]
    elseif 20 <= target_a <= _NOVA_MAX_A
        return ["ILI01", "JINAC", "BASEL", "JINAR", "JINAV", "NACRR", "NACRL", "NACRU", "KADON", "VITAL", "RVRSE"]
    else
        return ["JINAC", "BASEL", "JINAR", "JINAV", "KADON", "RVRSE", "VITAL"]
    end
end

"""
    SOURCE_PRESETS

Named, reaction-aware rate-source priority ladders — an alternative to the
flat `--priorities` list for when disambiguation genuinely depends on the
reaction's type and target mass, not just a single fixed source order.
Each entry is a function `(rtype, target_a) -> Vector{String}` (most
preferred first), tried by [`resolve_reaction`](@ref) via the `preset`
keyword. Add new presets here as the authors' own preferred configurations
are worked out — this is meant to grow, not just hold one entry.

`"iliadis2002"` ports NovaRunTools.jl's `source_priority`: NACRE below
A=20, ILI01 (Iliadis et al. 2002's own evaluation) from 20-40 for the
charged-particle "controlled" types, with separate ladders for weak
((-,g)/(+,g)) and (n,g) reactions, matching the source choices implied by
Iliadis et al. 2002's own network.
"""
const SOURCE_PRESETS = Dict{String,Function}(
    "iliadis2002" => _iliadis2002_source_priority,
)

"""
    resolve_reaction(net::Network, entry::ReactionPlanEntry; prefer_sources = String[],
                      preset = nothing) -> Reaction

Resolve `entry.name` to one `Reaction` in `net`, following the same rules as
NovaRunTools.jl's `resolve_reaction_index`:
- if `entry.index` is set, use that row (falling back to an unambiguous
  active alternative and warning, if the forced row turns out inactive);
- otherwise prefer an active row; if more than one is active, the priority
  list — `prefer_sources` tried first, then `preset`'s (a key into
  [`SOURCE_PRESETS`](@ref)) reaction-aware ladder if given — must narrow it
  to exactly one, or this throws; ambiguity is never silently guessed at;
- if no row is active, fall back to the single candidate if there's only one,
  otherwise throws.
"""
function resolve_reaction(net::Network, entry::ReactionPlanEntry; prefer_sources::Vector{String} = String[],
                           preset::Union{Nothing,AbstractString} = nothing)
    candidates = matching_reactions(net, entry.name)
    isempty(candidates) && throw(ArgumentError("$(entry.name): no matching row in networksetup.txt"))
    priorities = _resolve_priorities(entry.name, prefer_sources, preset)

    if entry.index !== nothing
        forced = filter(r -> r.index == entry.index, candidates)
        if !isempty(forced)
            found = only(forced)
            found.active && return found
            active = filter(r -> r.active, candidates)
            chosen = _preferred_active(active, priorities)
            if chosen !== nothing
                @warn "$(entry.name): configured index $(entry.index) is inactive; using $(chosen.index) ($(chosen.source)) instead"
                return chosen
            end
            return found  # no unambiguous active alternative; use the forced (inactive) row as-is
        end
        @warn "$(entry.name): configured index $(entry.index) not found; auto-resolving from $(length(candidates)) candidate(s)"
    end

    active = filter(r -> r.active, candidates)
    chosen = _preferred_active(active, priorities)
    chosen !== nothing && return chosen
    if length(active) > 1
        options = join(("$(r.index) ($(r.source))" for r in active), ", ")
        throw(ArgumentError("$(entry.name): ambiguous — $(length(active)) active rows ($options) — " *
                              "add an explicit \"index\" to the reaction plan, or pass --priorities/--preset"))
    elseif length(active) == 1
        return only(active)
    elseif length(candidates) == 1
        @warn "$(entry.name): selected row is inactive — factoring it may have no effect"
        return only(candidates)
    else
        options = join(("$(r.index) ($(r.source))" for r in candidates), ", ")
        throw(ArgumentError("$(entry.name): ambiguous — $(length(candidates)) inactive rows ($options) — " *
                              "add an explicit \"index\" to the reaction plan"))
    end
end

"Effective priority list for `name`: `prefer_sources` tried first, then `preset`'s reaction-aware ladder (if given)."
function _resolve_priorities(name::AbstractString, prefer_sources::Vector{String}, preset::Union{Nothing,AbstractString})
    preset === nothing && return prefer_sources
    ladder = get(SOURCE_PRESETS, preset, nothing)
    ladder === nothing && throw(ArgumentError("unknown preset \"$preset\" — known presets: $(join(keys(SOURCE_PRESETS), ", "))"))
    t = _parse_reaction_name(name)
    return vcat(prefer_sources, ladder(t.rtype, t.target.A))
end

function _preferred_active(active::Vector{Reaction}, prefer_sources::Vector{String})
    isempty(active) && return nothing
    length(active) == 1 && return only(active)
    for src in prefer_sources
        matched = filter(r -> r.source == src, active)
        length(matched) == 1 && return only(matched)
    end
    return nothing
end

"""
    resolve_reverse(net::Network, entry::ReactionPlanEntry, forward::Reaction) -> Union{Nothing,Reaction}

If `entry.reverse_index` is set, look it up and validate it's actually the
reverse of `forward`: an alpha-transfer type (`(p,a)`/`(a,p)`/`(n,a)`/`(a,n)`)
whose type is `forward`'s reverse, with reactant/product sets swapped.
"""
function resolve_reverse(net::Network, entry::ReactionPlanEntry, forward::Reaction)
    entry.reverse_index === nothing && return nothing
    reverse = findfirst(r -> r.index == entry.reverse_index, net.reactions)
    reverse === nothing && throw(ArgumentError("$(entry.name): reverse_index $(entry.reverse_index) not found"))
    r = net.reactions[reverse]
    forward.rtype in _ALPHA_REVERSE_RTYPES || throw(ArgumentError(
        "$(entry.name): reverse_index is only supported for alpha-transfer reactions, not $(forward.rtype)"))
    r.rtype == _REVERSE_RTYPE[forward.rtype] || throw(ArgumentError(
        "$(entry.name): reverse_index $(entry.reverse_index) has type $(r.rtype), expected $(_REVERSE_RTYPE[forward.rtype])"))
    (Set(r.reactants) == Set(forward.products) && Set(r.products) == Set(forward.reactants)) || throw(ArgumentError(
        "$(entry.name): reverse_index $(entry.reverse_index) is not the reverse of index $(forward.index)"))
    return r
end

# ---------------------------------------------------------------------------
# building run directories
# ---------------------------------------------------------------------------

const _STALE_OUTPUT_RE = r"^(iso_massf.*\.DAT|flux_.*\.DAT|x-time\.dat|OUT|fort\.6)$"i

"""
    copy_template!(template_dir, dest_dir)

Copy every file/directory from `template_dir` into `dest_dir` (created if
needed), skipping stale per-cycle output files so a run that doesn't
actually execute can't be mistaken for one that did. `NPDATA` is linked, not
copied (a full copy would duplicate a multi-GB rate library per run): if
`template_dir` (or its parent) has an `NPDATA` symlink, both `dest_dir/NPDATA`
and `dirname(dest_dir)/NPDATA` are symlinked at its real target, matching
`ppn_physics.input`'s `../NPDATA/...` relative paths regardless of whether
`dest_dir` sits directly under `out_dir` (baseline/) or one level deeper
(`<reaction>/fact_<factor>/`).
"""
function copy_template!(template_dir::AbstractString, dest_dir::AbstractString)
    mkpath(dest_dir)
    npdata_target = _find_npdata_target(template_dir)
    for name in readdir(template_dir)
        name == "NPDATA" && continue
        occursin(_STALE_OUTPUT_RE, name) && continue
        src = joinpath(template_dir, name)
        dst = joinpath(dest_dir, name)
        ispath(dst) && rm(dst; recursive = true, force = true)
        if islink(src)
            symlink(readlink(src), dst)
        else
            cp(src, dst; force = true)
        end
    end
    if npdata_target !== nothing
        for npdata_link in (joinpath(dest_dir, "NPDATA"), joinpath(dirname(dest_dir), "NPDATA"))
            ispath(npdata_link) || symlink(npdata_target, npdata_link)
        end
    end
end

function _find_npdata_target(template_dir::AbstractString)
    for candidate in (joinpath(template_dir, "NPDATA"), joinpath(template_dir, "..", "NPDATA"))
        islink(candidate) && return realpath(candidate)
        isdir(candidate) && return realpath(candidate)
    end
    return nothing
end

"""
    write_rate_factors!(ppn_physics_input_path, index_factor_pairs)

Insert `rate_index(i) = <index>` / `rate_factor(i) = <factor>` for each
`(index, factor)` pair, just before the `&ppn_physics` namelist's closing
`/`. Slots are numbered from 1 (NuPPN supports up to `num_rate_factors = 10`).
"""
function write_rate_factors!(ppn_physics_input_path::AbstractString, index_factor_pairs)
    lines = readlines(ppn_physics_input_path)
    terminator = findfirst(l -> strip(l) == "/", lines)
    terminator === nothing && throw(ArgumentError(
        "no namelist terminator '/' found in $ppn_physics_input_path"))
    new_lines = String[]
    for (i, (index, factor)) in enumerate(index_factor_pairs)
        push!(new_lines, "        rate_index($i) = $index")
        push!(new_lines, "        rate_factor($i) = $(factor)")
    end
    splice!(lines, terminator:(terminator - 1), new_lines)
    write(ppn_physics_input_path, join(lines, "\n") * "\n")
end

# ---------------------------------------------------------------------------
# top level
# ---------------------------------------------------------------------------

"""
    build_sweep(template_dir, reaction_plan_path, out_dir;
                jobs = 4, dry_run = false, prefer_sources = String[],
                preset = nothing) -> PPNSweep

Build a factored-rate sweep at `out_dir`: copy `template_dir` to
`out_dir/baseline/` and run it once (its `networksetup.txt` is the
authoritative reaction ordering every factored run's index is resolved
against — rate-set changes can shift indices, so this is never skipped),
then for every `(reaction, factor)` pair in the reaction plan at
`reaction_plan_path`, copy `template_dir` again to
`out_dir/<reaction>/fact_<factor>/` with that reaction's (and, if configured,
its `reverse_index`'s) `rate_index`/`rate_factor` written in. With
`dry_run = true`, directories are built and namelists written but `ppn.exe`
is never launched — useful for checking a reaction plan resolves before
committing real compute to it.

`prefer_sources`/`preset` are forwarded to [`resolve_reaction`](@ref) for
every entry — `preset` (a key into [`SOURCE_PRESETS`](@ref), e.g.
`"iliadis2002"`) picks a reaction-aware source ladder instead of (or ahead
of) a flat `prefer_sources` list.
"""
function build_sweep(template_dir::AbstractString, reaction_plan_path::AbstractString,
                      out_dir::AbstractString; jobs::Integer = 4, dry_run::Bool = false,
                      prefer_sources::Vector{String} = String[],
                      preset::Union{Nothing,AbstractString} = nothing)
    plan = read_reaction_plan(reaction_plan_path)
    baseline_dir = joinpath(out_dir, "baseline")

    copy_template!(template_dir, baseline_dir)
    dry_run || only(run_parallel([baseline_dir]; jobs = 1))[2] ||
        error("baseline run failed in $baseline_dir; see $(joinpath(baseline_dir, "run.log"))")

    net = network(PPNRun(baseline_dir))

    run_dirs = String[]
    unresolved = Tuple{String,String}[]  # (entry.name, error message)
    for entry in plan
        local reaction, reverse
        try
            reaction = resolve_reaction(net, entry; prefer_sources, preset)
            reverse = resolve_reverse(net, entry, reaction)
        catch err
            err isa ArgumentError || rethrow()
            push!(unresolved, (entry.name, err.msg))
            @warn "$(entry.name): could not resolve, skipping" exception = err.msg
            continue
        end
        pairs = reverse === nothing ? [(reaction.index, 0.0)] : [(reaction.index, 0.0), (reverse.index, 0.0)]
        for factor in entry.factors
            run_dir = joinpath(out_dir, entry.name, "fact_$(factor)")
            copy_template!(template_dir, run_dir)
            write_rate_factors!(joinpath(run_dir, "ppn_physics.input"),
                                  [(idx, factor) for (idx, _) in pairs])
            push!(run_dirs, run_dir)
        end
    end

    if !isempty(unresolved)
        println("$(length(unresolved))/$(length(plan)) reaction(s) could not be resolved and were skipped:")
        for (name, msg) in unresolved
            println("  $name: $msg")
        end
    end

    if !dry_run
        results = run_parallel(run_dirs; jobs)
        failed = [dir for (dir, ok) in results if !ok]
        isempty(failed) || @warn "some factored runs failed" failed
    end

    return PPNSweep(out_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    let jobs = 4, dry_run = false, prefer_sources = String[], preset = nothing, positional = String[], i = 1
        while i <= length(ARGS)
            a = ARGS[i]
            if a == "--jobs"
                jobs = parse(Int, ARGS[i + 1]); i += 2
            elseif a == "--dry-run"
                dry_run = true; i += 1
            elseif a == "--priorities"
                prefer_sources = String.(split(ARGS[i + 1], ',')); i += 2
            elseif a == "--preset"
                preset = ARGS[i + 1]; i += 2
            else
                push!(positional, a); i += 1
            end
        end
        length(positional) == 3 || error(
            "usage: julia build_sweep.jl <template_dir> <reaction_plan.json> <out_dir> " *
            "[--jobs N] [--dry-run] [--priorities SRC,...] [--preset NAME]")
        sweep = build_sweep(positional...; jobs, dry_run, prefer_sources, preset)
        println(sweep)
    end
end
