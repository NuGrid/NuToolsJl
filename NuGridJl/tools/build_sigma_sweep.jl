# build_sigma_sweep.jl — build a ±Nσ STARLIB rate-uncertainty sweep for one
# reaction, under both starlib_option=1 (MC10+MC13) and starlib_option=2
# (MC10+MC13 with ETR25 merged in), from a single compiled template ppn/ run
# directory.
#
# Unlike build_sweep.jl's factored-rate sweep (a single flat multiplicative
# rate_index/rate_factor applied via evaluate_rates.F90's apply_rate_factors,
# constant across the whole run), this tool preserves STARLIB's own
# T9-dependent uncertainty: each STARLIB reaction file tabulates a per-T9
# "factor" column alongside the median rate (`r0_sl`/`fact_sl` in
# starlib.F90) -- the reaction's 68% (1σ) confidence interval at each T9 is
# [median/factor, median*factor]. This tool multiplies the target reaction's
# tabulated median rate by factor^sigma at every T9 point (sigma=2 -> squared,
# per the standard lognormal assumption), in a full modified COPY of the
# underlying STARLIB data file -- not the sl_flag=2/3 mechanism already in
# starlib_read_reaction_data (confirmed unused: every row in the current
# starlib_mc10_mc13_082022.txt is sl_flag=1, and starlib_append_etr25 doesn't
# even read/apply sl_flag from ETR25_starlib.dat at all) and not
# apply_rate_factors (a single flat scalar, no T9 dependence). Requires no
# Fortran changes and no recompilation -- starlib_option is a runtime
# ppn_physics.input knob, so one compiled ppn.exe serves every option/sigma
# variant, just like build_decay_sweep.jl's decay_time.
#
# Which file actually needs modifying for a given (reaction, starlib_option)
# depends on STARLIB's own internal last-write-wins merge: starlib_option=2
# appends ETR25_starlib.dat's reactions *after* MC10+MC13's, so for a
# reaction present in both, ETR25's copy is what's actually live -- verified
# directly (option-1 vs option-2 test runs of a shared, non-ETR25-only
# reaction produced bit-identical final abundances until the reaction was
# genuinely ETR25-sourced). resolve_starlib_source below picks ETR25 first
# under option 2 when the reaction is actually there, falling back to
# MC10+MC13 (matching option 1) when it isn't (e.g. Cl-32(p,g)Ar-33, which
# only MC10+MC13 has).
#
# Each run gets its own NPDATA (not just a shared symlink like build_sweep.jl
# uses): starlib.F90 hardcodes the path '../NPDATA/starlib/...' directly in
# the compiled binary (not read from any input file), and since every sigma
# variant needs a *different* modified copy of the underlying data file, they
# can't share one NPDATA the way build_sweep.jl's rate_factor-only variants
# do. So each run is nested one level deeper than build_sweep.jl's
# baseline/<reaction>/fact_<factor>/ convention: out_dir/<run_name>/ppn/ holds
# the actual ppn.exe + inputs + outputs, out_dir/<run_name>/NPDATA holds that
# run's own NPDATA (a plain symlink to the shared one for baseline runs, or a
# "shadow" tree -- symlinks to every shared file except the one modified
# reaction file, which is a real, rewritten copy -- for sigma runs).
#
# Usage:
#   julia --project=<path to NuGridJl> tools/build_sigma_sweep.jl \
#       <template_dir> <out_dir> --reaction p,f18,he4,o15 \
#       [--sigmas -2,-1,1,2] [--options 1,2] [--jobs N] [--dry-run]

function usage()
    println("""
Usage:
  julia build_sigma_sweep.jl <template_dir> <out_dir> --reaction SPECIES,... [options]

Build a ±Nσ STARLIB rate-uncertainty sweep for one reaction: for each
starlib_option in --options, a baseline_opt<N>/ run (unmodified) plus one
run_<sigma>sigma_opt<N>/ per value in --sigmas, with that reaction's tabulated
median rate multiplied by STARLIB's own per-T9 uncertainty factor raised to
the sigma power, in a rewritten copy of the underlying STARLIB data file.

Arguments:
  template_dir   Directory with a compiled ppn.exe, ppn_physics.input,
                 isotopedatabase.txt, etc. -- copied for every run.
  out_dir        Where to build baseline_opt<N>/ and run_<sigma>sigma_opt<N>/
                 (each containing a ppn/ subdirectory with the actual run).

Options:
  --reaction SPECIES,...  (required) comma-separated species tokens exactly
                          as they appear in the STARLIB reaction file, e.g.
                          "p,f18,he4,o15" for 18F(p,a)15O. Matched as a set
                          against each reaction's non-blank e1..e6 fields --
                          must resolve to exactly one row.
  --sigmas N,...          Sigma levels to run (default: -2,-1,1,2). Each N
                          multiplies the tabulated median rate by factor^N at
                          every T9 point.
  --options N,...         starlib_option values to sweep (default: 1,2)
  --jobs N, -j N           Number of ppn.exe runs in parallel (default: 4)
  --dry-run                Build directories without launching ppn.exe
  -h, --help                Show this help

Example:
  julia build_sigma_sweep.jl co_nova_1.15_10_B_mixed/ppn 18F_pa_16O \\
      --reaction p,f18,he4,o15 --jobs 8
""")
end

if abspath(PROGRAM_FILE) == (@__FILE__) && !isempty(ARGS) && ARGS[1] in ("-h", "--help")
    usage(); exit(0)
end

include(joinpath(@__DIR__, "run_parallel.jl"))
include(joinpath(@__DIR__, "namelist_utils.jl"))

using Printf

# ---------------------------------------------------------------------------
# STARLIB reaction-file parsing/rewriting
# ---------------------------------------------------------------------------
# Both starlib_mc10_mc13_082022.txt and ETR25_starlib.dat (after the -SG
# doc-comment line added at its top) share the same block structure: exactly
# one header line, then repeating 61-line blocks (1 fixed-width reaction-info
# line + 60 free-format "T9 rate factor" lines). Verified directly against
# both files' actual bytes, not assumed from documentation.

const STARLIB_NT9 = 60
const STARLIB_HEADER_LINES = 1
const OPTION1_STARLIB_FILE = "starlib_mc10_mc13_082022.txt"
const ETR25_STARLIB_FILE = "ETR25_starlib.dat"

_field(line::AbstractString, a::Integer, b::Integer) = length(line) >= b ? strip(line[a:b]) : ""

"""
    StarlibBlock

One reaction's location within a STARLIB reaction-data file: `info_line`
(1-based line number of its fixed-width info line -- chapter + up to 6
species fields; the following `STARLIB_NT9` lines are its T9/rate/factor
data) and `species` (its non-blank e1..e6 tokens, order preserved).
"""
struct StarlibBlock
    info_line::Int
    species::Vector{String}
end

"""
    parse_starlib_blocks(lines) -> Vector{StarlibBlock}

Locate every reaction block in a STARLIB file already split into lines
(`readlines`), skipping `STARLIB_HEADER_LINES`. Column offsets (chapter at
1:2, six 5-wide species fields starting at column 6) verified directly
against both starlib_mc10_mc13_082022.txt and ETR25_starlib.dat.
"""
function parse_starlib_blocks(lines::Vector{<:AbstractString})
    blocks = StarlibBlock[]
    idx = STARLIB_HEADER_LINES + 1
    n = length(lines)
    while idx <= n
        line = lines[idx]
        species = String[]
        for k in 0:5
            tok = _field(line, 6 + k * 5, 10 + k * 5)
            isempty(tok) || push!(species, tok)
        end
        push!(blocks, StarlibBlock(idx, species))
        idx += STARLIB_NT9 + 1
    end
    return blocks
end

"""
    find_reaction_block(lines, target_species) -> StarlibBlock

The one block in `lines` whose non-blank species set exactly matches
`target_species` (order-independent). Throws if zero or more than one match
-- factoring the wrong (or an ambiguous) row silently would be worse than
refusing.
"""
function find_reaction_block(lines::Vector{<:AbstractString}, target_species::Vector{String})
    target = Set(target_species)
    matches = filter(b -> Set(b.species) == target, parse_starlib_blocks(lines))
    if length(matches) != 1
        throw(ArgumentError("expected exactly 1 reaction matching $(target_species), found $(length(matches))"))
    end
    return only(matches)
end

"""
    reaction_present(path, target_species) -> Bool

Whether `target_species` resolves to exactly one row in the STARLIB file at
`path`, without throwing -- used to decide which file actually needs
factoring for a given starlib_option (see `resolve_starlib_source`).
"""
function reaction_present(path::AbstractString, target_species::Vector{String})
    isfile(path) || return false
    try
        find_reaction_block(readlines(path), target_species)
        return true
    catch err
        err isa ArgumentError && return false
        rethrow()
    end
end

"""
    factor_starlib_file!(input_path, output_path, target_species, sigma) -> StarlibBlock

Copy the STARLIB reaction file at `input_path` to `output_path`, with the one
reaction matching `target_species` (see `find_reaction_block`) rewritten: at
every T9 point, its tabulated median rate `r0` becomes `r0 * factor^sigma`
(the file's own per-T9 uncertainty `factor` column, e.g. `sigma=2` for the
standard lognormal ±2σ assumption). The reaction's info line (chapter,
species, source label, Q-value) and every other reaction are left untouched.
Returns the matched block (for logging).
"""
function factor_starlib_file!(input_path::AbstractString, output_path::AbstractString,
                               target_species::Vector{String}, sigma::Real)
    lines = readlines(input_path)
    block = find_reaction_block(lines, target_species)
    for j in 1:STARLIB_NT9
        i = block.info_line + j
        parts = split(lines[i])
        length(parts) == 3 || throw(ArgumentError("unexpected T9 data line at $input_path:$i: $(lines[i])"))
        t9, r0, fact = parse.(Float64, parts)
        new_r0 = r0 * fact^sigma
        lines[i] = @sprintf("%.6E     %.6E     %.6E", t9, new_r0, fact)
    end
    mkpath(dirname(output_path))
    open(io -> foreach(l -> println(io, l), lines), output_path, "w")
    return block
end

"""
    resolve_starlib_source(shared_starlib_dir, option, target_species) -> String

Which STARLIB reaction file actually needs factoring for `target_species`
under a given `starlib_option`: `option == 2` prefers ETR25_starlib.dat (it
wins STARLIB's own internal last-write-wins merge over MC10+MC13 for any
reaction it carries -- confirmed via source tracing and an empirical
bit-identical-abundance check), falling back to
starlib_mc10_mc13_082022.txt when the reaction isn't in ETR25 at all (e.g.
Cl-32(p,g)Ar-33), same as option 1 always uses.
"""
function resolve_starlib_source(shared_starlib_dir::AbstractString, option::Integer, target_species::Vector{String})
    if option == 2
        etr25_path = joinpath(shared_starlib_dir, ETR25_STARLIB_FILE)
        reaction_present(etr25_path, target_species) && return ETR25_STARLIB_FILE
    end
    return OPTION1_STARLIB_FILE
end

# ---------------------------------------------------------------------------
# run directory construction
# ---------------------------------------------------------------------------

const _STALE_OUTPUT_RE = r"^(iso_massf.*\.DAT|flux_.*\.DAT|x-time\.dat|OUT|fort\.6|run\.log)$"i

function _find_npdata_dir(template_dir::AbstractString)
    for candidate in (joinpath(template_dir, "NPDATA"), joinpath(template_dir, "..", "NPDATA"))
        (islink(candidate) || isdir(candidate)) && return realpath(candidate)
    end
    throw(ArgumentError("no NPDATA directory found at or above $template_dir"))
end

"""
    copy_ppn_files!(template_dir, dest_ppn_dir)

Copy every file/directory from `template_dir` into `dest_ppn_dir` (created if
needed) except `NPDATA` (handled separately, see `link_plain_npdata!`/
`build_shadow_npdata!`) and stale per-cycle output files. `ppn.exe` is
symlinked, not copied -- starlib_option is a runtime knob, no recompilation
needed for any variant this tool builds.
"""
function copy_ppn_files!(template_dir::AbstractString, dest_ppn_dir::AbstractString)
    mkpath(dest_ppn_dir)
    for name in readdir(template_dir)
        name == "NPDATA" && continue
        occursin(_STALE_OUTPUT_RE, name) && continue
        src = joinpath(template_dir, name)
        dst = joinpath(dest_ppn_dir, name)
        ispath(dst) && rm(dst; recursive = true, force = true)
        if name == "ppn.exe"
            symlink(realpath(src), dst)
        elseif islink(src)
            symlink(readlink(src), dst)
        else
            cp(src, dst; force = true)
        end
    end
end

"""
    link_plain_npdata!(shared_npdata, dest_npdata)

Point `dest_npdata` (a baseline run's `<run_name>/NPDATA`) directly at the
shared NPDATA tree -- no modification, so no reason to duplicate anything.
"""
function link_plain_npdata!(shared_npdata::AbstractString, dest_npdata::AbstractString)
    ispath(dest_npdata) && rm(dest_npdata; recursive = true, force = true)
    symlink(shared_npdata, dest_npdata)
end

"""
    build_shadow_npdata!(shared_npdata, dest_npdata, modified_file_name, modified_file_path)

Build `dest_npdata` (a sigma run's `<run_name>/NPDATA`) as a real directory
symlinking every entry of the shared NPDATA tree, except that `starlib/` is
itself a real directory symlinking every file in the shared `NPDATA/starlib/`
except `modified_file_name`, which is a real copy of the already-factored
file at `modified_file_path` (see `factor_starlib_file!`). Every other data
file (reaclib tables, screening data, VITAL's own inputs, ...) stays a plain
symlink to the shared copy -- only the one reaction file this run actually
modifies is duplicated.
"""
function build_shadow_npdata!(shared_npdata::AbstractString, dest_npdata::AbstractString,
                               modified_file_name::AbstractString, modified_file_path::AbstractString)
    ispath(dest_npdata) && rm(dest_npdata; recursive = true, force = true)
    mkpath(dest_npdata)
    for name in readdir(shared_npdata)
        name == "starlib" && continue
        symlink(joinpath(shared_npdata, name), joinpath(dest_npdata, name))
    end
    shared_starlib = joinpath(shared_npdata, "starlib")
    dest_starlib = joinpath(dest_npdata, "starlib")
    mkpath(dest_starlib)
    for name in readdir(shared_starlib)
        dst = joinpath(dest_starlib, name)
        if name == modified_file_name
            cp(modified_file_path, dst; force = true)
        else
            symlink(joinpath(shared_starlib, name), dst)
        end
    end
end

"""
    build_run!(template_dir, run_dir, option; sigma_spec = nothing) -> String

Build one run at `run_dir/ppn/` (files from `template_dir`, `starlib_option`
patched to `option` in `ppn_physics.input`) with `run_dir/NPDATA` either a
plain symlink to the shared NPDATA (`sigma_spec === nothing`, i.e. a baseline
run) or a shadow tree with one reaction's rate factored by
`sigma_spec.sigma` (`sigma_spec = (target_species = [...], sigma = N)`).
Returns the `ppn/` directory (what `run_parallel` should actually launch
`ppn.exe` in).
"""
function build_run!(template_dir::AbstractString, run_dir::AbstractString, option::Integer;
                     sigma_spec::Union{Nothing,NamedTuple} = nothing)
    shared_npdata = _find_npdata_dir(template_dir)
    ppn_dir = joinpath(run_dir, "ppn")
    copy_ppn_files!(template_dir, ppn_dir)

    dest_npdata = joinpath(run_dir, "NPDATA")
    if sigma_spec === nothing
        link_plain_npdata!(shared_npdata, dest_npdata)
    else
        shared_starlib = joinpath(shared_npdata, "starlib")
        source_file = resolve_starlib_source(shared_starlib, option, sigma_spec.target_species)
        factored_path = joinpath(run_dir, "_factored_" * source_file)
        factor_starlib_file!(joinpath(shared_starlib, source_file), factored_path,
                              sigma_spec.target_species, sigma_spec.sigma)
        build_shadow_npdata!(shared_npdata, dest_npdata, source_file, factored_path)
        rm(factored_path)
    end

    physics_input = joinpath(ppn_dir, "ppn_physics.input")
    write(physics_input, update_namelist(read(physics_input, String), ["starlib_option" => option]))

    return ppn_dir
end

# ---------------------------------------------------------------------------
# top level
# ---------------------------------------------------------------------------

"""
    sigma_label(sigma) -> String

Directory-name fragment for a sigma level (`1` -> `"1sigma"`, `-2` ->
`"-2sigma"`).
"""
sigma_label(sigma::Real) = string(isinteger(sigma) ? Int(sigma) : sigma) * "sigma"

"""
    build_sigma_sweep(template_dir, out_dir, target_species;
                       sigmas = [-2, -1, 1, 2], options = [1, 2],
                       jobs = 4, dry_run = false) -> Vector{String}

Build `out_dir/baseline_opt<N>/ppn/` and `out_dir/run_<sigma>sigma_opt<N>/ppn/`
for every `N` in `options` and `sigma` in `sigmas` (see `build_run!`), then
(unless `dry_run`) run every one of them via `run_parallel`. Returns the
`ppn/` directories built, in build order.
"""
function build_sigma_sweep(template_dir::AbstractString, out_dir::AbstractString,
                            target_species::Vector{String}; sigmas::Vector{<:Real} = [-2, -1, 1, 2],
                            options::Vector{<:Integer} = [1, 2], jobs::Integer = 4, dry_run::Bool = false)
    run_dirs = String[]
    for option in options
        baseline_dir = joinpath(out_dir, "baseline_opt$option")
        push!(run_dirs, build_run!(template_dir, baseline_dir, option))

        for sigma in sigmas
            run_dir = joinpath(out_dir, "run_$(sigma_label(sigma))_opt$option")
            push!(run_dirs, build_run!(template_dir, run_dir, option;
                                        sigma_spec = (target_species = target_species, sigma = sigma)))
        end
    end

    if !dry_run
        results = run_parallel(run_dirs; jobs)
        failed = [dir for (dir, ok) in results if !ok]
        isempty(failed) || @warn "some sigma-sweep runs failed" failed
    end

    return run_dirs
end

function parse_int_list(text::AbstractString)
    return [parse(Int, strip(s)) for s in split(text, ',') if !isempty(strip(s))]
end

function parse_real_list(text::AbstractString)
    return [parse(Float64, strip(s)) for s in split(text, ',') if !isempty(strip(s))]
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    let jobs = 4, dry_run = false, reaction_arg = nothing,
        sigmas = Float64[-2, -1, 1, 2], options = Int[1, 2], positional = String[], i = 1

        if isempty(ARGS) || ARGS[1] in ("-h", "--help")
            usage(); exit(0)
        end
        while i <= length(ARGS)
            a = ARGS[i]
            if a in ("-h", "--help")
                usage(); exit(0)
            elseif a == "--reaction"
                i == length(ARGS) && error("--reaction requires a value")
                reaction_arg = String.(strip.(split(ARGS[i + 1], ','))); i += 2
            elseif a == "--sigmas"
                i == length(ARGS) && error("--sigmas requires a value")
                sigmas = parse_real_list(ARGS[i + 1]); i += 2
            elseif a == "--options"
                i == length(ARGS) && error("--options requires a value")
                options = parse_int_list(ARGS[i + 1]); i += 2
            elseif a in ("--jobs", "-j")
                i == length(ARGS) && error("--jobs requires a value")
                jobs = parse(Int, ARGS[i + 1]); i += 2
            elseif a == "--dry-run"
                dry_run = true; i += 1
            elseif startswith(a, "-")
                println(stderr, "Unknown option: $a\n")
                usage(); exit(1)
            else
                push!(positional, a); i += 1
            end
        end
        if length(positional) != 2
            println(stderr, "Expected 2 positional arguments (<template_dir> <out_dir>), got $(length(positional)).\n")
            usage(); exit(1)
        end
        reaction_arg === nothing && (println(stderr, "--reaction is required.\n"); usage(); exit(1))

        template_dir, out_dir = positional
        built = build_sigma_sweep(template_dir, out_dir, reaction_arg; sigmas, options, jobs, dry_run)
        println(length(built), " run(s) built at ", out_dir)
    end
end
