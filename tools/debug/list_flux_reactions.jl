#!/usr/bin/env julia
# list_flux_reactions.jl
#
# List reactions that contributed non-trivial flux during a nova trajectory,
# ranked by peak |dY/dt|. Reads all flux_*.DAT files from the baseline run
# and joins with networksetup.txt for reaction notation and rate source.
#
# Usage:
#   julia tools/list_flux_reactions.jl --nova ne_nova_1.15_12_X_weiss_mixed
#   julia tools/list_flux_reactions.jl --nova ne_nova_1.15_12_X_weiss_mixed --run runs_star_opt2_novital
#   julia tools/list_flux_reactions.jl --help

using Printf

# ── Paths ─────────────────────────────────────────────────────────────────────

const NOVA_CASES = joinpath(@__DIR__, "..", "nova_cases")

# ── Source descriptions ────────────────────────────────────────────────────────

const SOURCE_NAMES = Dict(
    "STL02" => "STARLIB (Sallaska et al. 2013, ETR25 Monte Carlo rates)",
    "NACRR" => "NACRE II recommended (Xu et al. 2013)",
    "NACRL" => "NACRE II lower limit",
    "NACRU" => "NACRE II upper limit",
    "ILI01" => "Iliadis et al. 2001 (shell-model, charged particle)",
    "JINAC" => "JINA REACLIB v2 (Cyburt et al. 2010)",
    "JINAR" => "JINA REACLIB v2 reverse",
    "JINAV" => "JINA REACLIB v2 (VITAL supplement)",
    "RVRSE" => "Reverse rate via detailed balance",
    "VITAL" => "VITAL internal rate (vital.F90; pp-chain / compound)",
    "NETB1" => "Nuclear structure (beta decays, Nubase/ENSDF)",
    "ODA94" => "Oda et al. 1994 (weak rates)",
    "LMP00" => "Langanke & Martínez-Pinedo 2000 (stellar weak rates)",
    "FFW85" => "Fuller, Fowler & Newman 1985 (stellar weak rates)",
    "KADON" => "KADoNiS neutron capture database",
    "BASEL" => "Basel (Hauser-Feshbach, Rauscher & Thielemann 2000)",
)

# ── Channel notation ───────────────────────────────────────────────────────────

const RTYPE_CHANNEL = Dict(
    "(p,g)"  => ("p",  "γ"),
    "(p,a)"  => ("p",  "α"),
    "(p,n)"  => ("p",  "n"),
    "(a,g)"  => ("α",  "γ"),
    "(a,p)"  => ("α",  "p"),
    "(a,n)"  => ("α",  "n"),
    "(n,g)"  => ("n",  "γ"),
    "(n,p)"  => ("n",  "p"),
    "(n,a)"  => ("n",  "α"),
    "(g,p)"  => ("γ",  "p"),
    "(g,a)"  => ("γ",  "α"),
    "(g,n)"  => ("γ",  "n"),
    "(+,g)"  => ("β+", "ν"),
    "(-,g)"  => ("β-", "ν̄"),
    "(e,v)"  => ("e-", "ν"),
)

function parse_sp(sp::AbstractString)
    sp = strip(sp)
    sp == "OOOOO" && return nothing
    sp == "NEUT"  && return (0, "n")
    sp == "PROT"  && return (1, "p")
    m = match(r"^([A-Z]+)\s+(\d+)$", sp)
    m === nothing && return nothing
    sym_raw, mass = m.captures[1], parse(Int, m.captures[2])
    sym = length(sym_raw) > 1 ? sym_raw[1] * lowercase(sym_raw[2:end]) : sym_raw
    return (mass, sym)
end

function fmt_sp(sp::AbstractString)
    p = parse_sp(sp)
    p === nothing && return ""
    mass, sym = p
    sym == "n" && return "n"
    sym == "p" && mass == 1 && return "p"
    return "$(mass)$(sym)"
end

function build_notation(n1, sp1, n2, sp2, n3, sp3, n4, sp4, rtype)
    ch = get(RTYPE_CHANNEL, rtype, nothing)

    if ch === nothing
        lhs = join(filter(!isempty, [
            (parse(Int, n1) > 0 && strip(sp1) != "OOOOO") ?
                (parse(Int, n1) > 1 ? "$(n1)×$(fmt_sp(sp1))" : fmt_sp(sp1)) : "",
            (parse(Int, n2) > 0 && strip(sp2) != "OOOOO") ?
                (parse(Int, n2) > 1 ? "$(n2)×$(fmt_sp(sp2))" : fmt_sp(sp2)) : "",
        ]), " + ")
        rhs = join(filter(!isempty, [
            (parse(Int, n3) > 0 && strip(sp3) != "OOOOO") ?
                (parse(Int, n3) > 1 ? "$(n3)×$(fmt_sp(sp3))" : fmt_sp(sp3)) : "",
            (parse(Int, n4) > 0 && strip(sp4) != "OOOOO") ?
                (parse(Int, n4) > 1 ? "$(n4)×$(fmt_sp(sp4))" : fmt_sp(sp4)) : "",
        ]), " + ")
        return "$(isempty(lhs) ? "?" : lhs) → $(isempty(rhs) ? "?" : rhs)"
    end

    proj, eject = ch
    target = fmt_sp(sp1)

    if proj == "γ"
        residual = fmt_sp(sp3); isempty(residual) && (residual = fmt_sp(sp4))
        return "$(target)(γ,$(eject))$(residual)"
    end
    if proj in ("β+", "β-", "e-")
        return "$(target)($(proj),$(eject))$(fmt_sp(sp3))"
    end
    residual = fmt_sp(sp3); isempty(residual) && (residual = fmt_sp(sp4))
    return "$(target)($(proj),$(eject))$(residual)"
end

# ── networksetup.txt parsing ───────────────────────────────────────────────────

function parse_fortran_float(s::AbstractString)
    v = tryparse(Float64, s)
    v !== nothing && return v
    m = match(r"^([+-]?\d+\.?\d*)([+-]\d+)$", s)
    m !== nothing && return parse(Float64, "$(m.captures[1])E$(m.captures[2])")
    return 0.0
end

const NETWORK_RE = r"^\s*(\d+)\s+([TF])\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+->\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)"

struct RxnMeta
    index    :: Int
    active   :: Bool
    Q_MeV    :: Float64
    source   :: String
    rtype    :: String
    notation :: String
end

function parse_networksetup(path::AbstractString)
    meta = Dict{Int, RxnMeta}()
    for line in eachline(path)
        m = match(NETWORK_RE, line)
        m === nothing && continue
        idx, act, n1, sp1, n2, sp2, n3, sp3, n4, sp4, q_str, src, rtype, _ = m.captures
        notation = build_notation(n1, sp1, n2, sp2, n3, sp3, n4, sp4, rtype)
        i = parse(Int, idx)
        meta[i] = RxnMeta(i, act == "T", parse_fortran_float(q_str), src, rtype, notation)
    end
    return meta
end

# ── Flux file parsing ──────────────────────────────────────────────────────────

# flux_*.DAT format (space-delimited):
#   idx  Z_k1 A_k1 Z_k3 A_k3  Z_k5 A_k5 Z_k7 A_k7  flux[dY/dt]  energy[erg/(g*s)]  timescale[s]
# Sentinel value 1e-99 means no flux; 1e+90 means infinite timescale.

const FLUX_SENTINEL = 1e-20   # flux below this is treated as zero
const ENERGY_SENTINEL = 1e+80 # timescale above this is effectively infinite

struct FluxStats
    max_flux   :: Float64   # peak |dY/dt| across all snapshots
    max_energy :: Float64   # peak |erg/(g*s)| across all snapshots
    n_active   :: Int       # number of snapshots with |flux| > FLUX_SENTINEL
end

function read_flux_files(flux_dir::AbstractString)
    files = sort(filter(f -> startswith(f, "flux_") && endswith(f, ".DAT"),
                        readdir(flux_dir)))
    isempty(files) && error("No flux_*.DAT files found in $flux_dir")

    stats = Dict{Int, Tuple{Float64, Float64, Int}}()  # idx → (max_flux, max_energy, n_active)

    for fname in files
        path = joinpath(flux_dir, fname)
        for line in eachline(path)
            startswith(strip(line), "#") && continue
            parts = split(line)
            length(parts) < 11 && continue
            idx_v = tryparse(Int, parts[1])
            idx_v === nothing && continue
            flux   = parse_fortran_float(parts[10])
            energy = parse_fortran_float(parts[11])
            aflux  = abs(flux)
            aenergy = abs(energy)

            prev = get(stats, idx_v, (0.0, 0.0, 0))
            n_active = prev[3] + (aflux > FLUX_SENTINEL ? 1 : 0)
            stats[idx_v] = (max(prev[1], aflux), max(prev[2], aenergy), n_active)
        end
    end

    return Dict(k => FluxStats(v...) for (k, v) in stats)
end

# ── Output ─────────────────────────────────────────────────────────────────────

function write_csv(rows, path, n_snapshots)
    open(path, "w") do f
        println(f, "index,reaction,rtype,Q_MeV,source,source_description,max_flux_dYdt,max_energy_erg_g_s,n_active_snapshots,fraction_active")
        for (meta, fstats) in rows
            desc = get(SOURCE_NAMES, meta.source, meta.source)
            frac = fstats.n_active / n_snapshots
            println(f,
                "$(meta.index),\"$(meta.notation)\",$(meta.rtype),$(meta.Q_MeV)," *
                "$(meta.source),\"$(desc)\"," *
                "$(fstats.max_flux),$(fstats.max_energy)," *
                "$(fstats.n_active),$(@sprintf("%.3f", frac))")
        end
    end
end

function write_text(rows, path, nova, run_dir, n_snapshots)
    src_counts = Dict{String,Int}()
    for (meta, _) in rows; src_counts[meta.source] = get(src_counts, meta.source, 0) + 1; end

    open(path, "w") do f
        println(f, "Active reaction flux summary")
        println(f, "Nova case  : $(nova)")
        println(f, "Run dir    : $(run_dir)")
        println(f, "Reactions  : $(length(rows)) with non-trivial flux (|dY/dt| > $(FLUX_SENTINEL))")
        println(f, "Snapshots  : $(n_snapshots) flux files")
        println(f, "Ranked by peak |dY/dt| across all trajectory snapshots")
        println(f, "=" ^ 105)
        println(f)
        @printf(f, "%5s  %-38s  %-7s  %9s  %-6s  %12s  %12s  %8s\n",
            "#", "Reaction", "Type", "Q (MeV)", "Source", "peak dY/dt", "peak erg/g/s", "active/N")
        println(f, "-" ^ 105)
        for (meta, fstats) in rows
            frac = "$(fstats.n_active)/$(n_snapshots)"
            @printf(f, "%5d  %-38s  %-7s  %9.4f  %-6s  %12.4e  %12.4e  %8s\n",
                meta.index, meta.notation, meta.rtype, meta.Q_MeV, meta.source,
                fstats.max_flux, fstats.max_energy, frac)
        end
        println(f)
        println(f, "=" ^ 105)
        println(f, "Rate source legend")
        println(f, "-" ^ 105)
        for (src, n) in sort(collect(src_counts); by = x -> -x[2])
            desc = get(SOURCE_NAMES, src, "(unknown)")
            @printf(f, "  %-6s  (%4d reactions)  %s\n", src, n, desc)
        end
    end
end

# ── CLI ────────────────────────────────────────────────────────────────────────

function usage()
    println("""
Usage:
  julia tools/list_flux_reactions.jl [options]

Options:
  --nova NAME        Nova case directory (required)
  --run  DIR         Run directory (default: runs)
  --out  STEM        Output file stem without extension
                     (default: <nova>/results/flux_reaction_list_<run>)
  --threshold FLUX   Min |dY/dt| to include a reaction (default: 1e-60)
  --include-vital    Include VITAL block reactions (index ≤ 117, omitted by default)
  --all              Include all active reactions regardless of flux (implies no threshold)
  -h, --help         Show this help

Output:
  <stem>.csv    CSV ranked by peak flux, with n_active_snapshots and fraction_active
  <stem>.txt    Human-readable table + source legend

Reactions are sorted by peak |dY/dt| across all flux_*.DAT snapshots, so the
most physically important reactions appear first.
""")
end

function main(args)
    nova      = ""
    run_dir   = "runs"
    out_stem  = ""
    threshold = FLUX_SENTINEL
    incvital  = false
    all_rxns  = false

    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-h", "--help")
            usage(); return
        elseif a == "--nova"
            i += 1; nova = args[i]
        elseif a == "--run"
            i += 1; run_dir = args[i]
        elseif a == "--out"
            i += 1; out_stem = args[i]
        elseif a == "--threshold"
            i += 1; threshold = parse(Float64, args[i])
        elseif a == "--include-vital"
            incvital = true
        elseif a == "--all"
            all_rxns = true; threshold = 0.0
        else
            error("Unknown argument: $a")
        end
        i += 1
    end

    isempty(nova) && error("--nova is required")

    nova_dir = joinpath(NOVA_CASES, nova)
    isdir(nova_dir) || error("Nova case not found: $nova_dir")

    baseline = joinpath(nova_dir, run_dir, "baseline")
    isdir(baseline) || error("Baseline dir not found: $baseline")

    nwsetup = joinpath(baseline, "networksetup.txt")
    isfile(nwsetup) || error("networksetup.txt not found: $nwsetup")

    println("Parsing networksetup from $nwsetup ...")
    meta_by_idx = parse_networksetup(nwsetup)

    println("Reading flux files from $baseline ...")
    flux_stats = read_flux_files(baseline)
    n_snapshots = length(filter(f -> startswith(f, "flux_") && endswith(f, ".DAT"),
                                readdir(baseline)))
    println("  $(n_snapshots) snapshot files, $(length(flux_stats)) reactions with flux data")

    # Build joined list: only active reactions with flux above threshold
    results = Tuple{RxnMeta, FluxStats}[]
    for (idx, fstats) in flux_stats
        meta = get(meta_by_idx, idx, nothing)
        meta === nothing && continue
        meta.active || continue
        incvital || idx > 117 || continue
        (all_rxns || fstats.max_flux > threshold) || continue
        push!(results, (meta, fstats))
    end

    # Sort by peak flux descending
    sort!(results; by = x -> -x[2].max_flux)

    results_dir = joinpath(nova_dir, "results")
    mkpath(results_dir)
    stem = isempty(out_stem) ? joinpath(results_dir, "flux_reaction_list_$(run_dir)") : out_stem

    write_csv(results, stem * ".csv", n_snapshots)
    write_text(results, stem * ".txt", nova, run_dir, n_snapshots)

    println("  $(length(results)) reactions with |dY/dt| > $(threshold)")
    println("  CSV : $(stem).csv")
    println("  Text: $(stem).txt")

    src_counts = Dict{String,Int}()
    for (meta, _) in results; src_counts[meta.source] = get(src_counts, meta.source, 0) + 1; end
    println("\nSource breakdown (flux-active reactions):")
    for (src, n) in sort(collect(src_counts); by = x -> -x[2])
        @printf("  %-6s  %4d  %s\n", src, n, get(SOURCE_NAMES, src, ""))
    end
end

main(ARGS)
