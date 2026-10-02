#!/usr/bin/env julia
# list_network_rates.jl
#
# List all active reactions and their rate sources for a nova case.
# Julia equivalent of list_network_rates.py.
#
# Usage:
#   julia tools/list_network_rates.jl --nova ne_nova_1.15_12_X_weiss_mixed
#   julia tools/list_network_rates.jl --nova ne_nova_1.15_12_X_weiss_mixed --run runs_star_opt2_novital
#   julia tools/list_network_rates.jl --help

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

# ── Species formatting ─────────────────────────────────────────────────────────

# sp: 5-char string from networksetup ("NEUT ", "PROT ", "C  12", "OOOOO", ...)
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

# ── Reaction notation ──────────────────────────────────────────────────────────

function build_notation(n1, sp1, n2, sp2, n3, sp3, n4, sp4, rtype)
    ch = get(RTYPE_CHANNEL, rtype, nothing)

    if ch === nothing
        # compound / (v,v): arrow notation
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
        residual = fmt_sp(sp3)
        isempty(residual) && (residual = fmt_sp(sp4))
        return "$(target)(γ,$(eject))$(residual)"
    end

    if proj in ("β+", "β-", "e-")
        daughter = fmt_sp(sp3)
        return "$(target)($(proj),$(eject))$(daughter)"
    end

    residual = fmt_sp(sp3)
    isempty(residual) && (residual = fmt_sp(sp4))
    return "$(target)($(proj),$(eject))$(residual)"
end

# ── networksetup.txt parsing ───────────────────────────────────────────────────

# Fortran may omit 'E' when exponent has 3 digits, e.g. "1.774-105"
function parse_fortran_float(s::AbstractString)
    v = tryparse(Float64, s)
    v !== nothing && return v
    m = match(r"^([+-]?\d+\.?\d*)([+-]\d+)$", s)
    m !== nothing && return parse(Float64, "$(m.captures[1])E$(m.captures[2])")
    return 0.0
end

const NETWORK_RE = r"^\s*(\d+)\s+([TF])\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+->\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)"

struct RxnRow
    index     :: Int
    active    :: Bool
    n1 :: String; sp1 :: String
    n2 :: String; sp2 :: String
    n3 :: String; sp3 :: String
    n4 :: String; sp4 :: String
    Q_MeV     :: Float64
    source    :: String
    rtype     :: String
    notation  :: String
end

function parse_networksetup(path::AbstractString)
    rows = RxnRow[]
    for line in eachline(path)
        m = match(NETWORK_RE, line)
        m === nothing && continue
        idx, act, n1, sp1, n2, sp2, n3, sp3, n4, sp4, q_str, src, rtype, _ = m.captures
        notation = build_notation(n1, sp1, n2, sp2, n3, sp3, n4, sp4, rtype)
        push!(rows, RxnRow(
            parse(Int, idx),
            act == "T",
            n1, strip(sp1), n2, strip(sp2),
            n3, strip(sp3), n4, strip(sp4),
            parse_fortran_float(q_str),
            src, rtype, notation,
        ))
    end
    return rows
end

# ── Output ─────────────────────────────────────────────────────────────────────

function write_csv(rows, path)
    open(path, "w") do f
        println(f, "index,reaction,rtype,Q_MeV,source,source_description")
        for r in rows
            desc = get(SOURCE_NAMES, r.source, r.source)
            println(f, "$(r.index),\"$(r.notation)\",$(r.rtype),$(r.Q_MeV),$(r.source),\"$(desc)\"")
        end
    end
end

function write_text(rows, path, nova, run_dir)
    src_counts = Dict{String,Int}()
    for r in rows; src_counts[r.source] = get(src_counts, r.source, 0) + 1; end

    open(path, "w") do f
        println(f, "Network rate source list")
        println(f, "Nova case : $(nova)")
        println(f, "Run dir   : $(run_dir)")
        println(f, "Reactions : $(length(rows)) active")
        println(f, "=" ^ 90)
        println(f)
        @printf(f, "%5s  %-38s  %-7s  %9s  %s\n", "#", "Reaction", "Type", "Q (MeV)", "Source")
        println(f, "-" ^ 90)
        for r in rows
            @printf(f, "%5d  %-38s  %-7s  %9.4f  %s\n",
                r.index, r.notation, r.rtype, r.Q_MeV, r.source)
        end
        println(f)
        println(f, "=" ^ 90)
        println(f, "Rate source legend")
        println(f, "-" ^ 90)
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
  julia tools/list_network_rates.jl [options]

Options:
  --nova NAME        Nova case directory (required)
  --run  DIR         Run directory (default: runs)
  --out  STEM        Output file stem without extension (default: <nova>/results/network_rate_list_<run>)
  --include-vital    Include VITAL block reactions (index ≤ 117, omitted by default)
  -h, --help         Show this help

Output:
  <stem>.csv    CSV: index, reaction, rtype, Q_MeV, source, source_description
  <stem>.txt    Human-readable table + source legend
""")
end

function main(args)
    nova     = ""
    run_dir  = "runs"
    out_stem = ""
    incvital = false

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
        elseif a == "--include-vital"
            incvital = true
        else
            error("Unknown argument: $a")
        end
        i += 1
    end

    isempty(nova) && error("--nova is required")

    nova_dir  = joinpath(NOVA_CASES, nova)
    isdir(nova_dir) || error("Nova case not found: $nova_dir")

    nwsetup = joinpath(nova_dir, run_dir, "baseline", "networksetup.txt")
    isfile(nwsetup) || error("networksetup.txt not found: $nwsetup")

    println("Parsing $nwsetup ...")
    all_rows = parse_networksetup(nwsetup)

    rows = filter(r -> r.active, all_rows)
    incvital || filter!(r -> r.index > 117, rows)

    results_dir = joinpath(nova_dir, "results")
    mkpath(results_dir)
    stem = isempty(out_stem) ? joinpath(results_dir, "network_rate_list_$(run_dir)") : out_stem

    write_csv(rows, stem * ".csv")
    write_text(rows, stem * ".txt", nova, run_dir)

    println("  $(length(rows)) active reactions written")
    println("  CSV : $(stem).csv")
    println("  Text: $(stem).txt")

    src_counts = Dict{String,Int}()
    for r in rows; src_counts[r.source] = get(src_counts, r.source, 0) + 1; end
    println("\nSource breakdown:")
    for (src, n) in sort(collect(src_counts); by = x -> -x[2])
        @printf("  %-6s  %4d  %s\n", src, n, get(SOURCE_NAMES, src, ""))
    end
end

main(ARGS)
