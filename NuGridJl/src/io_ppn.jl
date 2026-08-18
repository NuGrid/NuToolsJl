# io_ppn.jl — readers for nuppn per-cycle text output.
#
# Covers the abundance vector (`iso_massfNNNNN.DAT` and the decay variant),
# reaction fluxes (`flux_NNNNN.DAT`), the `x-time.dat` time series, and the
# `iniab_*.dat` initial-abundance / `trajectory.input` input files.  Every reader
# returns a typed value, not a bag of parallel arrays.

"Parse a Fortran real, tolerating `D`/`d` exponents (`\"1.00D-03\"` → `1.0e-3`)."
function parse_fortran(s::AbstractString)
    return parse(Float64, replace(strip(String(s)), r"[DdEe]" => "E"))
end

# ---------------------------------------------------------------------------
# Abundances (iso_massf / iso_massfdecay)
# ---------------------------------------------------------------------------

"""
    Abundances

Mass fractions for one cycle: `cycle`, `time` (s), `t9` (GK) and `rho`
(g cm⁻³) from the file header, plus parallel `iso::Vector{Isotope}` and
`X::Vector{Float64}`.  Index by isotope — `ab[Isotope("26Al")]` or `ab["c12"]`
— to get a mass fraction (0.0 if the species is absent).
"""
struct Abundances
    cycle::Int
    time::Float64
    t9::Float64
    rho::Float64
    iso::Vector{Isotope}
    X::Vector{Float64}
end

function _match_float(line, key)
    m = match(Regex(key * raw"\s*=?\s*([-+]?[\d.]+(?:[DdEe][-+]?\d+)?)"), line)
    m === nothing ? nothing : parse_fortran(m.captures[1])
end

"""
    read_abundances(path) -> Abundances

Read one `iso_massf*.DAT` file (the decay file has an identical layout).  Uses the
numeric Z/A/ISOM columns, so isotope identity never depends on re-parsing the
trailing text label.
"""
function read_abundances(path::AbstractString)
    cycle = 0
    time = NaN; t9 = NaN; rho = NaN
    isos = Isotope[]
    X = Float64[]

    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue

        if startswith(line, "#")
            if occursin("mod", line)
                m = match(r"mod\s+(-?\d+)", line)
                m !== nothing && (cycle = parse(Int, m.captures[1]))
                a = _match_float(line, "agej")
                a !== nothing && (time = a)
            elseif occursin("t9", line)
                v = _match_float(line, "t9"); v !== nothing && (t9 = v)
                v = _match_float(line, "rho"); v !== nothing && (rho = v)
            end
            continue
        end
        startswith(line, "H NUM") && continue

        parts = split(line)
        length(parts) < 6 && continue
        z = tryparse(Float64, parts[2]); z === nothing && continue
        a = tryparse(Float64, parts[3]); a === nothing && continue
        isom = tryparse(Float64, parts[4]); isom === nothing && continue
        x = tryparse(Float64, parts[5]); x === nothing && continue

        push!(isos, Isotope(Int(round(z)), Int(round(a)), max(Int(round(isom)) - 1, 0)))
        push!(X, x)
    end

    isempty(isos) && throw(ArgumentError("no abundance rows parsed from $path"))
    return Abundances(cycle, time, t9, rho, isos, X)
end

Base.length(ab::Abundances) = length(ab.iso)
isotopes(ab::Abundances) = ab.iso

function Base.getindex(ab::Abundances, iso::Isotope)
    idx = findfirst(==(iso), ab.iso)
    idx === nothing ? 0.0 : ab.X[idx]
end
Base.getindex(ab::Abundances, s::AbstractString) = ab[Isotope(s)]

"Mass fraction of `iso` in `ab` (0.0 if absent)."
mass_fraction(ab::Abundances, iso) = ab[iso isa Isotope ? iso : Isotope(iso)]

"""
    DataFrame(ab::Abundances) -> DataFrame

Tidy table with `:Z, :A, :N, :isomer, :isotope, :X` — the form the nuclear-chart
plots consume.
"""
function DataFrames.DataFrame(ab::Abundances)
    DataFrame(
        Z = [i.Z for i in ab.iso],
        A = [i.A for i in ab.iso],
        N = [neutron_number(i) for i in ab.iso],
        isomer = [i.isomer for i in ab.iso],
        isotope = [isotope_name(i) for i in ab.iso],
        X = ab.X,
    )
end

# ---------------------------------------------------------------------------
# Fluxes (flux_NNNNN.DAT)
# ---------------------------------------------------------------------------

"""
    read_fluxes(path) -> DataFrame

Reaction fluxes for one cycle.  Columns: `:index`, the start/end nuclide
coordinates (`:z_start, :a_start, :n_start, :z_end, :a_end, :n_end`), the start
and end `Isotope`s, and `:flux` (dY/dt), `:energy`, `:timescale`.  The flow runs
from the first reactant (k1) to the last product (k7), matching the file layout.
"""
function read_fluxes(path::AbstractString)
    rows = NamedTuple[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        p = split(line)
        length(p) < 12 && continue
        vals = tryparse.(Float64, p[1:9])
        any(isnothing, vals) && continue
        idx = Int(round(vals[1]))
        z1, a1 = Int(round(vals[2])), Int(round(vals[3]))
        z7, a7 = Int(round(vals[8])), Int(round(vals[9]))
        flux = tryparse(Float64, p[10]); flux === nothing && continue
        energy = tryparse(Float64, p[11])
        timescale = tryparse(Float64, p[12])
        push!(rows, (
            index = idx,
            z_start = z1, a_start = a1, n_start = a1 - z1,
            z_end = z7, a_end = a7, n_end = a7 - z7,
            reactant = Isotope(z1, a1), product = Isotope(z7, a7),
            flux = flux,
            energy = energy === nothing ? NaN : energy,
            timescale = timescale === nothing ? NaN : timescale,
        ))
    end
    isempty(rows) && throw(ArgumentError("no flux rows parsed from $path"))
    return DataFrame(rows)
end

# ---------------------------------------------------------------------------
# Time series (x-time.dat)
# ---------------------------------------------------------------------------

"""
    XTime

Whole-run time series from `x-time.dat`: parallel `cycle, time, t9, rho` vectors,
the species list `iso`, and a `nrows × nspecies` mass-fraction matrix `X`.  Index
by isotope (`xt[Isotope("c12")]`) for that species' column, or call
[`series`](@ref) to get an `(x, y)` pair against time, T9 or ρ.
"""
struct XTime
    cycle::Vector{Int}
    time::Vector{Float64}
    t9::Vector{Float64}
    rho::Vector{Float64}
    iso::Vector{Isotope}
    X::Matrix{Float64}
    _index::Dict{Isotope,Int}
end

function read_xtime(path::AbstractString)
    header = ""
    for raw in eachline(path)
        s = strip(raw)
        if startswith(s, "#")
            header = s
            break
        end
    end
    isempty(header) && throw(ArgumentError("no header line found in $path"))

    fields = strip.(split(lstrip(header, ['#']), '|'))
    fields = filter(!isempty, fields)
    length(fields) < 7 && throw(ArgumentError("unexpected x-time header in $path"))

    species = Isotope[]
    for f in fields[7:end]
        dash = findfirst('-', f)
        token = dash === nothing ? f : f[nextind(f, dash):end]
        iso = tryparse(Isotope, strip(token))
        iso === nothing || push!(species, iso)
    end

    data = readdlm(path; comments = true, comment_char = '#')
    ncols = size(data, 2)
    nspecies = min(length(species), ncols - 6)
    species = species[1:nspecies]

    Xmat = Float64.(data[:, 7:(6 + nspecies)])
    idx = Dict(iso => j for (j, iso) in enumerate(species))
    return XTime(
        Int.(round.(data[:, 1])),
        Float64.(data[:, 2]),
        Float64.(data[:, 3]),
        Float64.(data[:, 4]),
        species, Xmat, idx,
    )
end

Base.length(xt::XTime) = length(xt.cycle)

function Base.getindex(xt::XTime, iso::Isotope)
    j = get(xt._index, iso, nothing)
    j === nothing ? fill(0.0, length(xt.cycle)) : xt.X[:, j]
end
Base.getindex(xt::XTime, s::AbstractString) = xt[Isotope(s)]

"""
    series(xt, iso; x = :time) -> (xvals, yvals)

Mass fraction of `iso` versus `:time` (s), `:t9` (GK) or `:rho` (g cm⁻³).
"""
function series(xt::XTime, iso; x::Symbol = :time)
    isotope = iso isa Isotope ? iso : Isotope(iso)
    xvals = x === :time ? xt.time : x === :t9 ? xt.t9 : x === :rho ? xt.rho :
        throw(ArgumentError("x must be :time, :t9 or :rho (got :$x)"))
    return xvals, xt[isotope]
end

# ---------------------------------------------------------------------------
# Input files (iniab, trajectory)
# ---------------------------------------------------------------------------

"""
    read_iniab(path) -> Abundances

Initial abundances from an `iniab_*.dat` file (`Z SYMBOL A value`, or
`Z PROT value`).  Returned as `Abundances` with sentinel `cycle = -1` and `NaN`
thermodynamics so it drops straight into the abundance charts.
"""
function read_iniab(path::AbstractString)
    isos = Isotope[]
    X = Float64[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        p = split(line)
        length(p) < 3 && continue
        tryparse(Int, p[1]) === nothing && continue
        up = uppercase(p[2])
        if up == "PROT"
            push!(isos, Isotope(1, 1, 0)); push!(X, parse_fortran(p[3]))
        elseif up == "NEUT"
            push!(isos, Isotope(0, 1, 0)); push!(X, parse_fortran(p[3]))
        elseif length(p) >= 4
            z = parse(Int, p[1]); a = parse(Int, p[3])
            push!(isos, Isotope(z, a, 0)); push!(X, parse_fortran(p[4]))
        end
    end
    isempty(isos) && throw(ArgumentError("no abundances parsed from $path"))
    return Abundances(-1, NaN, NaN, NaN, isos, X)
end

"""
    read_trajectory(path) -> DataFrame

Thermodynamic trajectory (`time_s`, `temperature_T9`, `density_cgs`), applying the
`AGEUNIT`/`TUNIT`/`RHOUNIT` conversions declared in the file header.
"""
function read_trajectory(path::AbstractString)
    time = Float64[]; temp = Float64[]; dens = Float64[]
    ageunit = "SEC"; tunit = "T9K"; rhounit = "CGS"

    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        if occursin('=', line)
            key, value = strip.(split(line, '=', limit = 2))
            key = uppercase(key); value = uppercase(value)
            key == "AGEUNIT" && (ageunit = value; continue)
            key == "TUNIT" && (tunit = value; continue)
            key == "RHOUNIT" && (rhounit = value; continue)
            continue
        end
        p = split(line)
        length(p) < 3 && continue
        vals = tryparse.(Float64, p[1:3])
        any(isnothing, vals) && continue
        t, T, rho = vals

        if ageunit in ("YRS", "YR", "YEAR", "YEARS")
            t *= 365.25 * 24 * 3600
        elseif !(ageunit in ("SEC", "S"))
            throw(DomainError(ageunit, "unsupported trajectory AGEUNIT"))
        end
        if tunit == "T8K"
            T /= 10
        elseif tunit != "T9K"
            throw(DomainError(tunit, "unsupported trajectory TUNIT (expected T9K or T8K)"))
        end
        if rhounit == "LOG"
            rho = 10.0^rho
        elseif rhounit != "CGS"
            throw(DomainError(rhounit, "unsupported trajectory RHOUNIT (expected CGS or LOG)"))
        end

        push!(time, t); push!(temp, T); push!(dens, rho)
    end
    return DataFrame(time_s = time, temperature_T9 = temp, density_cgs = dens)
end
