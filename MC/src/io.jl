# io.jl: reading the ppn files a Monte Carlo needs (x-time.dat, rate_uncertainty_x.txt, iniab, namelists).

"""
    fnum(s) -> Float64

Parse a Fortran number, including the form with the `E` dropped from a 3-digit exponent
(`"2.24366-100"`). Unparseable strings give `NaN`.
"""
function fnum(s::AbstractString)
    v = tryparse(Float64, s)
    v !== nothing && return v
    m = match(r"^([-+]?[0-9.]+)([-+][0-9]{3})$", s)
    m === nothing ? NaN : parse(Float64, m[1] * "E" * m[2])
end

"""
    namelist_value(text, key) -> Union{String,Nothing}

The quoted string value of `key = '...'` in a Fortran namelist text, or `nothing`.
"""
namelist_value(text::AbstractString, key::AbstractString) =
    (m = match(Regex("$key\\s*=\\s*'([^']+)'"), text); m === nothing ? nothing : String(m[1]))

"""
    iso_name(column) -> Union{String,Nothing}

Isotope name of an x-time.dat column: `"38-AL 26"` → `"Al26"`, `"1103-AL26M"` → `"Al26m"`,
`"2-PROT"` → `"H1"`; non-isotope columns (`cycle`, `time`, ...) give `nothing`.
"""
function iso_name(col::AbstractString)
    m = match(r"^\s*\d+-([A-Z]+)\s*(\d+)([A-Z]?)\s*$", col)
    m !== nothing && return uppercasefirst(lowercase(m[1])) * string(parse(Int, m[2])) * lowercase(m[3])
    occursin(r"^\s*\d+-PROT", col) && return "H1"
    nothing
end

"""
    read_iniab(path) -> Dict{String,Float64}

Initial mass fractions from a ppn iniab file (either `1 PROT 7.0E-01` / `6 C 12 ...` or `6 C12 ...`
layout), keyed by isotope name (`"H1"`, `"C12"`, ...).
"""
function read_iniab(path::AbstractString)
    out = Dict{String,Float64}()
    for l in eachline(path)
        t = split(l)
        (length(t) < 3 || !occursin(r"^\d+$", t[1])) && continue
        if uppercase(t[2]) == "PROT"
            out["H1"] = fnum(t[3]); continue
        end
        m = match(r"^([A-Za-z]+)(\d+)$", t[2])
        el, a, x = m !== nothing ? (m[1], m[2], t[3]) : (t[2], t[3], t[4])
        out[uppercasefirst(lowercase(el)) * string(parse(Int, a))] = fnum(x)
    end
    out
end

"""
    run_iniab(run_dir) -> Dict{String,Float64}

The initial abundances of a ppn run directory, from the `ini_filename` in its `ppn_frame.input`.
"""
function run_iniab(run_dir::AbstractString)
    ini = namelist_value(read(joinpath(run_dir, "ppn_frame.input"), String), "ini_filename")
    ini === nothing && error("$run_dir/ppn_frame.input: no ini_filename")
    read_iniab(joinpath(run_dir, ini))
end

"""
    last_xtime_row(path) -> (names, values)

Column names and the raw value strings of the last data row of an x-time.dat, read from the end of the
file (x-time.dat can be large).
"""
function last_xtime_row(path::AbstractString)
    open(path) do io
        hdr = readline(io)
        sz = filesize(path)
        back = min(sz, 1 << 16)
        local tail
        while true
            seek(io, sz - back)
            tail = split(rstrip(read(io, String)), '\n')
            (length(tail) >= 2 || back == sz) && break
            back = min(sz, 2back)
        end
        String.(strip.(split(hdr, '|')[2:end])), String.(split(tail[end]))
    end
end

"""
    final_abundances(run_dir) -> Dict{String,Float64}

Mass fractions of the last x-time.dat row of a finished ppn run, keyed by isotope name. Isomer species
keep their own key (`"Al26m"`).
"""
function final_abundances(run_dir::AbstractString)
    names, vals = last_xtime_row(joinpath(run_dir, "x-time.dat"))
    out = Dict{String,Float64}()
    for (c, v) in zip(names, vals)
        nm = iso_name(c)
        nm === nothing || (out[nm] = fnum(v))
    end
    out
end

"""
    read_rate_x(path) -> Vector{NamedTuple}

The rows of a ppn `rate_uncertainty_x.txt`: `(row, label, reaction, x, how)` for every varied rate,
in file order.
"""
function read_rate_x(path::AbstractString)
    out = NamedTuple{(:row, :label, :reaction, :x, :how),Tuple{Int,String,String,Float64,String}}[]
    for l in eachline(path)
        (startswith(l, '#') || isempty(strip(l))) && continue
        m = match(r"^\s*(\d+)\s+(\S{5})\s+(.{46})\s*(\S+)\s*(.*)$", l)
        m === nothing && continue
        push!(out, (row=parse(Int, m[1]), label=String(m[2]), reaction=join(split(m[3]), ' '),
                    x=parse(Float64, m[4]), how=String(m[5])))
    end
    out
end
