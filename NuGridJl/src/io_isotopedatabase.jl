# io_isotopedatabase.jl — reader for isotopedatabase*.txt.
#
# Lists every isotope NuGrid knows how to add to a charged-particle network:
# proton/mass number, the mass number of that element's reference stable
# isotope, and whether the isotope is available (T) or excluded (F) from the
# current build. The element symbol column is redundant with Z and dropped —
# `Isotope` is built straight from the numeric Z/A columns, same as the other
# readers in this package.

"""
    IsotopeDatabaseEntry

One row of `isotopedatabase*.txt`: the isotope itself (`iso`), the mass number
of that element's reference stable isotope (`stable_a`), and whether it's
`active` (available to be added to a network).
"""
struct IsotopeDatabaseEntry
    iso::Isotope
    stable_a::Int
    active::Bool
end

"""
    read_isotopedatabase(path) -> Vector{IsotopeDatabaseEntry}

Parse an `isotopedatabase.txt` / `isotopedatabase_all.txt` / `isotopedatabase_cf.txt`
file: whitespace-delimited `Z  A  symbol  stable_A  T/F`, skipping blank and `#`
comment lines.
"""
function read_isotopedatabase(path::AbstractString)
    entries = IsotopeDatabaseEntry[]
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "#")) && continue
        p = split(line)
        length(p) < 5 && continue
        z = tryparse(Int, p[1]); z === nothing && continue
        a = tryparse(Int, p[2]); a === nothing && continue
        stable_a = tryparse(Int, p[4]); stable_a === nothing && continue
        flag = uppercase(p[5])
        (flag == "T" || flag == "F") || continue
        push!(entries, IsotopeDatabaseEntry(Isotope(z, a, 0), stable_a, flag == "T"))
    end
    isempty(entries) && throw(ArgumentError("no isotope database rows parsed from $path"))
    return entries
end
