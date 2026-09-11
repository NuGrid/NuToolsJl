# starlib_curve.jl — reading STARLIB's own tabulated (T9, rate, factor)
# reaction curves directly. Independent of REACLIB's exp(...)-formula
# rate_curve/rate_plot (npdata.jl): STARLIB reactions are tabulated
# point-by-point (60 T9 values per reaction), not a fitted analytic formula,
# so they need their own reader. Column offsets/block structure verified
# directly against starlib_mc10_mc13_082022.txt and ETR25_starlib.dat this
# session (see NuGrid-Tools/NuGridJl/tools/build_sigma_sweep.jl, which
# locates/rewrites these same reaction blocks for the σ-sweep tool).

const _STARLIB_CURVE_NT9 = 60
const _STARLIB_CURVE_HEADER_LINES = 1

"""
    starlib_rate_curve(path, target_species) -> DataFrame

Read one reaction's tabulated rate curve directly from a STARLIB reaction
file (`starlib_mc10_mc13_082022.txt`/`ETR25_starlib.dat`-format: one header
line, then repeating blocks of one fixed-width reaction-info line followed by
60 free-format "T9 rate factor" lines). `target_species` is a set of species
tokens exactly as they appear in the file (e.g. `["p", "f18", "he4", "o15"]`
for 18F(p,α)15O) — matched the same way `build_sigma_sweep.jl`'s σ-sweep tool
locates a reaction, and must resolve to exactly one row (throws otherwise).
Columns: `:T9` (GK), `:rate`, `:factor` (STARLIB's per-T9 lognormal
uncertainty spread — see [`Reaction`](@ref)/`build_sigma_sweep.jl`'s
docstrings for what it means).
"""
function starlib_rate_curve(path::AbstractString, target_species::Vector{<:AbstractString})
    lines = readlines(path)
    target = Set(target_species)
    idx = _STARLIB_CURVE_HEADER_LINES + 1
    n = length(lines)
    matches = Int[]
    while idx <= n
        line = lines[idx]
        species = String[]
        for k in 0:5
            a, b = 6 + k * 5, 10 + k * 5
            tok = length(line) >= b ? strip(line[a:b]) : ""
            isempty(tok) || push!(species, tok)
        end
        Set(species) == target && push!(matches, idx)
        idx += _STARLIB_CURVE_NT9 + 1
    end
    length(matches) == 1 || throw(ArgumentError(
        "expected exactly 1 reaction matching $target_species in $path, found $(length(matches))"))
    info_line = only(matches)

    T9 = Float64[]; rate = Float64[]; factor = Float64[]
    for j in 1:_STARLIB_CURVE_NT9
        line_idx = info_line + j
        parts = split(lines[line_idx])
        length(parts) == 3 || throw(ArgumentError("unexpected T9 data line at $path:$line_idx: $(lines[line_idx])"))
        t9, r0, fact = parse.(Float64, parts)
        push!(T9, t9); push!(rate, r0); push!(factor, fact)
    end
    return DataFrame(T9 = T9, rate = rate, factor = factor)
end
