# npdata.jl — REACLIB reaction-rate curves from NPDATA.
#
# Files under NPDATA/REACLIB/ are the standard JINA REACLIB ASCII format: a
# bare chapter number (1-8) followed by 3-line blocks (a species/label/Q-value
# line, then two lines of the 7 rate-formula coefficients). Multiple REACLIB
# vintages commonly live side by side in one NPDATA snapshot with no reliable
# way to guess which file a given networksetup.txt source label came from, so
# `rate_curve` searches every file under `npdata_root/REACLIB` and returns
# every match rather than silently picking one — compare the returned
# `:label`/`:file` columns against `Reaction.source` to see which curve is
# actually in use.

const _REACLIB_CHAPTER_COUNTS = Dict(1 => (1, 1), 2 => (1, 2), 3 => (1, 3), 4 => (2, 1),
                                      5 => (2, 2), 6 => (2, 3), 7 => (2, 4), 8 => (3, 1))

"REACLIB species token for `iso`, e.g. `Isotope(2, 4, 0) -> \"he4\"`; protons/neutrons become `\"p\"`/`\"n\"`."
function reaclib_species(iso::Isotope)
    iso == Isotope(1, 1, 0) && return "p"
    iso == Isotope(0, 1, 0) && return "n"
    return lowercase(string(element_symbol(iso.Z), iso.A))
end

_parse_reaclib_coefficients(line::AbstractString) = [
    parse(Float64, m.match)
    for m in eachmatch(r"[+-]?(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?", replace(String(line), 'D' => 'e', 'd' => 'e'))
]

function _parse_reaclib_species_line(chapter::Integer, line::AbstractString)
    tokens = split(line)
    length(tokens) < 4 && return nothing
    label = String(tokens[end - 1])
    qvalue = tryparse(Float64, tokens[end])
    qvalue === nothing && return nothing

    counts = get(_REACLIB_CHAPTER_COUNTS, chapter, nothing)
    counts === nothing && return nothing
    nreactants, nproducts = counts
    species = String.(tokens[1:(end - 2)])
    length(species) != nreactants + nproducts && return nothing

    return (reactants = sort(species[1:nreactants]), products = sort(species[(nreactants + 1):end]),
            label = label, qvalue = qvalue)
end

"One REACLIB rate block: the reaction it belongs to, its label/Q-value, and the 7 rate-formula coefficients."
struct ReaclibBlock
    file::String
    chapter::Int
    reactants::Vector{String}
    products::Vector{String}
    label::String
    qvalue::Float64
    coefficients::NTuple{7,Float64}
end

"""
    read_reaclib_blocks(path) -> Vector{ReaclibBlock}

Parse every rate block from a REACLIB-format file.
"""
function read_reaclib_blocks(path::AbstractString)
    lines = [strip(l) for l in readlines(path) if !isempty(strip(l))]
    blocks = ReaclibBlock[]
    chapter = nothing
    i = 1
    while i <= length(lines)
        line = lines[i]
        if occursin(r"^[1-8]$", line)
            chapter = parse(Int, line)
            i += 1
            continue
        end
        if chapter === nothing || i + 2 > length(lines)
            i += 1
            continue
        end
        parsed = _parse_reaclib_species_line(chapter, lines[i])
        coeffs = vcat(_parse_reaclib_coefficients(lines[i + 1]), _parse_reaclib_coefficients(lines[i + 2]))
        if parsed !== nothing && length(coeffs) == 7
            push!(blocks, ReaclibBlock(path, chapter, parsed.reactants, parsed.products,
                                         parsed.label, parsed.qvalue, Tuple(coeffs)))
            i += 3
        else
            i += 1
        end
    end
    return blocks
end

"""
    reaclib_rate(coefficients, T9) -> Float64

The REACLIB rate formula:
`exp(a1 + a2/T9 + a3/T9^(1/3) + a4*T9^(1/3) + a5*T9 + a6*T9^(5/3) + a7*ln(T9))`.
"""
function reaclib_rate(coefficients, T9::Real)
    T9 <= 0 && throw(DomainError(T9, "temperature must be positive and in GK"))
    T13 = T9^(1 / 3)
    a = coefficients
    exponent = a[1] + a[2] / T9 + a[3] / T13 + a[4] * T13 + a[5] * T9 + a[6] * T9 * T13 * T13 + a[7] * log(T9)
    exponent > log(floatmax(Float64)) && return Inf
    exponent < log(floatmin(Float64)) && return 0.0
    return exp(exponent)
end

"A log-spaced temperature grid in GK spanning the astrophysically relevant range by default."
rate_curve_temperature_grid(; Tmin = 0.01, Tmax = 10.0, n = 300) = 10 .^ range(log10(Tmin), log10(Tmax), length = n)

# NPDATA/REACLIB directories mix actual REACLIB rate tables with unrelated
# files of the same vintage (winvn*.dat species tables, readme*, factors.txt,
# history.txt, sunet*.dat) that are a different format entirely and, at a few
# MB to tens of MB each, not worth scanning line-by-line on every lookup.
const _REACLIB_SKIP_NAME = r"readme|history|factors|^sunet|^winvn"i

"""
    rate_curve(r::Reaction, npdata_root; temperatures = rate_curve_temperature_grid(),
               files = nothing) -> DataFrame

Every REACLIB rate curve under `npdata_root/REACLIB` whose reactants/products
match `r`, evaluated over `temperatures` (GK). Columns: `:T9, :rate, :label,
:file, :chapter, :qvalue`. Multiple files commonly carry the same reaction
under different labels (different rate evaluations, e.g. different REACLIB
vintages) — all matches are returned rather than guessed at; cross-check
`:label`/`:file` against `r.source` (from [`read_network`](@ref)) to see which
one is actually wired into the network. By default every file in the REACLIB
directory is searched (skipping obviously non-REACLIB files by name); pass
`files` (a collection of basenames) to restrict the search — useful since the
full search can take tens of seconds against a large NPDATA snapshot.
"""
function rate_curve(r::Reaction, npdata_root::AbstractString; temperatures = rate_curve_temperature_grid(),
                     files = nothing)
    reaclib_dir = joinpath(npdata_root, "REACLIB")
    isdir(reaclib_dir) || throw(ArgumentError("no REACLIB directory under $npdata_root"))

    target_reactants = sort(reaclib_species.(r.reactants))
    target_products = sort(reaclib_species.(r.products))

    candidates = files === nothing ?
        [f for f in readdir(reaclib_dir) if !occursin(_REACLIB_SKIP_NAME, f)] : collect(files)

    rows = NamedTuple[]
    for name in candidates
        file = joinpath(reaclib_dir, name)
        isfile(file) || continue
        for block in read_reaclib_blocks(file)
            (block.reactants == target_reactants && block.products == target_products) || continue
            for T in temperatures
                push!(rows, (T9 = T, rate = reaclib_rate(block.coefficients, T), label = block.label,
                              file = basename(file), chapter = block.chapter, qvalue = block.qvalue))
            end
        end
    end
    isempty(rows) && throw(ArgumentError("no REACLIB rate data for $(label(r)) under $reaclib_dir"))
    return DataFrame(rows)
end

"""
    rate_plot(r::Reaction, npdata_root; temperatures = rate_curve_temperature_grid(),
              title = label(r), figure_size = (900, 600)) -> CM.Figure

Log-log reaction-rate-vs-temperature plot for every REACLIB curve
[`rate_curve`](@ref) finds for `r`, one line per `(label, file)` pair — the
"what rate was used, from where" companion to [`describe_rate`](@ref).
"""
function rate_plot(r::Reaction, npdata_root::AbstractString; temperatures = rate_curve_temperature_grid(),
                    title = label(r), figure_size = (900, 600))
    curves = rate_curve(r, npdata_root; temperatures)
    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1]; xlabel = "T9 (GK)", ylabel = "rate", title, xscale = log10, yscale = log10)
        for (i, group) in enumerate(groupby(curves, [:label, :file]))
            positive = group[group.rate .> 0, :]
            isempty(positive) && continue
            g = first(group)
            CM.lines!(ax, positive.T9, positive.rate; label = "$(g.label) ($(g.file))",
                       color = NUGRID_PALETTE[mod1(i, length(NUGRID_PALETTE))])
        end
        CM.axislegend(ax; position = :lt, framevisible = false)
        fig
    end
end
