using Plots
import CairoMakie as CM

CM.set_theme!(CM.theme_latexfonts())

const ELEMENT_SYMBOLS = [
    "n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",
]

const ELEMENT_NAMES = [
    "neutron", "hydrogen", "helium", "lithium", "beryllium", "boron",
    "carbon", "nitrogen", "oxygen", "fluorine", "neon", "sodium",
    "magnesium", "aluminum", "silicon", "phosphorus", "sulfur",
    "chlorine", "argon", "potassium", "calcium", "scandium", "titanium",
    "vanadium", "chromium", "manganese", "iron", "cobalt", "nickel",
    "copper", "zinc", "gallium", "germanium", "arsenic", "selenium",
    "bromine", "krypton", "rubidium", "strontium", "yttrium",
    "zirconium", "niobium", "molybdenum", "technetium", "ruthenium",
    "rhodium", "palladium", "silver", "cadmium", "indium", "tin",
    "antimony", "tellurium", "iodine", "xenon", "cesium", "barium",
    "lanthanum", "cerium", "praseodymium", "neodymium", "promethium",
    "samarium", "europium", "gadolinium", "terbium", "dysprosium",
    "holmium", "erbium", "thulium", "ytterbium", "lutetium", "hafnium",
    "tantalum", "tungsten", "rhenium", "osmium", "iridium", "platinum",
    "gold", "mercury", "thallium", "lead", "bismuth", "polonium",
    "astatine", "radon", "francium", "radium", "actinium", "thorium",
    "protactinium", "uranium",
]

const REACTION_STYLES = Dict(
    "18F_pa_15O" => (color = "#1f77b4", marker = :circle),
    "3He_ag_7Be" => (color = "#ff7f0e", marker = :rect),
    "7Be_pg_8B" => (color = "#2ca02c", marker = :diamond),
    "13N_pg_14O" => (color = "#d62728", marker = :utriangle),
    "13C_pg_14N" => (color = "#9467bd", marker = :dtriangle),
    "14N_pg_15O" => (color = "#8c564b", marker = :star5),
    "15N_pa_12C" => (color = "#e377c2", marker = :hexagon),
    "16O_pg_17F" => (color = "#7f7f7f", marker = :cross),
    "17O_pg_18F" => (color = "#bcbd22", marker = :xcross),
    "17O_pa_14N" => (color = "#17becf", marker = :pentagon),
    "18F_pg_19Ne" => (color = "#393b79", marker = :circle),
    "18O_pa_15N" => (color = "#637939", marker = :rect),
    "19F_pg_20Ne" => (color = "#8c6d31", marker = :diamond),
    "19F_pa_16O" => (color = "#843c39", marker = :utriangle),
    "20Ne_pg_21Na" => (color = "#7b4173", marker = :dtriangle),
    "21Ne_pg_22Na" => (color = "#3182bd", marker = :star5),
    "22Ne_pg_23Na" => (color = "#e6550d", marker = :hexagon),
    "21Na_pg_22Mg" => (color = "#31a354", marker = :cross),
    "22Na_pg_23Mg" => (color = "#756bb1", marker = :xcross),
    "23Na_pg_24Mg" => (color = "#636363", marker = :pentagon),
    "23Na_pa_20Ne" => (color = "#6baed6", marker = :circle),
    "25Mg_pg_26Al" => (color = "#fd8d3c", marker = :rect),
    "26Mg_pg_27Al" => (color = "#74c476", marker = :diamond),
    "26Al_pg_27Si" => (color = "#9e9ac8", marker = :utriangle),
    "27Al_pg_28Si" => (color = "#969696", marker = :dtriangle),
    "28Si_pg_29P" => (color = "#9ecae1", marker = :star5),
    "29Si_pg_30P" => (color = "#fdae6b", marker = :hexagon),
    "30Si_pg_31P" => (color = "#a1d99b", marker = :cross),
    "30P_pg_31S" => (color = "#bcbddc", marker = :xcross),
    "31P_pa_28Si" => (color = "#bdbdbd", marker = :pentagon),
    "32S_pg_33Cl" => (color = "#08519c", marker = :circle),
    "33S_pg_34Cl" => (color = "#a63603", marker = :rect),
    "34S_pg_35Cl" => (color = "#006d2c", marker = :diamond),
    "36Ar_pg_37K" => (color = "#54278f", marker = :utriangle),
    "37Ar_pg_38K" => (color = "#252525", marker = :dtriangle),
    "38Ar_pg_39K" => (color = "#6b6ecf", marker = :star5),
    "38K_pg_39Ca" => (color = "#bd9e39", marker = :hexagon),
    "39K_pg_40Ca" => (color = "#ad494a", marker = :cross),
)

function reaction_styles(reactions = nothing)
    reactions === nothing && return copy(REACTION_STYLES)

    styles = Dict()
    for reaction in sort(unique(collect(skipmissing(reactions))))
        styles[string(reaction)] = reaction_style(REACTION_STYLES, reaction)
    end
    return styles
end

function reaction_style(styles, reaction)
    return get(styles, string(reaction), (color = "#777777", marker = :circle))
end

const MARKER_SYMBOLS = Dict(
    :circle => "●",
    :rect => "■",
    :diamond => "◆",
    :utriangle => "▲",
    :dtriangle => "▼",
    :star5 => "★",
    :hexagon => "⬢",
    :cross => "+",
    :xcross => "×",
    :pentagon => "⬟",
)

struct ReactionStyleTable
    rows::Vector{NamedTuple{(:reaction, :color, :marker, :marker_symbol), Tuple{String, String, Symbol, String}}}
end

function reaction_style_table(reactions = nothing)
    selected_reactions = if reactions === nothing
        sort(collect(keys(REACTION_STYLES)))
    else
        sort(string.(collect(skipmissing(reactions))))
    end

    rows = [
        (
            reaction = reaction,
            color = string(reaction_style(REACTION_STYLES, reaction).color),
            marker = reaction_style(REACTION_STYLES, reaction).marker,
            marker_symbol = get(MARKER_SYMBOLS, reaction_style(REACTION_STYLES, reaction).marker, "●"),
        )
        for reaction in selected_reactions
    ]

    return ReactionStyleTable(rows)
end

function Base.show(io::IO, table::ReactionStyleTable)
    println(io, "Reaction style table")
    for row in table.rows
        println(io, row.reaction, "  ", row.color, "  ", row.marker)
    end
end

function Base.show(io::IO, ::MIME"text/html", table::ReactionStyleTable)
    print(io, """
    <table style="border-collapse: collapse; font-family: sans-serif; font-size: 13px;">
      <thead>
        <tr>
          <th style="text-align: left; padding: 4px 8px; border-bottom: 1px solid #999;">Reaction</th>
          <th style="text-align: left; padding: 4px 8px; border-bottom: 1px solid #999;">Color</th>
          <th style="text-align: left; padding: 4px 8px; border-bottom: 1px solid #999;">Marker</th>
        </tr>
      </thead>
      <tbody>
    """)

    for row in table.rows
        print(io, """
        <tr>
          <td style="padding: 3px 8px;">$(row.reaction)</td>
          <td style="padding: 3px 8px;">
            <span style="display: inline-block; width: 1em; height: 1em; background: $(row.color); border: 1px solid #333; vertical-align: -0.12em;"></span>
            <code>$(row.color)</code>
          </td>
          <td style="padding: 3px 8px;">
            <span style="color: $(row.color); font-size: 1.35em; font-weight: bold;">$(row.marker_symbol)</span>
            <code>:$(row.marker)</code>
          </td>
        </tr>
        """)
    end

    print(io, """
      </tbody>
    </table>
    """)
end

function scatter_reactions!(p, df, x, y; styles, ylabel_prefix = "", markerstrokecolor = :auto, show_legend = false)
    for df_reaction in groupby(df, :reaction)
        reaction = first(df_reaction.reaction)
        style = reaction_style(styles, reaction)
        rows = parentindices(df_reaction)[1]

        scatter!(
            p,
            x[rows],
            y[rows],
            label = show_legend ? string(ylabel_prefix, reaction) : "",
            color = style.color,
            markercolor = style.color,
            markershape = style.marker,
            markerstrokecolor = markerstrokecolor,
            markersize = 6,
        )
    end
end

function plot_trajectory(filepath; title = "Trajectory")
    trajectory = read_trajectory(filepath)
    return plot_trajectory(trajectory; title = title)
end

function plot_trajectory(trajectory::DataFrame; title = "Trajectory")
    if nrow(trajectory) == 0
        throw(ArgumentError("trajectory contains no data points"))
    end

    p = plot(
        trajectory.time_s,
        trajectory.temperature_T9,
        xlabel = "Time (s)",
        ylabel = "Temperature (T9)",
        title = title,
        label = "Temperature",
        color = :red,
        linewidth = 2,
        yscale = :log10,
        legend = :topright,
        xguidefontcolor = :black,
        xtickfontcolor = :black,
        yguidefontcolor = :red,
        ytickfontcolor = :red,
    )

    plot!(
        twinx(p),
        trajectory.time_s,
        trajectory.density_cgs,
        ylabel = "Density (g cm^-3)",
        label = "Density",
        color = :blue,
        linewidth = 2,
        yscale = :log10,
        xguidefontcolor = :black,
        xtickfontcolor = :black,
        yguidefontcolor = :blue,
        ytickfontcolor = :blue,
        legend = :bottomright,
    )

    return p
end

function plot_dens_temp(filepath; title = "Density vs Temperature")
    trajectory = read_trajectory(filepath)
    return plot_dens_temp(trajectory; title = title)
end

function plot_dens_temp(trajectory::DataFrame; title = "Density vs Temperature")
    if nrow(trajectory) == 0
        throw(ArgumentError("trajectory contains no data points"))
    end

    return plot(
        trajectory.temperature_T9,
        trajectory.density_cgs,
        xlabel = "Temperature (T9)",
        ylabel = "Density (g cm^-3)",
        title = title,
        label = "Density",
        color = :blue,
        linewidth = 2,
        xscale = :log10,
        yscale = :log10,
        xguidefontcolor = :black,
        xtickfontcolor = :black,
        yguidefontcolor = :black,
        ytickfontcolor = :black,
        legend = :topright,
    )
end

function plot_rate_curve(
    reaction_name;
    networksetup_path,
    npdata_path,
    temperatures = rate_curve_temperature_grid(),
    source = nothing,
    match_source_label = true,
    max_files = 4,
    title = nothing,
)
    curves = read_rate_curve_data(
        reaction_name;
        networksetup_path = networksetup_path,
        npdata_path = npdata_path,
        temperatures = temperatures,
        source = source,
        match_source_label = match_source_label,
        max_files = max_files,
    )

    return plot_rate_curve(curves; reaction_name = reaction_name, title = title)
end

function plot_rate_curve(curves::DataFrame; reaction_name = "reaction", title = nothing)
    if nrow(curves) == 0
        throw(ArgumentError("rate curve data contains no rows"))
    end

    p = plot(
        xlabel = "Temperature T9 (GK)",
        ylabel = "Rate",
        title = title === nothing ? "Rate curve: $reaction_name" : title,
        xscale = :log10,
        yscale = :log10,
        legend = :outertopright,
        xguidefontcolor = :black,
        xtickfontcolor = :black,
        yguidefontcolor = :black,
        ytickfontcolor = :black,
    )

    for group in groupby(curves, [:label, :file, :network_source])
        plot!(
            p,
            group.T9,
            group.rate,
            linewidth = 2,
            label = string(first(group.label), " from ", first(group.file), " (", first(group.network_source), ")"),
        )
    end

    return p
end

function element_z(element_limit)
    if element_limit isa Integer
        return Int(element_limit)
    end

    value = string(element_limit)
    symbol = uppercasefirst(lowercase(value))
    symbol_idx = findfirst(==(symbol), ELEMENT_SYMBOLS)
    symbol_idx !== nothing && return symbol_idx - 1

    name = lowercase(value)
    name_idx = findfirst(==(name), ELEMENT_NAMES)
    name_idx !== nothing && return name_idx - 1

    throw(ArgumentError("unknown element limit: $element_limit"))
end

function element_symbol(Z::Integer)
    if 0 <= Z < length(ELEMENT_SYMBOLS)
        return ELEMENT_SYMBOLS[Z + 1]
    end
    return string("Z=", Z)
end

function parse_network_species(text)
    value = strip(text)

    if value == "PROT"
        return (A = 1, Z = 1, symbol = "H", label = "p")
    elseif value == "NEUT"
        return (A = 1, Z = 0, symbol = "n", label = "n")
    elseif value == "OOOOO"
        return (A = 0, Z = 0, symbol = "", label = "")
    end

    m = match(r"^([A-Z]+)\s*(\d+)$", value)
    m === nothing && return nothing

    symbol = uppercasefirst(lowercase(m.captures[1]))
    zval = element_z(symbol)
    aval = parse(Int, m.captures[2])
    return (A = aval, Z = zval, symbol = symbol, label = string(aval, symbol))
end

function reaction_network_label(reactant, product, rtype)
    if reactant === nothing || product === nothing
        return string(rtype)
    end

    return string(reactant.label, rtype, product.label)
end

function parse_network_float(value)
    text = strip(string(value))
    fixed = occursin(r"[EeDd]", text) ? replace(text, 'D' => 'E', 'd' => 'E') : replace(text, r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$" => s"\1E\2")
    return parse(Float64, fixed)
end

function read_network_reactions(networksetup_path)
    reactions = Dict{Int, NamedTuple}()
    isfile(networksetup_path) || return reactions

    network_re = r"^\s*(\d+)\s+([TF])\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+->\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)"

    for line in eachline(networksetup_path)
        m = match(network_re, line)
        m === nothing && continue

        idx = parse(Int, m.captures[1])
        reactant_count = parse(Int, m.captures[3])
        reactant = parse_network_species(m.captures[4])
        projectile_count = parse(Int, m.captures[5])
        projectile = parse_network_species(m.captures[6])
        product_count = parse(Int, m.captures[7])
        product = parse_network_species(m.captures[8])
        ejectile_count = parse(Int, m.captures[9])
        ejectile = parse_network_species(m.captures[10])
        rate = parse_network_float(m.captures[11])
        source = m.captures[12]
        rtype = m.captures[13]
        chapter = parse(Int, m.captures[14])
        multiplier = parse_network_float(m.captures[15])
        qvalue = parse_network_float(m.captures[16])

        reactions[idx] = (
            index = idx,
            active = m.captures[2] == "T",
            reactant_count = reactant_count,
            reactant = reactant,
            projectile_count = projectile_count,
            projectile = projectile,
            product_count = product_count,
            product = product,
            ejectile_count = ejectile_count,
            ejectile = ejectile,
            rate = rate,
            source = source,
            rtype = rtype,
            chapter = chapter,
            multiplier = multiplier,
            qvalue = qvalue,
            label = reaction_network_label(reactant, product, rtype),
            raw_line = line,
        )
    end

    return reactions
end

function default_networksetup_path(filepath)
    candidate = joinpath(dirname(filepath), "networksetup.txt")
    return isfile(candidate) ? candidate : nothing
end

function value_in_region(value, region)
    region === nothing && return true

    if region isa Tuple && length(region) == 2
        return first(region) <= value <= last(region)
    elseif region isa AbstractRange
        return value in region
    elseif region isa AbstractVector
        return value in region
    elseif region isa Integer
        return value == region
    end

    throw(ArgumentError("region must be nothing, an integer, a range, a vector, or a two-value tuple"))
end

function endpoint_in_region(a_start, z_start, a_end, z_end, a_range, z_range)
    return value_in_region(a_start, a_range) &&
        value_in_region(a_end, a_range) &&
        value_in_region(z_start, z_range) &&
        value_in_region(z_end, z_range)
end

function parse_rate_coefficients(line::AbstractString)
    normalized = replace(String(line), 'D' => 'e', 'd' => 'e')
    return [parse(Float64, m.match) for m in eachmatch(r"[+-]?(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?", normalized)]
end

function normalize_reaclib_species(species)
    value = lowercase(strip(string(species)))
    value in ("prot", "h1") && return "p"
    value in ("neut", "n1") && return "n"
    value in ("ooooo", "g", "") && return ""
    return value
end

function isotope_to_reaclib_species(text)
    value = strip(string(text))
    m = match(r"^(\d+)([A-Za-z]+)([gm]?)$", value)
    m === nothing && throw(ArgumentError("could not parse isotope label: $text"))
    return lowercase(string(m.captures[2], m.captures[1]))
end

function network_species_to_reaclib_species(species)
    species === nothing && return ""
    species.label == "p" && return "p"
    species.label == "n" && return "n"
    species.A == 0 && return ""
    return lowercase(string(species.symbol, species.A))
end

function reaclib_chapter_counts(chapter::Integer, nspecies::Integer)
    if chapter == 1
        return (1, nspecies - 1)
    elseif chapter == 2
        return (1, 2)
    elseif chapter == 3
        return (1, 3)
    elseif chapter == 4
        return (2, 1)
    elseif chapter == 5
        return (2, 2)
    elseif chapter == 6
        return (2, 3)
    elseif chapter == 7
        return (2, 4)
    elseif chapter == 8
        return (3, 1)
    end

    throw(ArgumentError("unsupported REACLIB chapter: $chapter"))
end

function parse_reaclib_species_line(chapter::Integer, line::AbstractString)
    tokens = split(line)
    length(tokens) < 4 && return nothing

    label = String(tokens[end - 1])
    qvalue = try
        parse(Float64, tokens[end])
    catch
        return nothing
    end

    species = normalize_reaclib_species.(tokens[1:end - 2])
    nreactants, nproducts = reaclib_chapter_counts(chapter, length(species))
    nreactants < 0 && return nothing
    length(species) != nreactants + nproducts && return nothing

    return (
        reactants = sort(species[1:nreactants]),
        products = sort(filter(!isempty, species[nreactants + 1:end])),
        label = label,
        qvalue = qvalue,
    )
end

function read_reaclib_rate_blocks(filepath)
    lines = [strip(line) for line in readlines(filepath) if !isempty(strip(line))]
    blocks = NamedTuple[]
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

        parsed = parse_reaclib_species_line(chapter, lines[i])
        coeffs = vcat(parse_rate_coefficients(lines[i + 1]), parse_rate_coefficients(lines[i + 2]))

        if parsed !== nothing && length(coeffs) == 7
            push!(
                blocks,
                (
                    filepath = filepath,
                    chapter = chapter,
                    reactants = parsed.reactants,
                    products = parsed.products,
                    label = parsed.label,
                    qvalue = parsed.qvalue,
                    coefficients = Tuple(coeffs),
                ),
            )
            i += 3
        else
            i += 1
        end
    end

    return blocks
end

function reaclib_rate(coefficients, T9)
    T9 <= 0 && throw(DomainError(T9, "temperature must be positive and in GK"))
    T13 = T9^(1 / 3)
    T53 = T9 * T13 * T13
    a = coefficients
    exponent = a[1] + a[2] / T9 + a[3] / T13 + a[4] * T13 + a[5] * T9 + a[6] * T53 + a[7] * log(T9)

    exponent > log(floatmax(Float64)) && return Inf
    exponent < log(floatmin(Float64)) && return 0.0
    return exp(exponent)
end

function rate_curve_temperature_grid(; Tmin = 0.01, Tmax = 10.0, n = 300)
    Tmin <= 0 && throw(ArgumentError("Tmin must be positive"))
    Tmax <= Tmin && throw(ArgumentError("Tmax must be greater than Tmin"))
    return 10 .^ range(log10(Tmin), log10(Tmax), length = n)
end

function parse_reaction_name_for_rate(reaction_name)
    m = match(r"^(\d+[A-Za-z]+[gm]?)_([a-z]{2})_(\d+[A-Za-z]+[gm]?)$", string(reaction_name))
    m === nothing && throw(ArgumentError("reaction name must look like 18F_pa_15O or 14N_pg_15O"))

    target = isotope_to_reaclib_species(m.captures[1])
    code = lowercase(m.captures[2])
    product = isotope_to_reaclib_species(m.captures[3])

    projectile = if startswith(code, "p")
        "p"
    elseif startswith(code, "a")
        "he4"
    elseif startswith(code, "n")
        "n"
    else
        throw(ArgumentError("unsupported projectile code in reaction name: $reaction_name"))
    end

    products = if endswith(code, "g")
        [product]
    elseif endswith(code, "a")
        sort([product, "he4"])
    elseif endswith(code, "p")
        sort([product, "p"])
    elseif endswith(code, "n")
        sort([product, "n"])
    else
        throw(ArgumentError("unsupported ejectile code in reaction name: $reaction_name"))
    end

    return (reactants = sort([target, projectile]), products = sort(products), code = code)
end

function reaction_matches_network_row(target, row)
    reactants = sort(filter(!isempty, [
        network_species_to_reaclib_species(row.reactant),
        network_species_to_reaclib_species(row.projectile),
    ]))
    products = sort(filter(!isempty, [
        network_species_to_reaclib_species(row.product),
        network_species_to_reaclib_species(row.ejectile),
    ]))

    return reactants == target.reactants && products == target.products
end

function network_species_a(species)
    species === nothing && return missing
    return species.A == 0 ? missing : species.A
end

function network_species_z(species)
    species === nothing && return missing
    return species.A == 0 && species.Z == 0 ? missing : species.Z
end

function network_species_label(species)
    species === nothing && return ""
    return species.label
end

function matching_network_reactions(reaction_name, networksetup_path)
    target = parse_reaction_name_for_rate(reaction_name)
    reactions = read_network_reactions(networksetup_path)
    matches = [(index = idx, row = row) for (idx, row) in reactions if reaction_matches_network_row(target, row)]
    return sort(matches, by = match -> match.index)
end

function network_match_table(reaction_name, networksetup_path)
    matches = matching_network_reactions(reaction_name, networksetup_path)
    n_matches = length(matches)
    active_count = count(match -> match.row.active, matches)

    return DataFrame([
        (
            reaction = string(reaction_name),
            network_index = match.index,
            active = match.row.active,
            source = match.row.source,
            rtype = match.row.rtype,
            chapter = match.row.chapter,
            multiplier = match.row.multiplier,
            rate = match.row.rate,
            qvalue = match.row.qvalue,
            label = match.row.label,
            n_matching_rows = n_matches,
            n_active_matching_rows = active_count,
            reactant = network_species_label(match.row.reactant),
            projectile = network_species_label(match.row.projectile),
            product = network_species_label(match.row.product),
            ejectile = network_species_label(match.row.ejectile),
            a_start = network_species_a(match.row.reactant),
            z_start = network_species_z(match.row.reactant),
            a_end = network_species_a(match.row.product),
            z_end = network_species_z(match.row.product),
            networksetup_path = string(networksetup_path),
            raw_line = match.row.raw_line,
        )
        for match in matches
    ])
end

"""
    reaction_audit(reaction_name; networksetup_path)

Return the `networksetup.txt` rows that match a reaction name such as
`"30P_pg_31S"`. The table includes active status, source, reaction type,
multiplier, Q value, and duplicate-row counts.
"""
function reaction_audit(reaction_name; networksetup_path)
    isfile(networksetup_path) || throw(ArgumentError("could not find networksetup file: $networksetup_path"))
    return network_match_table(reaction_name, networksetup_path)
end

function factor_audit_rows(
    reaction_name;
    reaction_run_path,
    factors,
    networksetup_filename,
    include_baseline,
)
    factor_entries = [(factor = 1.0, label = "baseline", run_path = joinpath(reaction_run_path, "baseline"))]

    if !include_baseline
        factor_entries = NamedTuple[]
    end

    for factor in factors
        factor_label = factor_to_folder(factor)
        push!(
            factor_entries,
            (
                factor = Float64(factor),
                label = factor_label,
                run_path = joinpath(reaction_run_path, string(reaction_name), "fact_$(factor_label)"),
            ),
        )
    end

    rows = NamedTuple[]

    for entry in factor_entries
        networksetup_path = joinpath(entry.run_path, networksetup_filename)

        if !isfile(networksetup_path)
            push!(
                rows,
                (
                    reaction = string(reaction_name),
                    factor = entry.label,
                    expected_multiplier = entry.factor,
                    network_index = missing,
                    active = missing,
                    source = missing,
                    rtype = missing,
                    multiplier = missing,
                    multiplier_over_expected = missing,
                    factor_applied = false,
                    n_matching_rows = 0,
                    n_active_matching_rows = 0,
                    networksetup_path = string(networksetup_path),
                    raw_line = "",
                ),
            )
            continue
        end

        audit = network_match_table(reaction_name, networksetup_path)

        if nrow(audit) == 0
            push!(
                rows,
                (
                    reaction = string(reaction_name),
                    factor = entry.label,
                    expected_multiplier = entry.factor,
                    network_index = missing,
                    active = missing,
                    source = missing,
                    rtype = missing,
                    multiplier = missing,
                    multiplier_over_expected = missing,
                    factor_applied = false,
                    n_matching_rows = 0,
                    n_active_matching_rows = 0,
                    networksetup_path = string(networksetup_path),
                    raw_line = "",
                ),
            )
            continue
        end

        for row in eachrow(audit)
            ratio = row.multiplier / entry.factor
            push!(
                rows,
                (
                    reaction = string(reaction_name),
                    factor = entry.label,
                    expected_multiplier = entry.factor,
                    network_index = row.network_index,
                    active = row.active,
                    source = row.source,
                    rtype = row.rtype,
                    multiplier = row.multiplier,
                    multiplier_over_expected = ratio,
                    factor_applied = row.active && isapprox(row.multiplier, entry.factor; rtol = 1e-8, atol = 0.0),
                    n_matching_rows = row.n_matching_rows,
                    n_active_matching_rows = row.n_active_matching_rows,
                    networksetup_path = row.networksetup_path,
                    raw_line = row.raw_line,
                ),
            )
        end
    end

    return rows
end

"""
    factor_audit(reaction_name; reaction_run_path="../runs", factors=[100,10,2,0.5,0.1,0.01])

Check each factored run for `reaction_name` and report whether the active
matching `networksetup.txt` row has the expected multiplier.
"""
function factor_audit(
    reaction_name;
    reaction_run_path = "../runs",
    factors = [100, 10, 2, 0.5, 0.1, 0.01],
    networksetup_filename = "networksetup.txt",
    include_baseline = true,
)
    return DataFrame(
        factor_audit_rows(
            reaction_name;
            reaction_run_path = reaction_run_path,
            factors = factors,
            networksetup_filename = networksetup_filename,
            include_baseline = include_baseline,
        ),
    )
end

function read_flux_values_by_index(flux_path)
    flux_by_index = Dict{Int, Float64}()
    isfile(flux_path) || return flux_by_index

    for line in eachline(flux_path)
        text = strip(line)
        isempty(text) && continue
        startswith(text, "#") && continue

        parts = split(text)
        length(parts) < 10 && continue

        try
            idx = parse(Int, parts[1])
            flux = parse_network_float(parts[10])
            flux_by_index[idx] = get(flux_by_index, idx, 0.0) + flux
        catch
            continue
        end
    end

    return flux_by_index
end

"""
    flux_audit(reaction_name; flux_path, networksetup_path=default_networksetup_path(flux_path))

Return the flux carried by the `networksetup.txt` rows matching `reaction_name`.
Use `tolerance` to hide smaller fluxes.
"""
function flux_audit(
    reaction_name;
    flux_path,
    networksetup_path = default_networksetup_path(flux_path),
    tolerance = 0.0,
)
    networksetup_path === nothing && throw(ArgumentError("networksetup_path was not provided and could not be inferred from $flux_path"))
    isfile(flux_path) || throw(ArgumentError("could not find flux file: $flux_path"))
    audit = reaction_audit(reaction_name; networksetup_path = networksetup_path)
    flux_by_index = read_flux_values_by_index(flux_path)

    rows = NamedTuple[]
    for row in eachrow(audit)
        flux = get(flux_by_index, row.network_index, 0.0)
        abs(flux) < tolerance && continue

        push!(
            rows,
            (
                reaction = row.reaction,
                network_index = row.network_index,
                active = row.active,
                source = row.source,
                rtype = row.rtype,
                multiplier = row.multiplier,
                flux = flux,
                abs_flux = abs(flux),
                label = row.label,
                flux_path = string(flux_path),
                networksetup_path = row.networksetup_path,
            ),
        )
    end

    return DataFrame(rows)
end

function find_network_rate_row(reaction_name, networksetup_path)
    target = parse_reaction_name_for_rate(reaction_name)
    reactions = read_network_reactions(networksetup_path)
    matches = [(index = idx, row = row) for (idx, row) in reactions if reaction_matches_network_row(target, row)]
    isempty(matches) && throw(ArgumentError("could not find reaction $reaction_name in $networksetup_path"))

    active_matches = filter(match -> match.row.active, matches)
    selected = isempty(active_matches) ? sort(matches, by = match -> match.index)[1] : sort(active_matches, by = match -> match.index)[1]
    return (
        index = selected.index,
        source = selected.row.source,
        rtype = selected.row.rtype,
        label = selected.row.label,
        reactants = target.reactants,
        products = target.products,
    )
end

function candidate_reaclib_files(npdata_path, source)
    reaclib_dir = joinpath(npdata_path, "REACLIB")
    isdir(reaclib_dir) || throw(ArgumentError("could not find REACLIB directory under $npdata_path"))

    files = [joinpath(reaclib_dir, file) for file in readdir(reaclib_dir) if isfile(joinpath(reaclib_dir, file))]
    source_upper = uppercase(string(source))

    priority_names = if source_upper == "JINAC"
        ["reaclib.nosmo", "20120510ReaclibV1.1", "results01111258", "20081109ReaclibV0.5"]
    elseif startswith(source_upper, "NACR")
        ["reaclib.nosmo", "20120510ReaclibV1.1", "20081109ReaclibV0.5"]
    elseif source_upper == "ILI01"
        ["results01111258", "20120510ReaclibV1.1", "reaclib.nosmo"]
    else
        ["reaclib.nosmo", "20120510ReaclibV1.1", "results01111258", "20081109ReaclibV0.5"]
    end

    priority = Dict(name => i for (i, name) in enumerate(priority_names))
    return sort(files, by = file -> get(priority, basename(file), length(priority_names) + 1))
end

function source_matches_reaclib_label(source, label)
    source_upper = uppercase(string(source))
    label_lower = lowercase(string(label))

    if startswith(source_upper, "NACR")
        return startswith(label_lower, lowercase(source_upper))
    elseif source_upper == "ILI01"
        return startswith(label_lower, "il01") || startswith(label_lower, "il10")
    elseif source_upper == "JINAC"
        return true
    end

    return startswith(label_lower, lowercase(source_upper))
end

function read_rate_curve_data(
    reaction_name;
    networksetup_path,
    npdata_path,
    temperatures = rate_curve_temperature_grid(),
    source = nothing,
    match_source_label = true,
    max_files = 4,
)
    network_row = find_network_rate_row(reaction_name, networksetup_path)
    selected_source = source === nothing ? network_row.source : string(source)
    files = candidate_reaclib_files(npdata_path, selected_source)

    matched_blocks = NamedTuple[]
    for file in Iterators.take(files, max_files)
        blocks = read_reaclib_rate_blocks(file)
        append!(
            matched_blocks,
            [
                block for block in blocks
                if block.reactants == network_row.reactants &&
                    block.products == network_row.products &&
                    (!match_source_label || source_matches_reaclib_label(selected_source, block.label))
            ],
        )
        !isempty(matched_blocks) && break
    end

    if isempty(matched_blocks) && match_source_label
        return read_rate_curve_data(
            reaction_name;
            networksetup_path = networksetup_path,
            npdata_path = npdata_path,
            temperatures = temperatures,
            source = selected_source,
            match_source_label = false,
            max_files = max_files,
        )
    end

    isempty(matched_blocks) && throw(ArgumentError("could not find REACLIB rate data for $reaction_name under $npdata_path"))

    grouped = Dict{String, Vector{NamedTuple}}()
    for block in matched_blocks
        push!(get!(grouped, block.label, NamedTuple[]), block)
    end

    T9 = Float64[]
    rate = Float64[]
    label = String[]
    file = String[]
    network_source = String[]
    network_index = Int[]
    qvalue = Float64[]

    for (group_label, blocks) in grouped
        for T in temperatures
            total_rate = sum(reaclib_rate(block.coefficients, T) for block in blocks)
            push!(T9, T)
            push!(rate, total_rate)
            push!(label, group_label)
            push!(file, basename(first(blocks).filepath))
            push!(network_source, network_row.source)
            push!(network_index, network_row.index)
            push!(qvalue, first(blocks).qvalue)
        end
    end

    return DataFrame(
        T9 = T9,
        rate = rate,
        label = label,
        file = file,
        network_source = network_source,
        network_index = network_index,
        qvalue = qvalue,
    )
end

function read_abundance_chart_data(filepath; element_limit = "Ca", tolerance = 1e-10)
    max_z = element_z(element_limit)

    A = Int[]
    N = Int[]
    Z = Int[]
    abundance = Float64[]
    element = String[]
    isotope = String[]

    for line in eachline(filepath)
        line = strip(line)

        isempty(line) && continue
        startswith(line, "#") && continue
        startswith(line, "H NUM") && continue

        parts = split(line)
        length(parts) < 6 && continue

        try
            zval = Int(round(parse(Float64, parts[2])))
            aval = Int(round(parse(Float64, parts[3])))
            xval = parse(Float64, parts[5])

            zval < 1 && continue
            zval > max_z && continue
            xval < tolerance && continue

            symbol = if parts[end] == "PROT"
                "H"
            elseif length(parts) >= 7
                uppercasefirst(lowercase(parts[6]))
            else
                m = match(r"([A-Za-z]+)(\d+)", parts[end])
                m === nothing ? element_symbol(zval) : uppercasefirst(lowercase(m.captures[1]))
            end

            push!(A, aval)
            push!(N, aval - zval)
            push!(Z, zval)
            push!(abundance, xval)
            push!(element, symbol)
            push!(isotope, string(symbol, "-", aval))
        catch
            continue
        end
    end

    return DataFrame(A=A, N=N, Z=Z, abundance=abundance, element=element, isotope=isotope)
end

function read_flow_chart_data(filepath; element_limit = "Ca", tolerance = 1e-10)
    max_z = element_z(element_limit)
    flows = Dict{Tuple{Int, Int, Int, Int, Int, Int, Int, Int}, Float64}()

    for line in eachline(filepath)
        line = strip(line)

        isempty(line) && continue
        startswith(line, "#") && continue

        parts = split(line)
        length(parts) < 10 && continue

        try
            z_start = parse(Int, parts[2])
            a_start = parse(Int, parts[3])
            z_end = parse(Int, parts[8])
            a_end = parse(Int, parts[9])
            flux = parse(Float64, parts[10])

            z_start < 1 && continue
            z_end < 1 && continue
            z_start > max_z && continue
            z_end > max_z && continue
            flux < tolerance && continue

            n_start = a_start - z_start
            n_end = a_end - z_end
            n_start == n_end && z_start == z_end && continue

            key = (n_start, z_start, n_end, z_end, a_start, z_start, a_end, z_end)
            flows[key] = get(flows, key, 0.0) + flux
        catch
            continue
        end
    end

    n_start = Int[]
    z_start = Int[]
    n_end = Int[]
    z_end = Int[]
    a_start = Int[]
    a_end = Int[]
    flux = Float64[]

    for (key, value) in flows
        push!(n_start, key[1])
        push!(z_start, key[2])
        push!(n_end, key[3])
        push!(z_end, key[4])
        push!(a_start, key[5])
        push!(a_end, key[7])
        push!(flux, value)
    end

    return DataFrame(
        n_start = n_start,
        z_start = z_start,
        n_end = n_end,
        z_end = z_end,
        a_start = a_start,
        a_end = a_end,
        flux = flux,
    )
end

function read_flow_region_chart_data(
    filepath;
    a_range = nothing,
    z_range = nothing,
    A_range = nothing,
    Z_range = nothing,
    tolerance = 1e-10,
    networksetup_path = default_networksetup_path(filepath),
)
    selected_a_range = A_range === nothing ? a_range : A_range
    selected_z_range = Z_range === nothing ? z_range : Z_range
    reactions = networksetup_path === nothing ? Dict{Int, NamedTuple}() : read_network_reactions(networksetup_path)

    candidates = NamedTuple[]

    for line in eachline(filepath)
        line = strip(line)

        isempty(line) && continue
        startswith(line, "#") && continue

        parts = split(line)
        length(parts) < 10 && continue

        try
            idx = parse(Int, parts[1])
            z1 = parse(Int, parts[2])
            a1 = parse(Int, parts[3])
            z2 = parse(Int, parts[8])
            a2 = parse(Int, parts[9])
            flow = parse(Float64, parts[10])

            z1 < 1 && continue
            z2 < 1 && continue
            endpoint_in_region(a1, z1, a2, z2, selected_a_range, selected_z_range) || continue

            n1 = a1 - z1
            n2 = a2 - z2
            n1 == n2 && z1 == z2 && continue

            metadata = get(reactions, idx, nothing)
            label = metadata === nothing ? string("reaction ", idx) : metadata.label
            src = metadata === nothing ? "unknown" : metadata.source
            active = metadata === nothing ? true : metadata.active

            push!(
                candidates,
                (
                    reaction_index = idx,
                    n_start = n1,
                    z_start = z1,
                    n_end = n2,
                    z_end = z2,
                    a_start = a1,
                    a_end = a2,
                    flux = flow,
                    reaction_label = label,
                    source = src,
                    active = active,
                ),
            )
        catch
            continue
        end
    end

    duplicate_groups = Dict{Tuple{Int, Int, Int, Int, Int, Int, String}, Vector{NamedTuple}}()

    for candidate in candidates
        key = (
            candidate.n_start,
            candidate.z_start,
            candidate.n_end,
            candidate.z_end,
            candidate.a_start,
            candidate.a_end,
            candidate.reaction_label,
        )
        push!(get!(duplicate_groups, key, NamedTuple[]), candidate)
    end

    endpoint_groups = Dict{Tuple{Int, Int, Int, Int, Int, Int}, Vector{NamedTuple}}()

    for group in values(duplicate_groups)
        active_rows = filter(row -> row.active, group)
        selected_rows = isempty(active_rows) ? group : active_rows
        selected_flux = sum(row.flux for row in selected_rows)
        selected_flux < tolerance && continue

        first_row = first(selected_rows)
        sources = sort(unique([string(row.source) for row in selected_rows]))
        indices = sort(unique([row.reaction_index for row in selected_rows]))
        endpoint_key = (
            first_row.n_start,
            first_row.z_start,
            first_row.n_end,
            first_row.z_end,
            first_row.a_start,
            first_row.a_end,
        )

        push!(
            get!(endpoint_groups, endpoint_key, NamedTuple[]),
            (
                reaction_indices = join(indices, ","),
                flux = selected_flux,
                reaction_label = first_row.reaction_label,
                source = join(sources, "/"),
                annotation_part = string(first_row.reaction_label, "\n", join(sources, "/")),
            ),
        )
    end

    reaction_index = String[]
    n_start = Int[]
    z_start = Int[]
    n_end = Int[]
    z_end = Int[]
    a_start = Int[]
    a_end = Int[]
    flux = Float64[]
    reaction_label = String[]
    source = String[]
    annotation = String[]

    for (endpoint, rows) in endpoint_groups
        sorted_rows = sort(rows, by = row -> row.reaction_label)
        total_flux = sum(row.flux for row in sorted_rows)
        parts = [row.annotation_part for row in sorted_rows]

        push!(reaction_index, join([row.reaction_indices for row in sorted_rows], ","))
        push!(n_start, endpoint[1])
        push!(z_start, endpoint[2])
        push!(n_end, endpoint[3])
        push!(z_end, endpoint[4])
        push!(a_start, endpoint[5])
        push!(a_end, endpoint[6])
        push!(flux, total_flux)
        push!(reaction_label, join([row.reaction_label for row in sorted_rows], ", "))
        push!(source, join([row.source for row in sorted_rows], ", "))
        push!(annotation, join(parts, ",\n"))
    end

    return DataFrame(
        reaction_index = reaction_index,
        n_start = n_start,
        z_start = z_start,
        n_end = n_end,
        z_end = z_end,
        a_start = a_start,
        a_end = a_end,
        flux = flux,
        reaction_label = reaction_label,
        source = source,
        annotation = annotation,
    )
end

function chart_tile_color(value, color_limits; target_color = CM.RGBAf(1.0, 0.0, 0.0, 1.0))
    color_fraction = clamp((value - color_limits[1]) / (color_limits[2] - color_limits[1]), 0.0, 1.0)
    return CM.RGBAf(
        1.0 + color_fraction * (target_color.r - 1.0),
        1.0 + color_fraction * (target_color.g - 1.0),
        1.0 + color_fraction * (target_color.b - 1.0),
        1.0,
    )
end

function chart_colormap_color(value, color_limits, colormap)
    color_fraction = clamp((value - color_limits[1]) / (color_limits[2] - color_limits[1]), 0.0, 1.0)
    return CM.cgrad(colormap)[color_fraction]
end

function draw_isotope_tiles!(
    ax,
    isotopes::DataFrame;
    value_col,
    color_limits,
    mass_label_size = 8,
)
    for row in eachrow(isotopes)
        value = log10(max(row[value_col], 10.0^color_limits[1]))
        corners = CM.Point2f[
            (row.N - 0.5, row.Z - 0.5),
            (row.N + 0.5, row.Z - 0.5),
            (row.N + 0.5, row.Z + 0.5),
            (row.N - 0.5, row.Z + 0.5),
        ]

        CM.poly!(
            ax,
            corners,
            color = chart_tile_color(value, color_limits),
            strokecolor = :black,
            strokewidth = 1,
        )

        CM.text!(
            ax,
            string(row.A),
            position = (row.N, row.Z),
            align = (:center, :center),
            fontsize = mass_label_size,
            color = :black,
        )
    end
end

function add_element_labels!(
    ax,
    row_data::DataFrame,
    min_n,
    max_z;
    min_z = 1,
    element_label_size = 16,
)
    for z in min_z:max_z
        rows = filter(row -> row.Z == z, row_data)
        n_label = nrow(rows) == 0 ? min_n - 1.5 : minimum(rows.N) - 0.75

        CM.text!(
            ax,
            element_symbol(z),
            position = (n_label, z),
            align = (:right, :center),
            fontsize = element_label_size,
            font = :bold,
            color = :black,
        )
    end
end

#=
function abundance_chart(filepath; element_limit = "Ca", tolerance = 1e-10, title = "Abundance Chart")
    abundances = read_abundance_chart_data(
        filepath;
        element_limit = element_limit,
        tolerance = tolerance,
    )
    return abundance_chart(abundances; element_limit = element_limit, tolerance = tolerance, title = title)
end

function abundance_chart(abundances::DataFrame; element_limit = "Ca", tolerance = 1e-10, title = "Abundance Chart")
    if nrow(abundances) == 0
        throw(ArgumentError("no isotopes found at or above tolerance=$tolerance up to element_limit=$element_limit"))
    end

    max_z = element_z(element_limit)
    min_n = minimum(abundances.N)
    max_n = maximum(abundances.N)
    n_values = collect(min_n:max_n)
    z_values = collect(0:max_z)
    log_abundance_grid = fill(NaN, length(z_values), length(n_values))

    n_index = Dict(n => i for (i, n) in enumerate(n_values))
    z_index = Dict(z => i for (i, z) in enumerate(z_values))

    for row in eachrow(abundances)
        log_abundance_grid[z_index[row.Z], n_index[row.N]] = log10(max(row.abundance, tolerance))
    end

    color_limits = (log10(tolerance), 0.0)

    p = heatmap(
        n_values,
        z_values,
        log_abundance_grid,
        color = cgrad([:white, :red]),
        clims = color_limits,
        xlabel = "neutron number (A-Z)",
        ylabel = "proton number (Z)",
        title = title,
        colorbar = true,
        colorbar_title = "log10(X)",
        xlims = (min_n - 3, max_n + 1),
        ylims = (-0.5, max_z + 1.0),
        xticks = min_n:max_n,
        yticks = 0:2:max_z,
        aspect_ratio = :equal,
        grid = false,
        legend = false,
    )

    for row in eachrow(abundances)
        plot!(
            p,
            [row.N - 0.5, row.N + 0.5, row.N + 0.5, row.N - 0.5, row.N - 0.5],
            [row.Z - 0.5, row.Z - 0.5, row.Z + 0.5, row.Z + 0.5, row.Z - 0.5],
            color = :black,
            linewidth = 1,
            label = "",
        )
    end

    for row in eachrow(abundances)
        annotate!(p, row.N, row.Z, text(string(row.A), 6, :black, :center))
    end

    for z in 1:max_z
        row_abundances = filter(row -> row.Z == z, abundances)
        n_label = nrow(row_abundances) == 0 ? min_n - 1.5 : minimum(row_abundances.N) - 1.25
        annotate!(p, n_label + 0.6, z, text(element_symbol(z), 7, :black, :right, "bold"))
        #annotate!(p, n_label + 0.02, z, text(element_symbol(z), 10, :black, :right))
    end

    return p
end
=#

function abundance_chart(
    filepath;
    element_limit = "Ca",
    tolerance = 1e-10,
    title = "Abundance Chart",
    figure_size = (900, 650),
    element_label_size = 16,
    mass_label_size = 8,
)
    abundances = read_abundance_chart_data(
        filepath;
        element_limit = element_limit,
        tolerance = tolerance,
    )
    return abundance_chart(
        abundances;
        element_limit = element_limit,
        tolerance = tolerance,
        title = title,
        figure_size = figure_size,
        element_label_size = element_label_size,
        mass_label_size = mass_label_size,
    )
end

function abundance_chart(
    abundances::DataFrame;
    element_limit = "Ca",
    tolerance = 1e-10,
    title = "Abundance Chart",
    figure_size = (900, 650),
    element_label_size = 16,
    mass_label_size = 12,
)
    if nrow(abundances) == 0
        throw(ArgumentError("no isotopes found at or above tolerance=$tolerance up to element_limit=$element_limit"))
    end

    max_z = element_z(element_limit)
    min_n = minimum(abundances.N)
    max_n = maximum(abundances.N)
    min_x = min_n - 3.0
    max_x = max_n + 1.0
    min_y = -0.5
    max_y = max_z + 1.0
    color_limits = (log10(tolerance), 0.0)

    fig = CM.Figure(size = figure_size)
    ax = CM.Axis(
        fig[1, 1],
        xlabel = "neutron number (A-Z)",
        ylabel = "proton number (Z)",
        title = title,
        aspect = CM.DataAspect(),
        xgridvisible = false,
        ygridvisible = false,
        xticks = min_n:max_n,
        yticks = 0:2:max_z,
        limits = (min_x, max_x, min_y, max_y),
    )

    colormap = [:white, :red]

    for row in eachrow(abundances)
        value = log10(max(row.abundance, tolerance))
        color_fraction = clamp((value - color_limits[1]) / (color_limits[2] - color_limits[1]), 0.0, 1.0)
        tile_color = CM.RGBAf(1.0, 1.0 - color_fraction, 1.0 - color_fraction, 1.0)
        corners = CM.Point2f[
            (row.N - 0.5, row.Z - 0.5),
            (row.N + 0.5, row.Z - 0.5),
            (row.N + 0.5, row.Z + 0.5),
            (row.N - 0.5, row.Z + 0.5),
        ]

        CM.poly!(
            ax,
            corners,
            color = tile_color,
            strokecolor = :black,
            strokewidth = 1,
        )

        CM.text!(
            ax,
            string(row.A),
            position = (row.N, row.Z),
            align = (:center, :center),
            fontsize = mass_label_size,
            color = :black,
        )
    end

    for z in 1:max_z
        row_abundances = filter(row -> row.Z == z, abundances)
        n_label = nrow(row_abundances) == 0 ? min_n - 1.5 : minimum(row_abundances.N) - 0.75

        CM.text!(
            ax,
            element_symbol(z),
            position = (n_label, z),
            align = (:right, :center),
            fontsize = element_label_size,
            font = :bold,
            color = :black,
        )
    end

    CM.Colorbar(
        fig[1, 2],
        colormap = colormap,
        colorrange = color_limits,
        label = "log10(X)",
    )

    return fig
end

function flow_tile_data(flows::DataFrame)
    tiles = Dict{Tuple{Int, Int}, Int}()

    for row in eachrow(flows)
        tiles[(row.n_start, row.z_start)] = row.a_start
        tiles[(row.n_end, row.z_end)] = row.a_end
    end

    N = Int[]
    Z = Int[]
    A = Int[]

    for ((n, z), a) in tiles
        push!(N, n)
        push!(Z, z)
        push!(A, a)
    end

    return DataFrame(N=N, Z=Z, A=A)
end

function draw_empty_isotope_tiles!(ax, tiles::DataFrame; mass_label_size = 8)
    for row in eachrow(tiles)
        corners = CM.Point2f[
            (row.N - 0.5, row.Z - 0.5),
            (row.N + 0.5, row.Z - 0.5),
            (row.N + 0.5, row.Z + 0.5),
            (row.N - 0.5, row.Z + 0.5),
        ]

        CM.poly!(
            ax,
            corners,
            color = :white,
            strokecolor = :black,
            strokewidth = 1,
        )

        CM.text!(
            ax,
            string(row.A),
            position = (row.N, row.Z),
            align = (:center, :center),
            fontsize = mass_label_size,
            color = :black,
        )
    end
end

function draw_arrow!(
    ax,
    x1,
    y1,
    x2,
    y2;
    color,
    linewidth = 2,
    head_length = 0.22,
    head_width = 0.16,
)
    dx = x2 - x1
    dy = y2 - y1
    distance = hypot(dx, dy)
    distance == 0 && return

    ux = dx / distance
    uy = dy / distance
    px = -uy
    py = ux

    start_margin = min(0.35, 0.25 * distance)
    end_margin = min(0.45, 0.35 * distance)
    sx = x1 + start_margin * ux
    sy = y1 + start_margin * uy
    tip_x = x2 - end_margin * ux
    tip_y = y2 - end_margin * uy
    base_x = tip_x - head_length * ux
    base_y = tip_y - head_length * uy

    CM.lines!(
        ax,
        [sx, base_x],
        [sy, base_y],
        color = color,
        linewidth = linewidth,
    )

    CM.poly!(
        ax,
        CM.Point2f[
            (tip_x, tip_y),
            (base_x + 0.5 * head_width * px, base_y + 0.5 * head_width * py),
            (base_x - 0.5 * head_width * px, base_y - 0.5 * head_width * py),
        ],
        color = color,
        strokecolor = color,
    )
end

function draw_arrow_label!(
    ax,
    x1,
    y1,
    x2,
    y2,
    label;
    fontsize = 8,
    offset = 0.18,
)
    dx = x2 - x1
    dy = y2 - y1
    distance = hypot(dx, dy)
    distance == 0 && return

    px = -dy / distance
    py = dx / distance
    angle = atan(dy, dx)

    CM.text!(
        ax,
        label,
        position = ((x1 + x2) / 2 + offset * px, (y1 + y2) / 2 + offset * py),
        align = (:center, :center),
        rotation = angle,
        fontsize = fontsize,
        color = :black,
    )
end

function flow_chart(
    filepath;
    element_limit = "Ca",
    tolerance = 1e-10,
    title = "Flow Chart",
    figure_size = (900, 650),
    element_label_size = 16,
    mass_label_size = 8,
    arrow_linewidth = 2,
)
    flows = read_flow_chart_data(
        filepath;
        element_limit = element_limit,
        tolerance = tolerance,
    )
    return flow_chart(
        flows;
        element_limit = element_limit,
        tolerance = tolerance,
        title = title,
        figure_size = figure_size,
        element_label_size = element_label_size,
        mass_label_size = mass_label_size,
        arrow_linewidth = arrow_linewidth,
    )
end

function flow_chart(
    flows::DataFrame;
    element_limit = "Ca",
    tolerance = 1e-10,
    title = "Flow Chart",
    figure_size = (900, 650),
    element_label_size = 16,
    mass_label_size = 12,
    arrow_linewidth = 4,
)
    if nrow(flows) == 0
        throw(ArgumentError("no fluxes found at or above tolerance=$tolerance up to element_limit=$element_limit"))
    end

    max_z = element_z(element_limit)
    tiles = flow_tile_data(flows)
    min_n = minimum(tiles.N)
    max_n = maximum(tiles.N)
    min_x = min_n - 3.0
    max_x = max_n + 1.0
    min_y = -0.5
    max_y = max_z + 1.0
    log_flux = log10.(max.(flows.flux, tolerance))
    color_limits = (log10(tolerance), maximum(log_flux))
    color_limits[1] == color_limits[2] && (color_limits = (color_limits[1], color_limits[1] + 1.0))
    colormap = :heat

    fig = CM.Figure(size = figure_size)
    ax = CM.Axis(
        fig[1, 1],
        xlabel = "neutron number (A-Z)",
        ylabel = "proton number (Z)",
        title = title,
        aspect = CM.DataAspect(),
        xgridvisible = false,
        ygridvisible = false,
        xticks = min_n:max_n,
        yticks = 0:2:max_z,
        limits = (min_x, max_x, min_y, max_y),
    )

    draw_empty_isotope_tiles!(ax, tiles; mass_label_size = mass_label_size)

    for (row, value) in zip(eachrow(flows), log_flux)
        arrow_color = chart_colormap_color(value, color_limits, colormap)
        draw_arrow!(
            ax,
            row.n_start,
            row.z_start,
            row.n_end,
            row.z_end;
            color = arrow_color,
            linewidth = arrow_linewidth,
        )
    end

    for row in eachrow(tiles)
        CM.text!(
            ax,
            string(row.A),
            position = (row.N, row.Z),
            align = (:center, :center),
            fontsize = mass_label_size,
            color = :black,
        )
    end

    add_element_labels!(
        ax,
        tiles,
        min_n,
        max_z;
        element_label_size = element_label_size,
    )

    CM.Colorbar(
        fig[1, 2],
        colormap = colormap,
        colorrange = color_limits,
        label = "log10(flux)",
    )

    return fig
end

function flow_region_chart(
    filepath;
    a_range = nothing,
    z_range = nothing,
    A_range = nothing,
    Z_range = nothing,
    tolerance = 1e-10,
    networksetup_path = default_networksetup_path(filepath),
    title = "Regional Flow Chart",
    figure_size = (900, 650),
    element_label_size = 16,
    mass_label_size = 12,
    arrow_linewidth = 4,
    reaction_label_size = 8,
    show_reaction_labels = true,
    show_reaction_info = show_reaction_labels,
)
    selected_a_range = A_range === nothing ? a_range : A_range
    selected_z_range = Z_range === nothing ? z_range : Z_range

    flows = read_flow_region_chart_data(
        filepath;
        a_range = selected_a_range,
        z_range = selected_z_range,
        tolerance = tolerance,
        networksetup_path = networksetup_path,
    )

    return flow_region_chart(
        flows;
        a_range = selected_a_range,
        z_range = selected_z_range,
        tolerance = tolerance,
        title = title,
        figure_size = figure_size,
        element_label_size = element_label_size,
        mass_label_size = mass_label_size,
        arrow_linewidth = arrow_linewidth,
        reaction_label_size = reaction_label_size,
        show_reaction_labels = show_reaction_info,
    )
end

function flow_region_chart(
    flows::DataFrame;
    a_range = nothing,
    z_range = nothing,
    A_range = nothing,
    Z_range = nothing,
    tolerance = 1e-10,
    title = "Regional Flow Chart",
    figure_size = (900, 650),
    element_label_size = 16,
    mass_label_size = 12,
    arrow_linewidth = 4,
    reaction_label_size = 8,
    show_reaction_labels = true,
    show_reaction_info = show_reaction_labels,
)
    if nrow(flows) == 0
        throw(ArgumentError("no fluxes found at or above tolerance=$tolerance in the requested A/Z region"))
    end

    tiles = flow_tile_data(flows)
    min_n = minimum(tiles.N)
    max_n = maximum(tiles.N)
    min_z = minimum(tiles.Z)
    max_z = maximum(tiles.Z)
    min_x = min_n - 1.0
    max_x = max_n + 1.0
    min_y = min_z - 0.75
    max_y = max_z + 0.75
    log_flux = log10.(max.(flows.flux, tolerance))
    color_limits = (log10(tolerance), maximum(log_flux))
    color_limits[1] == color_limits[2] && (color_limits = (color_limits[1], color_limits[1] + 1.0))
    colormap = :heat

    fig = CM.Figure(size = figure_size)
    ax = CM.Axis(
        fig[1, 1],
        xlabel = "neutron number (A-Z)",
        ylabel = "proton number (Z)",
        title = title,
        aspect = CM.DataAspect(),
        xgridvisible = false,
        ygridvisible = false,
        xticks = min_n:max_n,
        yticks = min_z:max_z,
        limits = (min_x, max_x, min_y, max_y),
    )

    draw_empty_isotope_tiles!(ax, tiles; mass_label_size = mass_label_size)

    sorted_flows = sort(flows, :flux)
    sorted_log_flux = log10.(max.(sorted_flows.flux, tolerance))

    for (row, value) in zip(eachrow(sorted_flows), sorted_log_flux)
        arrow_color = chart_colormap_color(value, color_limits, colormap)
        draw_arrow!(
            ax,
            row.n_start,
            row.z_start,
            row.n_end,
            row.z_end;
            color = arrow_color,
            linewidth = arrow_linewidth,
        )
    end

    if show_reaction_info
        for row in eachrow(sorted_flows)
            draw_arrow_label!(
                ax,
                row.n_start,
                row.z_start,
                row.n_end,
                row.z_end,
                row.annotation;
                fontsize = reaction_label_size,
            )
        end
    end

    for row in eachrow(tiles)
        CM.text!(
            ax,
            string(row.A),
            position = (row.N, row.Z),
            align = (:center, :center),
            fontsize = mass_label_size,
            color = :black,
        )
    end

    add_element_labels!(
        ax,
        tiles,
        min_n,
        max_z;
        min_z = min_z,
        element_label_size = element_label_size,
    )

    CM.Colorbar(
        fig[1, 2],
        colormap = colormap,
        colorrange = color_limits,
        label = "log10(flux)",
    )

    return fig
end

function analyze_factor(df_compare, f)

    # -------------------------
    # COLUMN NAMES
    # -------------------------
    col_ppn = Symbol(f * "_ppn")
    col_iliadis = Symbol(f * "_iliadis")
    col_ratio = Symbol(f * "_ratio")

    # -------------------------
    # FILTER VALID ROWS
    # -------------------------
    df_valid = filter(row ->
        !ismissing(row[col_ppn]) &&
        !ismissing(row[col_iliadis]),
        df_compare
    )

    # -------------------------
    # COMPUTE METRICS
    # -------------------------
    df_valid.ratio = df_valid[!, col_ppn] ./ df_valid[!, col_iliadis]

    df_valid.logratio = log10.(df_valid.ratio)
    df_valid.dev = abs.(df_valid.logratio)
    styles = reaction_styles(df_valid.reaction)

    # -------------------------
    # SCATTER: PPN vs Iliadis
    # -------------------------
    p1 = plot(
        xlabel = "Iliadis",
        ylabel = "PPN",
        title = "PPN vs Iliadis (factor = $f)",
        xscale = :log10,
        yscale = :log10,
        legend = false,
    )

    scatter_reactions!(
        p1,
        df_valid,
        df_valid[!, col_iliadis],
        df_valid[!, col_ppn];
        styles = styles,
        markerstrokecolor = :auto,
    )
    plot!(p1, [1e-3, 10], [1e-3, 10], linestyle=:dash, label="")

    # -------------------------
    # GROUP BY ISOTOPE
    # -------------------------
    isotopes = unique(df_valid.isotope)
    iso_index = Dict(iso => i for (i, iso) in enumerate(isotopes))
    x = [iso_index[iso] for iso in df_valid.isotope]

    # -------------------------
    # LOG RATIO PLOT
    # -------------------------
    p2 = plot(
        xticks = (1:length(isotopes), isotopes),
        xlabel = "Isotope",
        ylabel = "log10(PPN / Iliadis)",
        title = "Log Ratio (factor = $f)",
        xrotation = 60,
        size = (800,600),
        ylims = (-1, 3),
        legend = false,
    )

    scatter_reactions!(
        p2,
        df_valid,
        x,
        df_valid.logratio;
        styles = styles,
        markerstrokecolor = :auto,
    )
    hline!(p2, [0.0], linestyle=:dash, label = "")

    # -------------------------
    # DEVIATION PLOT
    # -------------------------
    p3 = plot(
        xticks = (1:length(isotopes), isotopes),
        xlabel = "Isotope",
        ylabel = "|log10(PPN / Iliadis)|",
        title = "Deviation (factor = $f)",
        xrotation = 60,
        size = (800, 600),
        legend = false,
    )

    scatter_reactions!(
        p3,
        df_valid,
        x,
        df_valid.dev;
        styles = styles,
        markerstrokecolor = :auto,
    )

    return p1, p2, p3, df_valid
end
