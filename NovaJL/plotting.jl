using Plots
import CairoMakie as CM

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

function scatter_reactions!(p, df, x, y; styles, ylabel_prefix = "", markerstrokecolor = :auto)
    for df_reaction in groupby(df, :reaction)
        reaction = first(df_reaction.reaction)
        style = reaction_style(styles, reaction)
        rows = parentindices(df_reaction)[1]

        scatter!(
            p,
            x[rows],
            y[rows],
            label = string(ylabel_prefix, reaction),
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
        yguidefontcolor = :blue,
        ytickfontcolor = :blue,
        legend = :topright,
    )
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
    mass_label_size = 8,
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
        legend = :outerright,
    )

    scatter_reactions!(
        p1,
        df_valid,
        df_valid[!, col_iliadis],
        df_valid[!, col_ppn];
        styles = styles,
        markerstrokecolor = :auto,
    )
    plot!(p1, [1e-3, 10], [1e-3, 10], linestyle=:dash, label="y = x")

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
        size = (900, 400),
        ylims = (-1, 3),
        legend = :outerright,
    )

    scatter_reactions!(
        p2,
        df_valid,
        x,
        df_valid.logratio;
        styles = styles,
        markerstrokecolor = :auto,
    )
    hline!(p2, [0.0], linestyle=:dash, label = "0")

    # -------------------------
    # DEVIATION PLOT
    # -------------------------
    p3 = plot(
        xticks = (1:length(isotopes), isotopes),
        xlabel = "Isotope",
        ylabel = "|log10(PPN / Iliadis)|",
        title = "Deviation (factor = $f)",
        xrotation = 60,
        size = (900, 400),
        legend = :outerright,
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
