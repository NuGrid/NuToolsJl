using Plots

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
