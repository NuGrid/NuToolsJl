# flux_chart.jl — the N/Z chart of reaction fluxes between isotopes.

"""
    flux_chart(flux_df::DataFrame, ab::Abundances; element_limit = "Ca", tolerance = 1e-10,
               show_abundance = true, n_range = nothing, z_range = nothing,
               title = "Flux Chart", figure_size = (900, 650), element_label_size = 16,
               mass_label_size = 8, arrow_linewidth = 2) -> CM.Figure

N/Z chart of `flux_df` (as returned by [`read_fluxes`](@ref)/[`fluxes`](@ref)):
one tile per isotope touched by a reaction, connected by arrows colored
`log10(flux)` on [`FLUX_COLORMAP`](@ref). Tiles are colored by `ab`'s mass
fraction when `show_abundance = true` (so flux and abundance read together),
or left white otherwise. `z_range`/`n_range` (each a `(lo, hi)` tuple) zoom to
a region instead of the whole chart up to `element_limit`.
"""
function flux_chart(flux_df::DataFrame, ab::Abundances; element_limit = "Ca", tolerance = 1e-10,
                     show_abundance = true, n_range = nothing, z_range = nothing,
                     title = "Flux Chart", figure_size = (900, 650), element_label_size = 16,
                     mass_label_size = 8, arrow_linewidth = 2)
    max_z = proton_number(element_limit)
    max_z === nothing && throw(ArgumentError("unknown element_limit \"$element_limit\""))
    z_lo, z_hi = z_range === nothing ? (1, max_z) : z_range

    df = filter(row -> row.flux >= tolerance && z_lo <= row.z_start <= z_hi && z_lo <= row.z_end <= z_hi, flux_df)
    if n_range !== nothing
        n_lo, n_hi = n_range
        df = filter(row -> n_lo <= row.n_start <= n_hi && n_lo <= row.n_end <= n_hi, df)
    end
    nrow(df) == 0 && throw(ArgumentError(
        "no fluxes at or above tolerance=$tolerance in the requested region"))

    tiles = flow_tile_data(df)
    min_n, max_n = extrema(tiles.N)
    log_flux = log10.(max.(df.flux, tolerance))
    flux_limits = (log10(tolerance), maximum(log_flux))
    flux_limits[1] == flux_limits[2] && (flux_limits = (flux_limits[1], flux_limits[1] + 1.0))

    zlo_axis = z_range === nothing ? 0 : z_lo

    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1];
            xlabel = "neutron number (A-Z)", ylabel = "proton number (Z)", title,
            aspect = CM.DataAspect(), xgridvisible = false, ygridvisible = false,
            xticks = min_n:max_n, yticks = zlo_axis:2:z_hi,
            limits = (min_n - 3.0, max_n + 1.0, zlo_axis - 0.5, z_hi + 1.0))

        if show_abundance
            abund_lookup = Dict((row.N, row.Z) => row.X for row in eachrow(DataFrame(ab)) if row.isomer == 0)
            tiles.X = [get(abund_lookup, (row.N, row.Z), tolerance) for row in eachrow(tiles)]
            abund_limits = (log10(tolerance), 0.0)
            draw_isotope_tiles!(ax, tiles; value_col = :X, color_limits = abund_limits,
                                 colormap = ABUNDANCE_COLORMAP, mass_label_size)
            CM.Colorbar(fig[1, 3]; colormap = ABUNDANCE_COLORMAP, colorrange = abund_limits, label = "log10(X)")
        else
            draw_empty_isotope_tiles!(ax, tiles; mass_label_size)
        end

        for (row, value) in zip(eachrow(df), log_flux)
            draw_arrow!(ax, row.n_start, row.z_start, row.n_end, row.z_end;
                        color = chart_colormap_color(value, flux_limits, FLUX_COLORMAP),
                        linewidth = arrow_linewidth)
        end
        add_element_labels!(ax, tiles, min_n, z_hi; min_z = max(1, zlo_axis), element_label_size)
        CM.Colorbar(fig[1, 2]; colormap = FLUX_COLORMAP, colorrange = flux_limits, label = "log10(flux)")
        fig
    end
end

"""
    flux_chart(run::PPNRun, cycle = :final; kwargs...) -> CM.Figure

Flux chart for `run` at `cycle` (a cycle number, or `:initial`/`:final`),
overlaid on that cycle's abundances.
"""
flux_chart(run::PPNRun, cycle = :final; kwargs...) = flux_chart(fluxes(run, cycle), abundances(run, cycle); kwargs...)
