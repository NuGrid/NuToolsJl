# flux_chart.jl — the N/Z chart of reaction fluxes between isotopes.

"""
    flux_chart(flux_df::DataFrame, ab::Abundances; element_limit = "Ca", tolerance = 1e-10,
               show_abundance = true, n_range = nothing, z_range = nothing,
               title = "Flux Chart", figure_size = (900, 650), element_label_size = 16,
               mass_label_size = 8, arrow_linewidth = 2,
               show_labels = false, net = nothing, label_fontsize = 8) -> CM.Figure

N/Z chart of `flux_df` (as returned by [`read_fluxes`](@ref)/[`fluxes`](@ref)):
one tile per isotope `ab` tracks within `element_limit`/`z_range`/`n_range`
— every tracked isotope gets a tile and a mass-number label regardless of
whether it carries any flux, so a blank patch always means "not tracked,"
never "tracked but quiet" — connected by arrows wherever flux `>= tolerance`,
colored `log10(flux)` on [`FLUX_COLORMAP`](@ref). Tiles are colored by `ab`'s
mass fraction when `show_abundance = true` (so flux and abundance read
together), or left white otherwise. `z_range`/`n_range` (each a `(lo, hi)`
tuple) zoom to a region instead of the whole chart up to `element_limit`.

With `show_labels = true`, each arrow is annotated with its reaction label
and source (e.g. `"12C(p,g)13N  JINAC"`, via [`label`](@ref) and
[`Reaction.source`](@ref Reaction)), which needs `net` (a [`Network`](@ref),
e.g. `network(run)`) to look reactions up by index — only sensible on a
tightly zoomed chart (e.g. `z_range = (5, 10)`), since every arrow on a
whole-chart flux plot would be unreadable clutter.
"""
function flux_chart(flux_df::DataFrame, ab::Abundances; element_limit = "Ca", tolerance = 1e-10,
                     show_abundance = true, n_range = nothing, z_range = nothing,
                     title = "Flux Chart", figure_size = (900, 650), element_label_size = 16,
                     mass_label_size = 8, arrow_linewidth = 2,
                     show_labels = false, net = nothing, label_fontsize = 8)
    show_labels && net === nothing && throw(ArgumentError("show_labels=true needs `net` (a Network) to look reactions up by index"))
    max_z = proton_number(element_limit)
    max_z === nothing && throw(ArgumentError("unknown element_limit \"$element_limit\""))
    z_lo, z_hi = z_range === nothing ? (1, max_z) : z_range

    tiles = filter(:Z => z -> z_lo <= z <= z_hi, DataFrame(ab))
    filter!(:isomer => ==(0), tiles)
    if n_range !== nothing
        n_lo, n_hi = n_range
        filter!(:N => n -> n_lo <= n <= n_hi, tiles)
    end
    nrow(tiles) == 0 && throw(ArgumentError("no isotopes tracked in the requested region"))
    min_n, max_n = extrema(tiles.N)

    arrows = filter(row -> row.flux >= tolerance && z_lo <= row.z_start <= z_hi && z_lo <= row.z_end <= z_hi, flux_df)
    if n_range !== nothing
        n_lo, n_hi = n_range
        filter!(row -> n_lo <= row.n_start <= n_hi && n_lo <= row.n_end <= n_hi, arrows)
    end
    log_flux = isempty(arrows) ? Float64[] : log10.(max.(arrows.flux, tolerance))
    flux_limits = isempty(log_flux) ? (log10(tolerance), log10(tolerance) + 1.0) : (log10(tolerance), maximum(log_flux))
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
            abund_limits = (log10(tolerance), 0.0)
            draw_isotope_tiles!(ax, tiles; value_col = :X, color_limits = abund_limits,
                                 colormap = ABUNDANCE_COLORMAP, mass_label_size)
            CM.Colorbar(fig[1, 3]; colormap = ABUNDANCE_COLORMAP, colorrange = abund_limits, label = "log10(X)")
        else
            draw_empty_isotope_tiles!(ax, tiles; mass_label_size)
        end

        reactions_by_index = show_labels ? Dict(r.index => r for r in net.reactions) : nothing

        for (row, value) in zip(eachrow(arrows), log_flux)
            draw_arrow!(ax, row.n_start, row.z_start, row.n_end, row.z_end;
                        color = chart_colormap_color(value, flux_limits, FLUX_COLORMAP),
                        linewidth = arrow_linewidth)
            if show_labels
                r = get(reactions_by_index, row.index, nothing)
                r === nothing && continue
                draw_arrow_label!(ax, row.n_start, row.z_start, row.n_end, row.z_end,
                                   "$(label(r))  $(r.source)"; fontsize = label_fontsize)
            end
        end
        add_element_labels!(ax, tiles, min_n, z_hi; min_z = max(1, zlo_axis), element_label_size)
        CM.Colorbar(fig[1, 2]; colormap = FLUX_COLORMAP, colorrange = flux_limits, label = "log10(flux)")
        fig
    end
end

"""
    flux_chart(run::PPNRun, cycle = :final; kwargs...) -> CM.Figure

Flux chart for `run` at `cycle` (a cycle number, or `:initial`/`:final`),
overlaid on that cycle's abundances. `net` (needed for `show_labels = true`)
is supplied automatically from `run` — no need to pass it yourself.
"""
flux_chart(run::PPNRun, cycle = :final; kwargs...) =
    flux_chart(fluxes(run, cycle), abundances(run, cycle); net = network(run), kwargs...)
