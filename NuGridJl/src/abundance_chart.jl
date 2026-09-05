# abundance_chart.jl — the N/Z nuclear-chart heatmap of mass fractions.

"""
    abundance_chart(ab::Abundances; element_limit = "Ca", tolerance = 1e-10,
                     hide_below_tolerance = false,
                     title = "Abundance Chart", figure_size = (900, 650),
                     element_label_size = 16, mass_label_size = 8) -> CM.Figure

N/Z tile chart of `ab`'s mass fractions, colored `log10(X)` on
[`ABUNDANCE_COLORMAP`](@ref), mass number in each tile, element symbols along
the left edge. `element_limit` (an element symbol) is a hard cutoff — heavier
isotopes aren't drawn at all.

By default (`hide_below_tolerance = false`) `tolerance` is only a *color*
floor, not a visibility cutoff: every isotope `ab` tracks within
`element_limit` gets a tile and a mass-number label, even if its abundance is
below `tolerance` (it just colors at the floor, indistinguishable from zero)
— so a blank patch on the chart always means "not tracked," never "tracked
but negligible." Pass `hide_below_tolerance = true` to instead drop those
tiles entirely (no tile, no label) — useful for decluttering a chart down to
only the isotopes that actually matter at a given threshold, at the cost of
that "blank = not tracked" guarantee (a blank patch can now also mean
"tracked, but below tolerance").
"""
function abundance_chart(ab::Abundances; element_limit = "Ca", tolerance = 1e-10,
                          hide_below_tolerance::Bool = false,
                          title = "Abundance Chart", figure_size = (900, 650),
                          element_label_size = 16, mass_label_size = 8)
    max_z = proton_number(element_limit)
    max_z === nothing && throw(ArgumentError("unknown element_limit \"$element_limit\""))

    df = filter(:Z => z -> 1 <= z <= max_z, DataFrame(ab))
    hide_below_tolerance && filter!(:X => x -> x >= tolerance, df)
    nrow(df) == 0 && throw(ArgumentError(
        hide_below_tolerance ?
            "no isotopes at or above tolerance=$tolerance up to element_limit=$element_limit" :
            "no isotopes tracked up to element_limit=$element_limit"))

    min_n, max_n = extrema(df.N)
    color_limits = (log10(tolerance), 0.0)

    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1];
            xlabel = "neutron number (A-Z)", ylabel = "proton number (Z)", title,
            aspect = CM.DataAspect(), xgridvisible = false, ygridvisible = false,
            xticks = min_n:max_n, yticks = 0:2:max_z,
            limits = (min_n - 3.0, max_n + 1.0, -0.5, max_z + 1.0))

        draw_isotope_tiles!(ax, df; value_col = :X, color_limits, colormap = ABUNDANCE_COLORMAP, mass_label_size)
        add_element_labels!(ax, df, min_n, max_z; element_label_size)
        CM.Colorbar(fig[1, 2]; colormap = ABUNDANCE_COLORMAP, colorrange = color_limits, label = "log10(X)")
        fig
    end
end

"""
    abundance_chart(run::PPNRun, cycle = :final; kwargs...) -> CM.Figure

Abundance chart for `run` at `cycle` (a cycle number, or `:initial`/`:final`/`:decay`).
"""
abundance_chart(run::PPNRun, cycle = :final; kwargs...) = abundance_chart(abundances(run, cycle); kwargs...)
