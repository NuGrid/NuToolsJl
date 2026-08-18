# xtime_plot.jl — isotope mass fraction vs. time / T9 / density.

const _XTIME_AXIS_LABELS = Dict(
    :time => "time (s)",
    :t9 => "T9 (GK)",
    :rho => "density (g/cm³)",
)

"""
    abundance_vs_time(run::PPNRun, isos; x = :time, logy = true, abundance_floor = 1e-99,
                       title = "Abundance evolution", figure_size = (900, 550)) -> CM.Figure

Mass fraction of each isotope in `isos` (an `Isotope`, an isotope-name string,
or a collection of either) against `:time` (s), `:t9` (GK) or `:rho`
(g cm⁻³), read from `run`'s `x-time.dat`. With `logy = true` (the default),
values are floored at `abundance_floor` before plotting so zero abundances
don't break the log axis.
"""
function abundance_vs_time(run::PPNRun, isos; x::Symbol = :time, logy::Bool = true,
                            abundance_floor::Real = 1e-99,
                            title = "Abundance evolution", figure_size = (900, 550))
    xlabel = get(_XTIME_AXIS_LABELS, x) do
        throw(ArgumentError("x must be :time, :t9 or :rho (got :$x)"))
    end
    isotope_list = isos isa Union{Isotope,AbstractString} ? [isos] : collect(isos)
    isempty(isotope_list) && throw(ArgumentError("no isotopes given"))

    xt = xtime(run)

    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1]; xlabel, ylabel = "mass fraction X", title,
                      yscale = logy ? log10 : identity)
        for (i, iso) in enumerate(isotope_list)
            xv, yv = series(xt, iso; x)
            logy && (yv = max.(yv, abundance_floor))
            isotope = iso isa Isotope ? iso : Isotope(iso)
            CM.lines!(ax, xv, yv; label = isotope_name(isotope),
                       color = NUGRID_PALETTE[mod1(i, length(NUGRID_PALETTE))])
        end
        length(isotope_list) > 1 && CM.axislegend(ax; position = :rb)
        fig
    end
end
