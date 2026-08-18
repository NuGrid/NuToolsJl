# theme.jl — a shared CairoMakie look, applied locally.
#
# Unlike the old NovaJL code, importing NuGridJl must never mutate global Makie
# state.  We expose a theme *value* and a `with_nugrid_theme` wrapper; every plot
# builds its figure inside that wrapper, so a user's own `set_theme!` survives.

"Qualitative colour cycle (colour-blind friendly), reused across line/scatter plots."
const NUGRID_PALETTE = [
    CM.RGBf(0.00, 0.45, 0.70),   # blue
    CM.RGBf(0.84, 0.37, 0.00),   # vermillion
    CM.RGBf(0.00, 0.62, 0.45),   # green
    CM.RGBf(0.80, 0.47, 0.65),   # reddish purple
    CM.RGBf(0.34, 0.71, 0.91),   # sky blue
    CM.RGBf(0.94, 0.89, 0.26),   # yellow
    CM.RGBf(0.55, 0.34, 0.29),   # brown
    CM.RGBf(0.35, 0.35, 0.35),   # gray
]

"Sequential colormap for abundance mass fractions (white → deep red)."
const ABUNDANCE_COLORMAP = CM.cgrad([:white, CM.RGBf(0.70, 0.0, 0.0)])

"Diverging colormap for log ratios, centred at 0 (blue → white → red)."
const RATIO_COLORMAP = CM.cgrad([:steelblue, :white, :firebrick])

"Sequential colormap for reaction fluxes."
const FLUX_COLORMAP = :inferno

"""
    nugrid_theme(; fontsize = 15, resolution = (900, 650))

The house style: `theme_latexfonts()` plus tuned fonts, gridlines and a
qualitative colour cycle.  Returns a `Makie.Theme` you can pass to `with_theme`
or [`set_nugrid_theme!`](@ref); it is never applied on import.
"""
function nugrid_theme(; fontsize::Integer = 15, resolution = (900, 650))
    CM.merge(
        CM.Theme(
            fontsize = fontsize,
            size = resolution,
            figure_padding = 12,
            palette = (color = NUGRID_PALETTE,),
            Axis = (
                xgridvisible = true,
                ygridvisible = true,
                xgridcolor = (:gray, 0.25),
                ygridcolor = (:gray, 0.25),
                xtickalign = 1,
                ytickalign = 1,
                titlefont = :bold,
            ),
            Legend = (framevisible = false, patchsize = (18, 14)),
            Colorbar = (tickalign = 1,),
        ),
        CM.theme_latexfonts(),
    )
end

"""
    set_nugrid_theme!(; kwargs...)

Opt in to the house style globally for the session (calls `set_theme!`).  Plots
already look right without this; it only matters for figures built by hand.
"""
set_nugrid_theme!(; kwargs...) = CM.set_theme!(nugrid_theme(; kwargs...))

"Reset Makie back to its default theme."
reset_theme!() = CM.set_theme!()

"""
    with_nugrid_theme(f; kwargs...)

Run `f()` with the house theme active, restoring the previous theme afterwards.
Every plotting function in NuGridJl wraps its body in this.
"""
with_nugrid_theme(f; kwargs...) = CM.with_theme(f, nugrid_theme(; kwargs...))

"""
    savefig(fig, path; kwargs...)

Save `fig` to `path`, creating parent directories as needed.  Thin wrapper over
`CairoMakie.save` so notebooks and scripts have one obvious verb.
"""
function savefig(fig, path::AbstractString; kwargs...)
    mkpath(dirname(abspath(path)))
    CM.save(path, fig; kwargs...)
    return path
end
