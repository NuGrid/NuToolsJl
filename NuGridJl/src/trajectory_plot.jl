# trajectory_plot.jl — plots of the thermodynamic trajectory that drives a run.

"""
    plot_trajectory(traj::DataFrame; title = "Trajectory", figure_size = (900, 550)) -> CM.Figure

Temperature (T9, left axis) and density (right axis) vs. time, both on a log
scale, from a trajectory `DataFrame` (see [`read_trajectory`](@ref)/
[`trajectory`](@ref)). Time is left linear since trajectories commonly start
at `t = 0`, which a log axis can't represent.
"""
function plot_trajectory(traj::DataFrame; title = "Trajectory", figure_size = (900, 550))
    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax1 = CM.Axis(fig[1, 1]; xlabel = "time (s)", ylabel = "T9 (GK)", title,
                       yscale = log10, ylabelcolor = :firebrick, yticklabelcolor = :firebrick)
        ax2 = CM.Axis(fig[1, 1]; ylabel = "density (g/cm³)", yscale = log10,
                       yaxisposition = :right, ylabelcolor = :steelblue, yticklabelcolor = :steelblue,
                       xgridvisible = false, ygridvisible = false)
        CM.hidespines!(ax2)
        CM.hidexdecorations!(ax2)
        CM.linkxaxes!(ax1, ax2)
        CM.lines!(ax1, traj.time_s, traj.temperature_T9; color = :firebrick)
        CM.lines!(ax2, traj.time_s, traj.density_cgs; color = :steelblue)
        fig
    end
end

"""
    plot_density_temperature(traj::DataFrame; title = "Density vs Temperature",
                              figure_size = (700, 550)) -> CM.Figure

Log-log density vs. temperature — a quick check of how adiabatic (or not)
the trajectory is, independent of its time axis.
"""
function plot_density_temperature(traj::DataFrame; title = "Density vs Temperature", figure_size = (700, 550))
    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1]; xlabel = "T9 (GK)", ylabel = "density (g/cm³)", title,
                      xscale = log10, yscale = log10)
        CM.lines!(ax, traj.temperature_T9, traj.density_cgs; color = NUGRID_PALETTE[1])
        fig
    end
end

"""
    plot_trajectory(path_or_run; kwargs...) -> CM.Figure
    plot_density_temperature(path_or_run; kwargs...) -> CM.Figure

Convenience overloads taking a `trajectory.input` path or a [`PPNRun`](@ref)
directly.
"""
plot_trajectory(path::AbstractString; kwargs...) = plot_trajectory(read_trajectory(path); kwargs...)
plot_trajectory(run::PPNRun; kwargs...) = plot_trajectory(trajectory(run); kwargs...)
plot_density_temperature(path::AbstractString; kwargs...) = plot_density_temperature(read_trajectory(path); kwargs...)
plot_density_temperature(run::PPNRun; kwargs...) = plot_density_temperature(trajectory(run); kwargs...)
