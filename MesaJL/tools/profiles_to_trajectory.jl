#!/usr/bin/env julia
# profiles_to_trajectory.jl — turn a MESA LOGS/ directory into a single-zone
# ppn trajectory.input file, following the same (t, T, rho) format used by
# the existing NovaSensitivityStudy/single-zone/nova_cases/*/ppn/trajectory*.input
# runs (see nuppn/utils/source/trajectories.F90 for the format nuppn's `ppn`
# frame actually parses).
#
# Usage:
#   julia profiles_to_trajectory.jl [logs_dir] [output_path] [target_mass]
#
# All arguments are optional and default to the PostTNR co_nova_1.1_fiducial
# run. target_mass (in Msun) picks which Lagrangian zone to track; if
# omitted, it's chosen automatically as the mass coordinate that reaches
# peak temperature during the TNR (see extract_zone_trajectory).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MesaJL
using DataFrames

const POSTTNR = joinpath(@__DIR__, "..", "..", "..")

logs_dir = length(ARGS) >= 1 ? ARGS[1] :
    joinpath(POSTTNR, "mesa", "co_nova_1.1_fiducial", "LOGS")
output_path = length(ARGS) >= 2 ? ARGS[2] :
    joinpath(POSTTNR, "trajectories", "trajectory_co_nova_1.1_fiducial.input")
target_mass = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : nothing

traj = extract_zone_trajectory(logs_dir; target_mass = target_mass)

mkpath(dirname(output_path))
write_ppn_trajectory(traj, output_path)

age_days = (traj.age_yr[end] - traj.age_yr[1]) * 365.25
println("Tracked zone mass   = ", traj.mass_Msun[1], " Msun",
        " (drift over run: ", extrema(traj.mass_Msun), ")")
println("Timesteps written   = ", nrow(traj))
println("Age span covered    = ", round(age_days, digits = 3), " days")
println("Peak temperature    = ", maximum(traj.temperature_K), " K",
        " (T9 = ", round(maximum(traj.temperature_K) / 1e9, digits = 3), ")")
println("Peak density        = ", maximum(traj.density_cgs), " g/cc")
println("Wrote: ", output_path)
