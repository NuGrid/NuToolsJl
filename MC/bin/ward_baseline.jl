#!/usr/bin/env julia
#=
Ward et al. (2025) Table 3/4 quantities of single finished ppn runs (see `?MC.ward_baseline`).

    julia MC/bin/ward_baseline.jl [--model co|one] RUN [RUN ...]

RUN: a ppn run directory (with x-time.dat) or a case directory with a ppn/ subdirectory.
=#
include(joinpath(@__DIR__, "_activate.jl"))

o, pos = cli_args(ARGS, Dict("model" => "co"))
isempty(pos) && error("usage: julia ward_baseline.jl [--model co|one] RUN [RUN ...]")
ward_baseline(pos; model=Symbol(o["model"]))
