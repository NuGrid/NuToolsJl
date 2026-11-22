#!/usr/bin/env julia
#=
Gather the per-seed output of run_mc into OUT/X.f64, OUT/xr.f64 and OUT/seeds.txt (see `?MC.collect_mc`).

    julia MC/bin/collect_mc.jl OUT
=#
include(joinpath(@__DIR__, "_activate.jl"))

length(ARGS) == 1 || error("usage: julia collect_mc.jl OUT")
collect_mc(ARGS[1])
