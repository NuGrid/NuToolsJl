#!/usr/bin/env julia
#=
Run a ppn reaction-rate Monte Carlo (see `?MC.run_mc`).

    julia MC/bin/run_mc.jl --template CASE/ppn --out OUT --seeds 1-10000 --jobs 32 [--timeout 1800]

CASE/ppn: a finished-setup ppn run directory whose parent holds the NPDATA link. --seeds takes "a-b",
"a,b,c" or a mix. Rerunning the same command continues where it stopped; several processes (SLURM array
tasks, see slurm_array.sh) can share one OUT. Then: julia MC/bin/collect_mc.jl OUT
=#
include(joinpath(@__DIR__, "_activate.jl"))

o, pos = cli_args(ARGS, Dict("jobs" => string(Sys.CPU_THREADS), "timeout" => "1800"))
isempty(pos) || error("unexpected argument $(pos[1])")
for k in ("template", "out", "seeds")
    haskey(o, k) || error("missing --$k (see the header of run_mc.jl)")
end
run_mc(o["template"], o["out"], parse_seeds(o["seeds"]); jobs=parse(Int, o["jobs"]),
       timeout=parse(Float64, o["timeout"]))
