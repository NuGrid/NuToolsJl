#!/usr/bin/env julia
#=
Sensitivity analysis of a collected Monte Carlo (see `?MC.sensitivity_report`, `?MC.ward_report`).

    julia -t auto MC/bin/analyze_mc.jl OUT [--model co|one|none] [--iniab FILE] [--xi XI]
                                         [--pairs C12/C13,Na22/Na23] [--isotopes Na22,Al26]
                                         [--top 10] [--mi-candidates 40] [--outdir DIR]

--model co|one: Ward et al. (2025) ratios with their Table 3 (CO, xi = 9) or Table 4 (ONe, xi = 17)
values and Table 5 top reactions alongside; needs --iniab (the runs' initial abundances) for the mixed
ratio. --model none (the default when --pairs is given): only the ratios in --pairs, mixed ratio only if
--iniab and --xi are given. --isotopes adds single-isotope targets (log10 X). Writes mc_ratios.txt,
mc_sensitivity.txt and mc_sensitivity.csv in --outdir (default OUT). -t auto threads the MI/Spearman loops.
=#
include(joinpath(@__DIR__, "_activate.jl"))

o, pos = cli_args(ARGS, Dict("top" => "10", "mi-candidates" => "40"))
length(pos) == 1 || error("usage: julia -t auto analyze_mc.jl OUT [--model co|one|none] [--iniab FILE] ...")
model = Symbol(get(o, "model", haskey(o, "pairs") ? "none" : "co"))
ens = MCEnsemble(pos[1])
println(ens)
kw = (top=parse(Int, o["top"]), candidates=parse(Int, o["mi-candidates"]), outdir=get(o, "outdir", pos[1]),
      isotopes=haskey(o, "isotopes") ? String.(split(o["isotopes"], ',')) : String[])
if model in (:co, :one)
    haskey(o, "iniab") || error("--model $model needs --iniab")
    haskey(o, "pairs") && error("--pairs is for --model none")
    ward_report(ens; model, iniab=o["iniab"], xi=haskey(o, "xi") ? parse(Float64, o["xi"]) : WARD_XI[model], kw...)
elseif model == :none
    pairs = haskey(o, "pairs") ? [Tuple(String.(split(p, '/'))) for p in split(o["pairs"], ',')] : WARD_PAIRS
    X0 = haskey(o, "iniab") ? read_iniab(o["iniab"]) : nothing
    sensitivity_report(ens; pairs, X0, xi=haskey(o, "xi") ? parse(Float64, o["xi"]) : 0.0, kw...)
else
    error("--model co|one|none")
end
