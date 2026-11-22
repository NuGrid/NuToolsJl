"""
    MC

Reaction-rate Monte Carlo sensitivity studies with nuppn's single-zone `ppn.exe`:

- [`run_mc`](@ref): run one ppn per seed (`factor_uncertainty_mc_seed`), many at a time, keeping the
  final abundances and the drawn rate factors; restartable, SLURM-array friendly;
- [`collect_mc`](@ref) and [`MCEnsemble`](@ref): gather and load the samples;
- [`abundance`](@ref), [`ratio`](@ref), [`spread`](@ref): the quantities and their uncertainty factors;
- [`rank_reactions`](@ref): Spearman + KSG mutual-information ranking of the varied rates;
- [`sensitivity_report`](@ref), [`ward_report`](@ref), [`ward_baseline`](@ref): reports, with the
  Ward et al. (2025) nova values for comparison.

Standard library only. Command-line front ends are in `MC/bin/`.
"""
module MC

using Printf, Random, Statistics

export fnum, iso_name, read_iniab, run_iniab, last_xtime_row, final_abundances, read_rate_x,
       parse_seeds, physics_with_seed, run_seed, run_mc,
       collect_mc, MCEnsemble, RateInfo, nsamples, abundance, ratio, uncertainty_source,
       ranks, spearman, spearman_all, ksg_mi, spread,
       rank_reactions, sensitivity_report,
       WARD_PAIRS, WARD_XI, WARD2025, WARD2025_TOP, ward_report, ward_baseline

include("io.jl")
include("runner.jl")
include("ensemble.jl")
include("stats.jl")
include("ward2025.jl")
include("sensitivity.jl")

end
