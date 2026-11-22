# MC: reaction-rate Monte Carlo sensitivity studies with ppn

Runs a reaction-rate Monte Carlo with nuppn's single-zone `ppn.exe` and finds which rates control the
abundances. Every sample is one ppn run in which every varied rate is multiplied by its uncertainty
factor raised to its own x ~ N(0, 1). nuppn draws the x itself (`factor_uncertainty_mc = .true.`,
`factor_uncertainty_mc_seed = seed`); this package runs the samples, collects them and ranks the rates
by Spearman correlation and mutual information (the KSG estimator of Ward et al. 2025, ApJ 986, 109).

Standard library only (Julia ≥ 1.9). The first run of a `bin/` script sets up the environment by itself.

| | |
|---|---|
| `bin/run_mc.jl` | runs one ppn per seed, `--jobs` at a time, keeping only the final abundances and the drawn x. Restartable, and SLURM-array friendly (`bin/slurm_array.sh`). |
| `bin/collect_mc.jl` | gathers the seeds into `X.f64`, `xr.f64` and `seeds.txt` |
| `bin/analyze_mc.jl` | ratio/isotope uncertainties and the per-target reaction ranking; optionally with Ward et al. (2025)'s CO/ONe values alongside |
| `bin/ward_baseline.jl` | Ward Table 3/4 quantities of single finished ppn runs (no MC), to check a setup |
| `src/` | the `MC` package behind them, for your own analyses (`using MC`) |

## What gets varied
This is set by nuppn (`physics/source/rate_uncertainty.F90`), not by this package:

| rates | factor |
|---|---|
| STARLIB (`starlib_option`) | the file's own f.u.(T) |
| NACRE II, and NACRE I where NACRE II has no rate | the paper limits (x = ±1 lands on a limit) |
| every other nuclear rate with ≥ 2 reactants | `factor_uncertainty_default` (e.g. ×10 / ×0.1 at x = ±1) |
| weak, reverse (these follow their forward rate), photodisintegration, pp/pep, KADoNiS | not varied |

Each sample's `rate_uncertainty_x.txt` lists every varied rate with its x and the source of its
uncertainty. The same seed and the same compiler give the same sample, bit for bit.

## 1. Make a template run
A template is an ordinary ppn run directory with a finished setup: `ppn.exe`, the three `*.input`
files, `isotopedatabase*`, the iniab file and the trajectory. The `NPDATA` link sits in its parent
directory (`CASE/NPDATA`, `CASE/ppn/`). Run it once on its own and check it ends with
`simulation complete`. If you compare with Ward et al. (2025), also run:
```
julia MC/bin/ward_baseline.jl --model co CASE        # or --model one
```
Use the last trajectory point you care about: the analysis takes the last row of `x-time.dat`.

## 2. Smoke test
```
julia MC/bin/run_mc.jl --template CASE/ppn --out smoke20 --seeds 1-20 --jobs 4
julia MC/bin/collect_mc.jl smoke20
julia -t auto MC/bin/analyze_mc.jl smoke20 --model co --iniab CASE/ppn/iniab_co_nova.dat
```
To check a second machine against the first, run the same seeds on both and compare `X.f64`:
```
julia -e 'a=reinterpret(Float64,read("smoke20/X.f64")); b=reinterpret(Float64,read("X_other.f64")); println(maximum(abs.(a./b .- 1)))'
```
The result should be 0, or at round-off level (≲ 1e-6) with a different compiler.

## 3. Full run
Choose `--jobs` to match the cores you want to use. ppn is serial, so one job uses one core.
```
nohup julia MC/bin/run_mc.jl --template CASE/ppn --out mc_10k --seeds 1-10000 --jobs NCORES > mc_10k.log 2>&1 &
```
If the run stops, rerun the same command and it continues where it left off.

With SLURM, edit the `#SBATCH` lines and the paths in `bin/slurm_array.sh`, then `sbatch` it. All
array tasks write into one `OUT`; to resubmit a task that died, use `sbatch --array=<id> slurm_array.sh`.

Failed seeds are listed in `OUT/failed.txt`, with their logs in `OUT/logs/`. A few failures (ppn stopping
when its solver fails) are normal. `--timeout` (default 1800 s) kills runs that hang.

Cost: about 1.3–1.5 s of wall time per sample per core for a nova trajectory on a 2026 desktop, and about
100 kB per sample (the 5,900 drawn x dominate), plus about 15 MB per running worker while it runs.

## 4. Collect and analyze
```
julia MC/bin/collect_mc.jl mc_10k
julia -t auto MC/bin/analyze_mc.jl mc_10k --model co --iniab CASE/ppn/iniab_co_nova.dat --top 10
julia -t auto MC/bin/analyze_mc.jl mc_10k --pairs Na22/Na23,Al26/Al27 --isotopes Na22,Al26,F18
```
- `--model co|one`: Ward's 11 ratios, with their Table 3 (CO, ξ = 9) or Table 4 (ONe, ξ = 17) values and
  Table 5 top reactions alongside. The mixed ratio (Ward Eq. 4: 1 part ejecta to ξ parts initial matter)
  needs `--iniab`.
- `--pairs a/b,...` (implies `--model none`): your own ratios. `--isotopes`: single-isotope targets.
  Add `--iniab FILE --xi XI` for a mixed ratio.
- `--mi-candidates N` (default 40): MI is computed for the N rates with the largest |Spearman| per target.

Written to `OUT` (or `--outdir`):
- `mc_ratios.txt`: per ratio, the median, p16, p84 and factor1 = √(p84/p16), and the same for the
  mixed ratio as factor2; per isotope, the same for X.
- `mc_sensitivity.txt` / `.csv`: per target, the rates ranked by MI (nats) with log10(target), with the
  Spearman correlation and whether the rate's uncertainty is its own or the default factor. Ward et al.
  count MI ≥ 0.3 as important. Below about 1000 samples, MI is biased and noisy.

26Al is the ground state plus the isomer species (`AL26M`) when the network has both.

## In Julia
```julia
using Pkg; Pkg.activate("NuToolsJl/MC"); using MC
ens = MCEnsemble("mc_10k")                       # X (species × n), xr (rates × n), rates, seeds
r = ratio(ens, "Al26", "Al27")                   # per sample; abundance(ens, "Na22") likewise
spread(r)                                        # (median, p16, p84, factor)
rank_reactions(ens, log10.(r))[1:5]              # (index, row, label, reaction, source, mi, spearman)
i = findfirst(q -> q.row == 1287, ens.rates)     # a rate by its networksetup.txt row
ksg_mi(ens.xr[i, :], log10.(r))                  # any rate against any quantity
```
`X.f64` and `xr.f64` are raw Float64 arrays with rows named in `species.txt` and `rates.txt`. In Python,
read them with `numpy.fromfile("X.f64").reshape(n, nsp)`.

## Tests
```
julia --project=MC -t 4 -e 'using Pkg; Pkg.test()'
```
The tests use a stand-in `ppn.exe` (a shell script), so no nuppn build is needed. They cover the
runner (failures, restart, cleanup), collection, the Al26 isomer sum, the mixed ratio, the KSG estimator
against the exact Gaussian MI, the ranking and the reports.
