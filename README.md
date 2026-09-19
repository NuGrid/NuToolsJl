# NuToolsJl

A Julia toolkit for analyzing [NuPPN](https://github.com/NuGrid/NuPPN) (nuclear
post-processing) and MESA output: reading run results, nuclear-chart and
sensitivity plots, and orchestration scripts for building rate-factor sweeps.

This repo bundles two packages:

- **[`NuGridJl/`](NuGridJl)** — the active package. Reads completed `nuppn`
  output (abundances, fluxes, `x-time.dat`, the reaction network, the isotope
  database, the input trajectory) into typed, DataFrame-friendly structures,
  and plots them (abundance/flux/ratio/residual nuclear charts, time
  evolution, rate curves, Iliadis (2002)-style rate-sensitivity sweeps).
  `NuGridJl/tools/` holds the run-orchestration scripts (building sweeps,
  launching `ppn.exe`) that write inputs and launch runs — deliberately kept
  out of the package itself, which only ever reads finished output.
- **[`NovaJL/`](NovaJL)** — an earlier, narrower package for nova
  post-processing comparison work (decay solving, comparison plots). Kept
  for existing workflows that depend on it; `NuGridJl` is where new work
  happens.

## Quick start

```bash
git clone https://github.com/NuGrid/NuToolsJl
cd NuToolsJl/NuGridJl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Then run the first demo notebook, which works standalone against fixture
data bundled in `test/data/` — no real `nuppn` install or external data
needed:

```bash
jupyter nbconvert --to notebook --execute --inplace demos/01_reading_ppn_output.ipynb
```

or, from a Julia REPL/script:

```julia
using Pkg
Pkg.activate(joinpath(@__DIR__, "NuGridJl"))
using NuGridJl, DataFrames

run = PPNRun(joinpath(@__DIR__, "NuGridJl", "test", "data", "nuppn_data"))

# abundances at the final cycle, top 10 by mass fraction
sort(DataFrame(abundances(run, :final)), :X; rev = true)[1:10, :]

# the reaction network itself
net = network(run)
(n_isotopes = length(net.isotopes), n_reactions = length(net.reactions))
```

## Learn more

- [`NuGridJl/demos/`](NuGridJl/demos) — five notebooks walking through the
  package in order: reading output, nuclear-chart plotting, comparing runs
  over time, the rate-sensitivity sweep, and reaction reporting. See
  [`NuGridJl/demos/README.md`](NuGridJl/demos/README.md).
- [`NuGridJl/tools/`](NuGridJl/tools) — scripts for building and running
  rate-factor sweeps against a real `nuppn` install. See
  [`NuGridJl/tools/README.md`](NuGridJl/tools/README.md).
- [`NuGridJl/test/runtests.jl`](NuGridJl/test/runtests.jl) — the test suite,
  runnable with `julia --project=NuGridJl -e 'using Pkg; Pkg.test()'`.

## License

BSD-3-Clause (see [LICENSE](LICENSE)), matching the rest of the NuGrid
project.
