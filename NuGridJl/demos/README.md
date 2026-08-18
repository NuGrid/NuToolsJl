# demos/

Jupyter notebooks walking through `NuGridJl`, in order:

1. **`01_reading_ppn_output.ipynb`** — `PPNRun`, abundances, fluxes, the
   `x-time.dat` series, the reaction network, the isotope database.
2. **`02_abundance_and_flux_charts.ipynb`** — the nuclear-chart plots.
3. **`03_ratio_charts_and_time_evolution.ipynb`** — comparing two abundance
   sets, and plotting isotopes over time.
4. **`04_iliadis_sensitivity_sweep.ipynb`** — `PPNSweep` and the Iliadis
   (2002)-style factored-rate sensitivity analysis, plus NPDATA rate curves.
5. **`05_reaction_reporting.ipynb`** — the pre-study (flux) and post-study
   (sensitivity) reaction lists.

All five run standalone against the fixture data bundled in `test/data/` —
no real NuPPN install or external data needed. Each notebook's first code
cell activates the package's own environment
(`Pkg.activate(joinpath(@__DIR__, ".."))`), so it works regardless of which
Julia kernel launched it.

A Monte Carlo ensemble-analysis notebook is planned but not yet written
(deferred alongside `NuGridJl`'s Monte Carlo analysis module and
`tools/build_ensemble.jl`).

## Running these yourself

A Jupyter kernel for this project's environment is registered as
`Julia-NuGridJl 1.12` (via `IJulia.installkernel`). To (re-)execute a
notebook from the command line:

```sh
jupyter nbconvert --to notebook --execute --inplace demos/01_reading_ppn_output.ipynb
```
