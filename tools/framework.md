# Nova Sensitivity Study — Workflow

This framework runs a PPN (Post-Processing Network) sensitivity study on a nova model: each reaction rate in `config/reaction_plan.json` is perturbed by a set of multiplicative factors, and the resulting isotopic yields are compared across rate sets and against the Iliadis et al. 2002 reference sensitivity study.

All tools live in `tools/`. Notebooks live in `nova_cases/<case>/analysis/`. All Julia scripts accept `--help` for full option documentation.

---

## Rate sets

| Label  | `starlib_option` | Source file                              | Notes                                      |
|--------|-----------------|------------------------------------------|--------------------------------------------|
| REACLIB | — (baseline)   | built into PPN                           | Iliadis-era rates; used by `run_ppn_sweep.jl` |
| STL01  | 1               | `starlib_mc10_mc13_082022.txt`           | MC10 + MC13 combined, 59 reactions         |
| STL02  | 2               | `starlib_etr25_2025.txt`                 | ETR25 2025 evaluation, 76 reactions        |

---

## Step 0 — Build PPN (once, or after any Fortran change)

```bash
cd nova_cases/<case>/ppn
make distclean && make
```

Required after any change to `physics/source/` (e.g. `ppn_physics.F90`, `array_sizes.F90`).

---

## Step 1 — Network setup (only if rates or network changed)

```bash
julia tools/setup_network.jl --nova <case>
```

Edits `ppn/networksetup.txt` according to `config/network_edits.json`. Skip if nothing changed.

---

## Step 2 — Nucleosynthesis sweep

### Old REACLIB rates

```bash
julia tools/run_ppn_sweep.jl --nova <case> --jobs 8
```

Outputs: `nova_cases/<case>/runs/`

### STARLIB rate sets (run once per option)

```bash
julia tools/run_ppn_sweep_starlib.jl --nova <case> --starlib-option 1 --runs-name runs_star_opt1 --jobs 8
julia tools/run_ppn_sweep_starlib.jl --nova <case> --starlib-option 2 --runs-name runs_star_opt2 --jobs 8
```

Outputs: `nova_cases/<case>/runs_star_opt{1,2}/`

### How the sweep works (three phases)

Both sweep scripts follow the same three-phase logic:

1. **Phase 1 — baseline**: copies `ppn/` to `runs/baseline/`, (for STARLIB: writes `starlib_option=N` into `ppn_physics.input`), runs `ppn.exe`. This generates `runs/baseline/networksetup.txt` with the actual reaction ordering for this configuration.
2. **Phase 2 — index resolution**: reads `runs/baseline/networksetup.txt`, resolves each reaction's network index from that file, and builds all factored run directories with the correct indices.
3. **Phase 3 — parallel runs**: runs all factored directories in parallel (baseline excluded).

> **Why this matters for STL01/STL02:** STARLIB inserts additional reactions into the network, shifting all subsequent indices relative to old REACLIB. Phase 1 must run first so that Phase 2 reads the correct post-insertion ordering for each option.

### Reaction index resolution notes

- A reaction may appear twice in `networksetup.txt` — once active (`T`) and once inactive (`F`). VITAL entries are hand-coded rates that are superseded by a newer evaluation set to `T`. The sweep always uses the active (`T`) row and displays all candidates with their status.
- If a reaction's configured `index` in `reaction_plan.json` is not found in the current network (e.g. the index came from a different rate set), the sweep warns and auto-selects the active row.
- If a reaction has no active row (all `F`), the sweep proceeds with a warning — the perturbation will likely have no effect.
- If a reaction is completely absent from the network, it is skipped with a warning.

---

## Step 3 — Decay time calibration

```bash
julia tools/decay_time_scan.jl --nova <case>
```

Runs real PPN decay from the baseline across a grid of decay times. Results go to `analysis/results/baseline_decay_checker/`.

Then open and run:

```
nova_cases/<case>/analysis/00_decay_time_calibration.ipynb
```

This notebook (Julia kernel, NovaJL) compares each scanned decay time against Iliadis 2002 Table 4 (JCH1 reference abundances) and scores them. The best decay time is written into the scan CSV and used automatically by the analysis notebooks downstream.

---

## Step 4 — Decay sweep

Apply the chosen decay time to every forward run:

```bash
# Old REACLIB rates
julia tools/decay_ppn_sweep.jl --nova <case> \
  --runs-name runs --decay-runs-name decay_runs --jobs 8 --decay-time <chosen>

# STARLIB options
julia tools/decay_ppn_sweep.jl --nova <case> \
  --runs-name runs_star_opt1 --decay-runs-name decay_runs_star_opt1_8hrs --jobs 8 --decay-time <chosen>
julia tools/decay_ppn_sweep.jl --nova <case> \
  --runs-name runs_star_opt2 --decay-runs-name decay_runs_star_opt2_8hrs --jobs 8 --decay-time <chosen>
julia tools/decay_ppn_sweep.jl --nova <case> \
  --runs-name runs_star_opt3 --decay-runs-name decay_runs_star_opt3_8hrs --jobs 8 --decay-time <chosen>
```

Outputs mirror the forward runs structure: `decay_runs_*/baseline/iso_massfdecay.DAT` and one per reaction×factor.

---

## Step 5 — Analysis notebooks

Open and run in order from `nova_cases/<case>/analysis/`. All use the Julia kernel and the `NovaJL` package at `analysis/NovaJL/`.

| Notebook | What it does |
|----------|-------------|
| `00_decay_time_calibration.ipynb` | Calibrates decay time vs. Iliadis 2002 Table 4 |
| `01_Iliadis02_table_8_rec_eval.ipynb` | Compares PPN (old REACLIB) sensitivity against Iliadis 2002 Table 8 |
| `02_PPN_oldrates_vs_starlib_opt1.ipynb` | Old REACLIB rates vs. STL01 (MC10/MC13) |
| `03_STARLIB_opt1_vs_opt2.ipynb` | STL01 (MC10/MC13) vs. STL02 (STARLIB v6.10) |
| `04_STARLIB_opt2_vs_opt3.ipynb` | STL02 (STARLIB v6.10) vs. STL03 (ETR25 2025) |
| `05_rate_evolution_summary.ipynb` | Full 25-year rate evolution summary across all rate sets |

Scored comparison outputs (`.md` score tables, `.csv` detail files) are written to `analysis/results/`.

---

## Configuration

### `config/reaction_plan.json`

Defines which reactions to perturb and with what factors. Each entry:

```json
{
  "name": "23Na_pg_24Mg",
  "article_reaction": "23Na(p,γ)24Mg",
  "isotopes": ["24Mg", "23Na"],
  "factors": [0.01, 0.1, 0.5, 2.0, 10.0, 100.0],
  "index": 392
}
```

- `index` is the reaction's line number in `networksetup.txt`. For rate sets that shift indices (e.g. STL02), the sweep's Phase 1/2 resolves the correct runtime index automatically — the configured `index` is used as a fallback hint only.
- `reverse_index` (optional): include for alpha-transfer reactions to perturb forward and reverse simultaneously.

---

## NovaJL package

Located at `analysis/NovaJL/src/`. Key modules:

| File | Contents |
|------|----------|
| `io.jl` | Reading `iso_massf*.DAT`, `networksetup.txt`, Iliadis data tables |
| `processing.jl` | Sensitivity table assembly, factor/reaction audits |
| `comparison.jl` | Iliadis scoring (`run_iliadis_comparison`), rate-set scoring (`run_rate_set_comparison`), `build_ppn_table8` |
| `plotting.jl` | Abundance charts, flow tiles, sensitivity scatter plots (Makie + Plots) |
| `decay_solve.jl` | Lightweight post-decay solver |
| `utils.jl` | Shared helpers |

---

## Known caveats

- **30P(p,γ)31S Table 8**: PPN vs Iliadis shows a large ratio (~2 dex) for some isotopes for this reaction. This is an anomaly in the Iliadis Table 8 entry itself, not a code bug. Documented in notebook 01.
- **STARLIB v6.10 index shifting**: STL02 inserts ~42 reactions before many existing ones. Always re-run from scratch when switching between rate sets; never copy a `reaction_plan.json` index from one rate set and expect it to be valid for another.
- **Julia cold-start latency**: first `julia` invocation takes ~30–60 s for precompilation. Subsequent calls in the same session are fast.
