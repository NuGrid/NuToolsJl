# tools/

Run-orchestration scripts. These are deliberately **not** part of the
`NuGridJl` package — `NuGridJl` (`src/`) only ever reads already-completed
`nuppn` output; anything that writes `ppn_physics.input` or launches
`ppn.exe` lives here instead.

Run with the `NuGridJl` package's own environment, e.g. from this directory's
parent:

```sh
julia --project=. tools/build_sweep.jl <template_dir> <reaction_plan.json> <out_dir> [--jobs N] [--dry-run]
```

## `build_sweep.jl`

Builds an Iliadis (2002)-style one-reaction-at-a-time rate sweep: a
`baseline/` copy of `template_dir` (run once so its `networksetup.txt` is the
authoritative reaction ordering — rate-set changes can shift indices), then a
`<reaction>/fact_<factor>/` directory per entry in a reaction-plan JSON file,
with that reaction's `rate_index`/`rate_factor` written into
`ppn_physics.input`. Read the result back with `NuGridJl.PPNSweep(out_dir)`.

This ports the reaction-resolution logic (`resolve_reaction_index`,
`matching_rows_for_reaction`, `validate_reverse_index`, `copy_ppn`) from
`NovaSensitivityStudy/single-zone/tools/NovaRunTools.jl` onto NuGridJl's
typed `Network`/`Reaction`/`Isotope`, generalized off that module's
`nova_cases/`-specific paths. Deliberately **not** ported:
`source_priority`/`enforce_activation_policy`/`network_edits.json` — that's
network *curation* (making sure each physical reaction has exactly one
active row before a sweep is ever built), a separate concern from
sweep-building that `resolve_reaction_index` doesn't depend on either.
Worth adding as its own `setup_network.jl`-equivalent tool if it turns out
to be needed here too.

Reaction-plan schema — a JSON array (or `{"reactions": [...]}`) of:

```json
{"name": "13N_pg_14O", "factors": [0.5, 2.0]}
```

optionally with `"index"` (force a specific network row instead of
resolving by name — needed when a reaction has more than one active row and
`--priorities` doesn't disambiguate it) and/or `"reverse_index"` (also
factor the reverse reaction by the same amount, for alpha-transfer
reactions — validated to actually be the reverse, not just assumed).

`name` is `<target><channel><product>`, matching the convention already used
throughout `NovaSensitivityStudy`'s `reaction_plan.json` files and sweep
directory names (e.g. `"3He_ag_7Be"`, `"25Mg_pg_26Alg"` for isomer-tagged
products). Channel codes: `pg`, `pa`, `pn`, `pp`, `ag`, `an`, `ap`, `ng`,
`na`, `np`.

`--priorities SRC,...` is a source-label tie-break list, tried in order,
used only on reactions where more than one row is active. If it doesn't
narrow a reaction to one row, `build_sweep` throws rather than guessing —
use an explicit `"index"` for that reaction instead.

Use `--dry-run` to build the directory tree and resolve every reaction name
against the network (throws if any name doesn't match, or matches
ambiguously) without actually launching `ppn.exe` — a cheap way to check a
reaction plan before committing real compute to it.

## `run_parallel.jl`

Shared helper: run a command in N directories concurrently, up to a job
limit, logging each to `<dir>/run.log`. `include`d by `build_sweep.jl`.

## Monte Carlo ensemble building

Not yet implemented here — the STARLIB Monte Carlo sampling itself already
exists natively in `nuppn`'s Fortran physics engine (`starlib_mc_rfac`/
`etr25_mc_rfac` + seed knobs in `ppn_physics.input`), so an ensemble builder
just needs to vary those seeds across N run directories. Coming alongside the
Monte Carlo analysis side of `NuGridJl` (deferred for now).
