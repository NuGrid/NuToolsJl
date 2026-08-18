"""
    NuGridJl

A typed, CairoMakie-based successor to NuGridPy's nucleosynthesis analysis
side, focused on single-zone `nuppn` output and reaction sensitivity studies.

`using NuGridJl` gives you a flat, curated API: a `PPNRun` directory handle,
value types (`Isotope`, `Abundances`, `Reaction`), latex-themed abundance/flux/
ratio/rate charts, an Iliadis (2002)-style factored-rate sensitivity sweep
built on `PPNSweep`, and the readers behind all of it. It only reads
already-completed run output — it never launches `ppn.exe` itself; the
scripts that build sweep/Monte Carlo run directories live in `tools/`,
outside this package.
"""
module NuGridJl

using Printf
using Statistics
using Random
using DelimitedFiles
using DataFrames
using CSV
using JSON
using LaTeXStrings
import CairoMakie as CM

include("elements.jl")
include("theme.jl")
include("io_ppn.jl")
include("io_network.jl")
include("io_isotopedatabase.jl")
include("run.jl")
include("runset.jl")
include("chart_primitives.jl")
include("abundance_chart.jl")
include("flux_chart.jl")
include("ratio_chart.jl")
include("xtime_plot.jl")
include("reaction_lookup.jl")
include("npdata.jl")
include("tables.jl")
include("sensitivity_iliadis.jl")
include("reaction_report.jl")

export
    # elements
    Isotope, proton_number, mass_number, neutron_number,
    element_symbol, isotope_name, latex_label, is_stable, is_magic,
    ELEMENT_SYMBOLS, MAGIC_NUMBERS, STABLE_ISOTOPES,
    # theme
    nugrid_theme, set_nugrid_theme!, reset_theme!, with_nugrid_theme, savefig,
    NUGRID_PALETTE, ABUNDANCE_COLORMAP, RATIO_COLORMAP, FLUX_COLORMAP,
    # io + types
    Abundances, XTime, Reaction, Network,
    read_abundances, read_fluxes, read_xtime, read_iniab, read_trajectory,
    read_network, read_namelist,
    isotopes, mass_fraction, series, label,
    # isotope database
    IsotopeDatabaseEntry, read_isotopedatabase,
    # run
    PPNRun, abundances, fluxes, xtime, network, inputs, isotopedatabase,
    # sweep / ensemble
    FactoredRun, PPNSweep, MCSample, PPNEnsemble,
    reactions, factors, sweep_run, sample,
    # charts
    abundance_chart, flux_chart, ratio_chart, changed_isotopes, abundance_vs_time,
    # network lookups
    reactions_for_isotope, describe_rate,
    # NPDATA / REACLIB rate curves
    ReaclibBlock, read_reaclib_blocks, reaclib_species, reaclib_rate,
    rate_curve_temperature_grid, rate_curve, rate_plot,
    # table export
    RenderedHTML, dataframe_to_markdown, dataframe_to_html, save_table,
    # Iliadis-style sensitivity
    sensitivity, sensitivity_table, iliadis_table, rank_reactions, sensitivity_plot,
    # reaction reports
    flux_reaction_list, sensitivity_reaction_report

end # module
