module NovaJL

using DataFrames
using Plots
using CSV
using JSON

# -------------------------
# INCLUDE FILES
# -------------------------
include("io.jl")
include("utils.jl")
include("processing.jl")
include("plotting.jl")
include("comparison.jl")
include("decay_solve.jl")

# -------------------------
# EXPORTS
# -------------------------
export read_iso_massf,
       read_trajectory,
       factor_to_folder,
       dots_to_missing!,
       DEFAULT_REACTION_ISOTOPES,
       REACTION_ISOTOPES,
       output_sensitivity_table,
       REACTION_STYLES,
       reaction_styles,
       reaction_style_table,
       reaction_audit,
       factor_audit,
       flux_audit,
       plot_trajectory,
       plot_dens_temp,
       read_rate_curve_data,
       plot_rate_curve,
       abundance_chart,
       flow_chart,
       flow_region_chart,
       analyze_factor,
       nova_case_paths,
       latest_validated_decay_run,
       final_cycle_label,
       load_reaction_plan,
       MODEL_TABLE8_FILES,
       read_iliadis_table8,
       read_iliadis_table3,
       build_ppn_table8,
       compare_to_iliadis,
       classify_goodness,
       score_comparison,
       plot_ppn_vs_iliadis,
       dataframe_to_markdown,
       summary_markdown_table,
       dataframe_to_html,
       RenderedHTML,
       decode_isomer_label,
       isomer_reaction_table,
       run_iliadis_comparison,
       compare_ppn_rate_sets,
       plot_rate_set_comparison,
       run_rate_set_comparison,
       combined_score_table,
       overall_scores,
       isotope_z,
       read_iliadis_final_abundance,
       compare_baseline_to_iliadis,
       score_baseline,
       decay_time_scan_results,
       plot_decay_time_scan,
       plot_baseline_comparison,
       NOVA_HALF_LIVES_S,
       solve_isotope_decay_times,
       plot_decay_time_solve,
       decay_solve_summary,
       solve_isotope_decay_times_analytic,
       compare_decay_solve_methods,
       plot_method_comparison

end
