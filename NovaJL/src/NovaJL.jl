module NovaJL

using DataFrames
using Plots

# -------------------------
# INCLUDE FILES
# -------------------------
include("io.jl")
include("utils.jl")
include("processing.jl")
include("plotting.jl")

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
       analyze_factor

end
