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
       REACTION_STYLES,
       reaction_styles,
       plot_trajectory,
       plot_dens_temp,
       analyze_factor

end
