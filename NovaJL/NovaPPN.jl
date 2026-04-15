module NovaPPN

using DataFrames
using CairoMakie

# -------------------------
# ORDER MATTERS HERE
# -------------------------
include("utils.jl")        # helpers first
include("processing.jl")   # defines AbundanceSet
include("io.jl")           # parsers
include("plotting.jl")     # uses AbundanceSet

using .Plotting:
       plot_trajectory,
       abundance_chart,
       analyze_factor,
       reaction_stylepoint_table

# -------------------------
# EXPORTS
# -------------------------
export read_initial_abundance,
       read_iso_massf,
       read_trajectory,
       plot_trajectory,
       abundance_chart,
       analyze_factor,
       reaction_stylepoint_table,
       build_abundance_set,
       factor_to_folder,
       dots_to_missing!
end
