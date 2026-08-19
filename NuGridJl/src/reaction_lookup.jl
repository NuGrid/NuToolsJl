# reaction_lookup.jl — structural lookups over an already-parsed Network.

"""
    reactions_for_isotope(net::Network, iso::Isotope; active_only = true) -> Vector{Reaction}

Every reaction in `net` where `iso` appears as a reactant or product — plain
structural participation, not ranked by flux or significance. Pair with
[`changed_isotopes`](@ref) to go from "this isotope moved" to "these are the
candidate reactions to factor."
"""
function reactions_for_isotope(net::Network, iso::Isotope; active_only::Bool = true)
    filter(net.reactions) do r
        (!active_only || r.active) && (iso in r.reactants || iso in r.products)
    end
end

"""
    reactions_for_isotope(run::PPNRun, iso::Isotope, cycle_or_cycles; active_only = true,
                           threshold = 1e-60) -> DataFrame

The flux-aware companion to the `Network`-based method: structural
candidates (`iso` a reactant or product) that also carry flux `>= threshold`
somewhere in `cycle_or_cycles` — a single cycle for "at this point in the
run," or a vector of cycles (pass `run.cycles` for the whole trajectory,
initial to final) for "anywhere in this stellar environment's outburst."
Ranked by flux via [`flux_reaction_list`](@ref), most important first.
"""
function reactions_for_isotope(run::PPNRun, iso::Isotope, cycle_or_cycles; active_only::Bool = true,
                                threshold::Real = 1e-60)
    candidates = Set(r.index for r in reactions_for_isotope(network(run), iso; active_only))
    full = cycle_or_cycles isa AbstractVector ?
        flux_reaction_list(run, cycle_or_cycles; threshold) :
        flux_reaction_list(run; cycle = cycle_or_cycles, threshold)
    return filter(:index => in(candidates), full)
end

"""
    describe_rate(r::Reaction) -> String

One-line summary of `r`: its label, rate source, printed rate and applied
multiplier — already available from `networksetup.txt` via [`read_network`](@ref),
no NPDATA lookup needed.
"""
describe_rate(r::Reaction) = @sprintf("%-18s source=%-6s rate=%.4e  x%.4g", label(r), r.source, r.rate, r.multiplier)
