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
    describe_rate(r::Reaction) -> String

One-line summary of `r`: its label, rate source, printed rate and applied
multiplier — already available from `networksetup.txt` via [`read_network`](@ref),
no NPDATA lookup needed.
"""
describe_rate(r::Reaction) = @sprintf("%-18s source=%-6s rate=%.4e  x%.4g", label(r), r.source, r.rate, r.multiplier)
