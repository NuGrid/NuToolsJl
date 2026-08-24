# reaction_lookup.jl — structural lookups over an already-parsed Network.

"""
    reaction_by_index(net::Network, index::Integer) -> Reaction

The reaction at `index` — the direct answer to "I have a raw reaction index
from eyeballing a `flux_*.DAT`/`iso_massf*.DAT` file (they only carry `Z`/`A`
and the index, no source label — see [`read_fluxes`](@ref)), what actually
is it and which rate library (`.source`) supplied it?" Throws if `index`
isn't in `net`. For looking many reactions up at once, [`flux_reaction_list`](@ref)
already joins flux data against the network and includes `source` as a
column — reach for this only for a one-off index you're staring at by hand.
"""
function reaction_by_index(net::Network, index::Integer)
    i = findfirst(r -> r.index == index, net.reactions)
    i === nothing && throw(ArgumentError("no reaction with index $index in this network"))
    return net.reactions[i]
end

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
Ranked by flux via [`flux_reaction_list`](@ref), most important first. A
thin convenience wrapper over `flux_reaction_list`'s own `isotope` keyword —
equivalent to calling it directly with `isotope = iso`.
"""
function reactions_for_isotope(run::PPNRun, iso::Isotope, cycle_or_cycles; active_only::Bool = true,
                                threshold::Real = 1e-60)
    return cycle_or_cycles isa AbstractVector ?
        flux_reaction_list(run, cycle_or_cycles; threshold, isotope = iso, active_only) :
        flux_reaction_list(run; cycle = cycle_or_cycles, threshold, isotope = iso, active_only)
end

"""
    describe_rate(r::Reaction) -> String

One-line summary of `r`: its label, rate source, printed rate and applied
multiplier — already available from `networksetup.txt` via [`read_network`](@ref),
no NPDATA lookup needed. The printed rate is the fixed build-time snapshot
`r.rate` (see [`Reaction`](@ref) — evaluated once at T9≈0, not the rate at
any particular cycle); labeled `rate(T9≈0)=` here so it can't be mistaken for
one. For the real, temperature-dependent curve use [`rate_curve`](@ref).
"""
describe_rate(r::Reaction) = @sprintf("%-18s source=%-6s rate(T9≈0)=%.4e  x%.4g", label(r), r.source, r.rate, r.multiplier)

"""
    self_loop_reactions(net::Network; active_only = true) -> Vector{Reaction}

Reactions where the same isotope appears on both sides — a reactant that is
also a product. This is essentially always a `network_boundaries.F90`/
`natashamcloane` redirect artifact, not real physics: the reaction's true
product falls outside the tracked isotope list,
so the isobar-walk that's supposed to land on the nearest *tracked* neighbor
instead lands back on the reactant itself (e.g. a restricted network missing
Li-7 turns `7Be(n,p)7Li` into a recorded `7Be(n,p)7Be`).

Confirmed against `jac_rhs.F90`'s `calculate_dxdt`: a reactant/product-1
self-loop like this has *zero* net effect on the shared isotope's own
abundance (the consumption and production terms land in the same `dxdt` slot
with equal and opposite sign and cancel exactly) — but it silently swallows
whatever real destruction/production channel the redirect was standing in
for, and the *other* reactant/product in the reaction is still consumed/
produced for real (in the Be-7 example, a neutron genuinely disappears and a
proton genuinely appears, with no corresponding Li-7 ever created). Worth
running as a one-off audit of a restricted isotope list, or after any
[`reactions_for_isotope`](@ref)/[`flux_reaction_list`](@ref) result that looks
structurally odd.
"""
function self_loop_reactions(net::Network; active_only::Bool = true)
    filter(net.reactions) do r
        (!active_only || r.active) && !isempty(intersect(r.reactants, r.products))
    end
end
