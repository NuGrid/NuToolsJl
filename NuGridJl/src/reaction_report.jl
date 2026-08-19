# reaction_report.jl — reaction lists for experimentalists and astrophysicists.
#
# Two reports, matching the two questions a sensitivity study is meant to
# answer: (1) before any sweep, which reactions actually mattered in the
# baseline (from flux), and (2) after a sweep, which of those turned out to be
# sensitive enough to be worth a better rate. Neither cares whether the sweep
# was Iliadis-style factors or a sigma-based STARLIB sweep — both just work
# off `sensitivity_table`'s long-format output.

"""
    flux_reaction_list(run::PPNRun; cycle = :final, threshold = 1e-60) -> DataFrame

The "pre-study" reaction list: every reaction carrying flux `>= threshold` at
`cycle`, ranked by flux, with its label, rate source and printed rate (via
[`label`](@ref)/[`describe_rate`](@ref)'s underlying [`Reaction`](@ref)
fields) — what actually mattered in the baseline, before any reaction gets
factored or sampled.
"""
function flux_reaction_list(run::PPNRun; cycle = :final, threshold::Real = 1e-60)
    fx = fluxes(run, cycle)
    net = network(run)
    reactions_by_index = Dict(r.index => r for r in net.reactions)

    rows = NamedTuple[]
    for row in eachrow(fx)
        row.flux >= threshold || continue
        r = get(reactions_by_index, row.index, nothing)
        r === nothing && continue
        push!(rows, (index = row.index, reaction = label(r), source = r.source, rtype = r.rtype,
                      active = r.active, flux = row.flux, rate = r.rate))
    end
    isempty(rows) && return DataFrame(index = Int[], reaction = String[], source = String[],
                                        rtype = String[], active = Bool[], flux = Float64[], rate = Float64[])
    return sort(DataFrame(rows), :flux; rev = true)
end

"""
    flux_reaction_list(run::PPNRun, cycles; threshold = 1e-60) -> DataFrame

The whole-trajectory form: every reaction carrying flux `>= threshold` at
*any* cycle in `cycles` (pass `run.cycles` for the entire outburst, initial
to final), ranked by its peak flux across those cycles. Same columns as the
single-cycle form, plus `:peak_cycle` (which cycle the peak occurred at) —
`networksetup.txt` lists every reaction the compiled isotope database and
rate libraries make *possible*; this is the much shorter list of what
actually mattered somewhere along this particular trajectory, the list an
astrophysicist would actually want in a hydro sim's reaction network.

Cycles with no `flux_NNNNN.DAT` (commonly the very first cycle — dY/dt isn't
meaningful before any integration step has happened) are skipped rather than
throwing, so `run.cycles` itself is always a safe choice for `cycles`.
"""
function flux_reaction_list(run::PPNRun, cycles; threshold::Real = 1e-60)
    net = network(run)
    reactions_by_index = Dict(r.index => r for r in net.reactions)

    peak = Dict{Int,NamedTuple{(:flux, :cycle),Tuple{Float64,Int}}}()
    for cyc in cycles
        isfile(_flux_path(run, cyc)) || continue
        fx = fluxes(run, cyc)
        for row in eachrow(fx)
            row.flux >= threshold || continue
            prev = get(peak, row.index, nothing)
            (prev === nothing || row.flux > prev.flux) && (peak[row.index] = (flux = row.flux, cycle = Int(cyc)))
        end
    end

    rows = NamedTuple[]
    for (idx, best) in peak
        r = get(reactions_by_index, idx, nothing)
        r === nothing && continue
        push!(rows, (index = idx, reaction = label(r), source = r.source, rtype = r.rtype,
                      active = r.active, flux = best.flux, peak_cycle = best.cycle, rate = r.rate))
    end
    isempty(rows) && return DataFrame(index = Int[], reaction = String[], source = String[], rtype = String[],
                                        active = Bool[], flux = Float64[], peak_cycle = Int[], rate = Float64[])
    return sort(DataFrame(rows), :flux; rev = true)
end

"""
    sensitivity_reaction_report(table::DataFrame; threshold = 0.1) -> DataFrame

The "post-study" reaction list: collapses a long-format
[`sensitivity_table`](@ref) (one row per reaction×isotope×factor) to one row
per reaction — how many isotopes/factors it was tested against, its biggest
observed `|log10(ratio)|` swing (`max_abs_log_ratio`) and which
isotope/factor produced it, and whether that crosses `threshold` (`sensitive`).
Ranked most sensitive first — the "what to remeasure" list.
"""
function sensitivity_reaction_report(table::DataFrame; threshold::Real = 0.1)
    function summarize(g)
        log_ratios = log10.(g.ratio)
        mask = isfinite.(log_ratios)
        if !any(mask)
            return (n_factors = length(unique(g.factor)), n_isotopes = length(unique(g.isotope)),
                     max_abs_log_ratio = NaN, worst_isotope = missing, worst_factor = missing,
                     sensitive = false)
        end
        abs_lr = abs.(log_ratios[mask])
        i = argmax(abs_lr)
        (
            n_factors = length(unique(g.factor)),
            n_isotopes = length(unique(g.isotope)),
            max_abs_log_ratio = abs_lr[i],
            worst_isotope = g.isotope[mask][i],
            worst_factor = g.factor[mask][i],
            sensitive = abs_lr[i] >= threshold,
        )
    end
    report = combine(groupby(table, :reaction), summarize)
    return sort(report, :max_abs_log_ratio; rev = true)
end
