# sensitivity_iliadis.jl — Iliadis (2002)-style one-reaction-at-a-time
# factor-sweep sensitivity.
#
# The whole algorithm is a post-hoc ratio: each tested rate factor is already
# a completed, independent single-zone integration (built by tools/, read
# here via PPNSweep) — sensitivity is just `X(factor) / X(baseline)`. No
# perturbation theory, no correlation between reactions.

"""
    sensitivity(sweep::PPNSweep, reaction, iso; cycle = :final) -> Dict{Float64,Float64}

For every tested `factor` of `reaction` in `sweep`, the ratio
`X(iso, factor) / X(iso, baseline)` at `cycle`. A zero baseline abundance
gives `Inf` (grew from nothing) or `NaN` (stayed at nothing).
"""
function sensitivity(sweep::PPNSweep, reaction::AbstractString, iso; cycle = :final)
    isotope = iso isa Isotope ? iso : Isotope(iso)
    baseline_x = abundances(sweep.baseline, cycle)[isotope]
    ratios = Dict{Float64,Float64}()
    for f in factors(sweep, reaction)
        x = abundances(sweep_run(sweep, reaction, f), cycle)[isotope]
        ratios[f] = baseline_x == 0 ? (x == 0 ? NaN : Inf) : x / baseline_x
    end
    return ratios
end

"""
    sensitivity_table(sweep::PPNSweep, isos; cycle = :final) -> DataFrame

Long-format sensitivity table: one row per `(reaction, isotope, factor)` for
every reaction in `sweep` and every isotope in `isos` (an `Isotope`,
isotope-name string, or a collection of either). Columns: `:reaction,
:isotope, :factor, :ratio`.
"""
function sensitivity_table(sweep::PPNSweep, isos; cycle = :final)
    isotope_list = isos isa Union{Isotope,AbstractString} ? [isos] : collect(isos)
    rows = NamedTuple[]
    for reaction in reactions(sweep), iso in isotope_list
        isotope = iso isa Isotope ? iso : Isotope(iso)
        for (factor, ratio) in sensitivity(sweep, reaction, isotope; cycle)
            push!(rows, (reaction = reaction, isotope = isotope_name(isotope), factor = factor, ratio = ratio))
        end
    end
    return sort(DataFrame(rows), [:reaction, :isotope, :factor])
end

"""
    iliadis_table(sweep::PPNSweep, isos; cycle = :final) -> DataFrame

Wide-format "Iliadis (2002) Table 8"-style sensitivity table: one row per
`(reaction, isotope)`, one column per tested factor.
"""
iliadis_table(sweep::PPNSweep, isos; cycle = :final) =
    unstack(sensitivity_table(sweep, isos; cycle), [:reaction, :isotope], :factor, :ratio)

"""
    rank_reactions(table::DataFrame, iso; metric = :max_abs_log_ratio) -> DataFrame

Rank the reactions in a long-format [`sensitivity_table`](@ref) `table` by how
much they move `iso`'s abundance, most sensitive first. The only implemented
`metric` is `:max_abs_log_ratio`: `max(|log10(ratio)|)` across every tested
factor for that reaction — the biggest swing away from baseline, regardless of
which factor or direction caused it.
"""
function rank_reactions(table::DataFrame, iso; metric::Symbol = :max_abs_log_ratio)
    metric === :max_abs_log_ratio ||
        throw(ArgumentError("unknown metric :$metric (only :max_abs_log_ratio is implemented)"))
    isotope_label = iso isa Isotope ? isotope_name(iso) : string(iso)
    rows = filter(:isotope => ==(isotope_label), table)
    isempty(rows) && throw(ArgumentError("no rows for isotope $isotope_label in table"))
    scored = combine(groupby(rows, :reaction), :ratio => (r -> maximum(abs.(log10.(r)))) => :score)
    return sort(scored, :score; rev = true)
end

"""
    sensitivity_plot(sweep::PPNSweep, reaction, isos; cycle = :final,
                      title = reaction, figure_size = (700, 500)) -> CM.Figure

Log-log plot of sensitivity ratio vs. rate factor for `reaction`, one line per
isotope in `isos`, with a horizontal reference line at ratio = 1.
"""
function sensitivity_plot(sweep::PPNSweep, reaction::AbstractString, isos; cycle = :final,
                           title = reaction, figure_size = (700, 500))
    isotope_list = isos isa Union{Isotope,AbstractString} ? [isos] : collect(isos)
    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1]; xlabel = "rate factor", ylabel = "X(factor) / X(baseline)",
                      title, xscale = log10, yscale = log10)
        CM.hlines!(ax, [1.0]; color = :gray, linestyle = :dash)
        for (i, iso) in enumerate(isotope_list)
            isotope = iso isa Isotope ? iso : Isotope(iso)
            sens = sensitivity(sweep, reaction, isotope; cycle)
            fs = sort(collect(keys(sens)))
            ys = [sens[f] for f in fs]
            mask = isfinite.(ys) .& (ys .> 0)
            any(mask) || continue
            color = NUGRID_PALETTE[mod1(i, length(NUGRID_PALETTE))]
            CM.lines!(ax, fs[mask], ys[mask]; label = isotope_name(isotope), color)
            CM.scatter!(ax, fs[mask], ys[mask]; color)
        end
        length(isotope_list) > 1 && CM.axislegend(ax; position = :lt, framevisible = false)
        fig
    end
end
