# sensitivity.jl: rank the varied rates by how much they control a quantity, and write the reports.

"""
    rank_reactions(ens, y; candidates=40) -> Vector{NamedTuple}

Rank the varied rates of `ens` by their influence on the per-sample quantity `y` (e.g.
`log10.(ratio(ens, "C12", "C13"))`). Samples where `y` is not finite are dropped. The Spearman
correlation of `y` with every rate's x is computed, and the mutual information ([`ksg_mi`](@ref), nats)
for the `candidates` rates with the largest |Spearman|. Returned sorted by MI, one
`(index, row, label, reaction, source, mi, spearman)` per candidate; `index` is the rate's position in
`ens.rates`/`ens.xr`. Ward et al. (2025) count a reaction as important for MI ≥ 0.3.
"""
function rank_reactions(ens::MCEnsemble, y::AbstractVector; candidates::Integer=40)
    length(y) == nsamples(ens) || throw(DimensionMismatch("y has $(length(y)) values, ensemble $(nsamples(ens))"))
    ok = isfinite.(y)
    yo, xo = y[ok], ens.xr[:, ok]
    rho = spearman_all(xo, yo)
    cand = sortperm(abs.(rho); rev=true)[1:min(candidates, length(rho))]
    mi = zeros(length(cand))
    Threads.@threads for c in eachindex(cand)
        mi[c] = ksg_mi(xo[cand[c], :], yo)
    end
    order = sortperm(mi; rev=true)
    [(index=cand[c], row=ens.rates[cand[c]].row, label=ens.rates[cand[c]].label,
      reaction=ens.rates[cand[c]].reaction, source=uncertainty_source(ens.rates[cand[c]].how),
      mi=mi[c], spearman=rho[cand[c]]) for c in order]
end

"`log10` of a ratio or abundance, `NaN` where it is not positive"
_log10pos(v) = [x > 0 && isfinite(x) ? log10(x) : NaN for x in v]

"""
    sensitivity_report(ens; pairs=WARD_PAIRS, isotopes=String[], X0=nothing, xi=0, reference=nothing,
                       reference_top=nothing, title="", top=10, candidates=40, outdir=ens.dir, io=stdout)

The Monte Carlo summary of a set of isotopic ratios (`pairs`, e.g. `[("C12", "C13")]`) and single
isotopes (`isotopes`, e.g. `["Na22"]`), written to `outdir`:
- `mc_ratios.txt`: per ratio the median, p16, p84 and factor1 = √(p84/p16) of the ratio; with `X0` and
  `xi` also the mixed ratio (Ward Eq. 4) and its factor2; per isotope the same for X. `reference`
  (`Dict(pair => (Xa, Xb, burn, factor1, mix, factor2))`, e.g. `WARD2025[:co]`) adds published values;
- `mc_sensitivity.txt` / `.csv`: per target, the reactions ranked by MI with log10(target) (see
  [`rank_reactions`](@ref)); `reference_top` (`Dict(pair => text)`) prints the published top reaction.

Returns `Dict(target => ranking)`.
"""
function sensitivity_report(ens::MCEnsemble; pairs=WARD_PAIRS, isotopes=String[], X0=nothing, xi::Real=0,
                            reference=nothing, reference_top=nothing, title::AbstractString="",
                            top::Integer=10, candidates::Integer=40, outdir::AbstractString=ens.dir, io::IO=stdout)
    n = nsamples(ens)
    mix = X0 !== nothing && xi > 0
    lines = [@sprintf("# %d Monte Carlo samples%s%s", n, isempty(title) ? "" : "; $title",
                      mix ? @sprintf("; mixed ratio with xi = %g", xi) : ""),
             @sprintf("# %-10s %10s %10s %10s %8s %9s %6s %10s %8s %9s %6s", "ratio", "median", "p16", "p84", "factor1",
                      "Ward", "W f1", "mix med", "factor2", "Ward mix", "W f2")]
    sens = String[]
    csv = ["target,rank,row,label,reaction,source,MI_nats,spearman"]
    result = Dict{String,Any}()
    function rank!(name, y, note)
        ranking = rank_reactions(ens, y; candidates)
        result[name] = ranking
        push!(sens, "\n$name   ($note)")
        push!(sens, @sprintf("  %6s %8s  row label reaction  [source of uncertainty]", "MI", "Spearman"))
        for r in ranking[1:min(top, end)]
            push!(sens, @sprintf("  %6.2f %+8.2f  %d %s %s  [%s]", r.mi, r.spearman, r.row, r.label, r.reaction, r.source))
        end
        for (k, r) in enumerate(ranking)
            push!(csv, @sprintf("%s,%d,%d,%s,\"%s\",%s,%.4f,%.4f", name, k, r.row, r.label, r.reaction,
                                ens.rates[r.index].how, r.mi, r.spearman))
        end
    end
    nan6 = (NaN, NaN, NaN, NaN, NaN, NaN)
    for (A, B) in pairs
        r = ratio(ens, A, B)
        s = spread(r)
        m = mix ? spread(ratio(ens, A, B; X0, xi)) : (median=NaN, factor=NaN)
        w = reference === nothing ? nan6 : get(reference, (A, B), nan6)
        push!(lines, @sprintf("  %-10s %10.3g %10.3g %10.3g %8.3f %9.3g %6.2f %10.3g %8.3f %9.3g %6.2f", "$A/$B",
                              s.median, s.p16, s.p84, s.factor, w[3], w[4], m.median, m.factor, w[5], w[6]))
        note = reference_top === nothing ? "log10 ratio" : "Ward: $(get(reference_top, (A, B), "-"))"
        rank!("$A/$B", _log10pos(r), note)
    end
    if !isempty(isotopes)
        push!(lines, @sprintf("# %-10s %10s %10s %10s %8s", "isotope", "median X", "p16", "p84", "factor"))
        for iso in isotopes
            x = abundance(ens, iso)
            s = spread(x)
            push!(lines, @sprintf("  %-10s %10.3g %10.3g %10.3g %8.3f", iso, s.median, s.p16, s.p84, s.factor))
            rank!(iso, _log10pos(x), "log10 X")
        end
    end
    hdr = [@sprintf("# %d samples. MI: KSG k=3 on log10(target) vs x, for the %d reactions with the largest", n,
                    candidates),
           "# |Spearman| per target. With few samples MI is biased and noisy: below ~1000 samples use this only as a check."]
    mkpath(outdir)
    write(joinpath(outdir, "mc_ratios.txt"), join(lines, '\n') * "\n")
    write(joinpath(outdir, "mc_sensitivity.txt"), join(vcat(hdr, sens), '\n') * "\n")
    write(joinpath(outdir, "mc_sensitivity.csv"), join(csv, '\n') * "\n")
    println(io, join(lines, '\n'))
    println(io, join(vcat(hdr, sens), '\n'))
    result
end
