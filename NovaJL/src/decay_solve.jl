# Per-isotope decay time inversion: for each isotope in a reference table (e.g. Iliadis
# Table 4), find the decay time t_i at which the PPN baseline abundance matches the
# reference, by interpolating across the decay_time_scan grid.
#
# The spread of {t_i} across isotopes measures how consistently the reference used a
# single decay time. Outliers flag isotopes where disagreement reflects rate differences
# rather than timing.

function _log_interpolate_crossing(t1, x1, t2, x2, x_target)
    # Log-log interpolation between (t1,x1) and (t2,x2) to find t where x == x_target.
    # Falls back to linear when t1 == 0 (can't take log).
    if t1 == 0.0
        abs(x2 - x1) < 1e-300 && return 0.5 * t2
        frac = (x_target - x1) / (x2 - x1)
        return frac * t2
    end
    log_t1, log_t2 = log10(t1), log10(t2)
    log_x1 = x1 > 0 ? log10(x1) : -300.0
    log_x2 = x2 > 0 ? log10(x2) : -300.0
    log_xt = x_target > 0 ? log10(x_target) : -300.0
    abs(log_x2 - log_x1) < 1e-10 && return 10^((log_t1 + log_t2) / 2)
    frac = (log_xt - log_x1) / (log_x2 - log_x1)
    return 10^(log_t1 + frac * (log_t2 - log_t1))
end

function _find_crossings(times, curve, x_target)
    # Return all brackets (t_solved, t_lower, t_upper) where the curve crosses x_target.
    crossings = Tuple{Float64,Float64,Float64}[]
    for i in 1:length(times)-1
        x1, x2 = curve[i], curve[i+1]
        (x1 - x_target) * (x2 - x_target) <= 0 || continue
        t_s = _log_interpolate_crossing(times[i], x1, times[i+1], x2, x_target)
        push!(crossings, (t_s, times[i], times[i+1]))
    end
    return crossings
end

function _curve_direction(curve)
    length(curve) < 2 && return "flat"
    diffs = diff(curve)
    any_pos = any(d -> d > 1e-300, diffs)
    any_neg = any(d -> d < -1e-300, diffs)
    any_pos && any_neg && return "non_monotonic"
    any_pos && return "increasing"
    any_neg && return "decreasing"
    return "flat"
end

"""
    solve_isotope_decay_times(paths; model, min_reference_abundance) -> NamedTuple

For each isotope in the Iliadis `model` final abundance table, interpolate the
decay_time_scan grid to find t_i where X_PPN(t_i) = X_Iliadis.

Returned NamedTuple fields:
  - `per_isotope`        — DataFrame, one row per Iliadis isotope
  - `consensus_time`     — geometric mean of t_i over isotopes with a valid crossing (seconds)
  - `std_log10_t`        — std dev of log10(t_i) across those isotopes
  - `n_with_crossing`    — number of isotopes where a crossing was found
  - `n_without_crossing` — number of isotopes where no crossing exists in the scan range

Per-isotope DataFrame columns:
  isotope, X_iliadis, X_ppn_t0, X_ppn_tmax,
  t_solved, log10_t_solved,
  t_lower_bracket, t_upper_bracket, log10_bracket_width,
  direction (increasing/decreasing/flat/non_monotonic),
  has_crossing, n_crossings
"""
function solve_isotope_decay_times(
    paths;
    model = "JCH1",
    min_reference_abundance = 1.0e-30,
)
    isfile(paths.decay_time_scan_manifest) || error(
        "missing $(paths.decay_time_scan_manifest); run " *
        "`julia tools/decay_time_scan.jl --nova <case>` first",
    )

    manifest = CSV.read(paths.decay_time_scan_manifest, DataFrame)
    df_iliadis = read_iliadis_final_abundance(paths, model)

    # Load every successful scan time point
    time_frames = Dict{Float64,DataFrame}()
    for row in eachrow(manifest)
        row.status != "ok" && continue
        out = String(row.output)
        isfile(out) || continue
        time_frames[row.decay_time_seconds] = read_iso_massf(out)
    end
    isempty(time_frames) && error("no successful decay time scan results in $(paths.decay_time_scan_manifest)")

    times_sorted = sort(collect(keys(time_frames)))

    rows = NamedTuple[]
    for ili_row in eachrow(df_iliadis)
        iso   = ili_row.isotope
        X_ili = ili_row.X
        X_ili >= min_reference_abundance || continue

        # Build (time, abundance) curve for this isotope
        valid_t = Float64[]
        curve   = Float64[]
        for t in times_sorted
            df_t = time_frames[t]
            idx  = findfirst(==(iso), df_t.isotope)
            idx === nothing && continue
            X = df_t.X[idx]
            X > 0 || continue
            push!(valid_t, t)
            push!(curve, X)
        end
        length(curve) < 2 && continue

        direction = _curve_direction(curve)
        crossings = _find_crossings(valid_t, curve, X_ili)

        if isempty(crossings)
            push!(rows, (
                isotope              = iso,
                X_iliadis            = X_ili,
                X_ppn_t0             = curve[1],
                X_ppn_tmax           = curve[end],
                t_solved             = NaN,
                log10_t_solved       = NaN,
                t_lower_bracket      = NaN,
                t_upper_bracket      = NaN,
                log10_bracket_width  = NaN,
                direction            = direction,
                has_crossing         = false,
                n_crossings          = 0,
            ))
        else
            # Pick the first (earliest) crossing; flag non-monotonic cases via n_crossings > 1
            t_s, t_lo, t_hi = crossings[1]

            # log10_t_solved: use a floor of 1 s so log10 is finite for near-zero solves
            log10_t_s = log10(max(t_s, 1.0))

            # Bracket width in log10 space (precision of the solve given the grid density)
            log10_lo = log10(max(t_lo, 1.0))
            log10_hi = log10(max(t_hi, 1.0))
            bracket_width = abs(log10_hi - log10_lo)

            push!(rows, (
                isotope              = iso,
                X_iliadis            = X_ili,
                X_ppn_t0             = curve[1],
                X_ppn_tmax           = curve[end],
                t_solved             = t_s,
                log10_t_solved       = log10_t_s,
                t_lower_bracket      = t_lo,
                t_upper_bracket      = t_hi,
                log10_bracket_width  = bracket_width,
                direction            = direction,
                has_crossing         = true,
                n_crossings          = length(crossings),
            ))
        end
    end

    isempty(rows) && error("no isotopes found with sufficient abundance in both scan and Iliadis data")
    df = sort(DataFrame(rows), [:isotope])

    valid_rows  = filter(r -> r.has_crossing, df)
    valid_log10 = filter(isfinite, valid_rows[!, :log10_t_solved])
    n_valid = length(valid_log10)

    mean_log10 = n_valid > 0 ? sum(valid_log10) / n_valid : NaN
    std_log10  = n_valid >= 2 ?
        sqrt(sum((valid_log10 .- mean_log10) .^ 2) / (n_valid - 1)) : NaN

    return (
        per_isotope          = df,
        consensus_time       = isfinite(mean_log10) ? 10^mean_log10 : NaN,
        std_log10_t          = std_log10,
        n_with_crossing      = n_valid,
        n_without_crossing   = nrow(df) - n_valid,
    )
end

"""
    plot_decay_time_solve(result; title) -> Plot

Scatter plot of t_solved per isotope with the consensus time and ±1σ band overlaid.
Isotopes without a crossing are excluded from the plot but noted in the title.
"""
function plot_decay_time_solve(result; title = "Per-isotope solved decay time vs Iliadis")
    df_valid = filter(r -> r.has_crossing && isfinite(r.log10_t_solved), result.per_isotope)
    isempty(df_valid) && error("no isotopes with valid crossing to plot")

    t_consensus = result.consensus_time
    std_log10   = result.std_log10_t
    n_no        = result.n_without_crossing

    subtitle = n_no > 0 ? " ($n_no isotopes without crossing not shown)" : ""

    p = scatter(
        df_valid.isotope,
        df_valid.t_solved;
        yscale     = :log10,
        xlabel     = "Isotope",
        ylabel     = "Solved decay time (s)",
        title      = title * subtitle,
        label      = "t_i (solved)",
        markersize = 5,
        xrotation  = 45,
    )

    hline!(
        p, [t_consensus];
        label     = "Consensus: $(round(t_consensus / 3600; digits = 1)) h",
        linestyle = :dash,
        color     = :red,
        linewidth = 2,
    )

    if isfinite(std_log10)
        t_lo = 10^(log10(t_consensus) - std_log10)
        t_hi = 10^(log10(t_consensus) + std_log10)
        hline!(p, [t_lo, t_hi]; label = "±1σ (log10)", linestyle = :dot, color = :orange, linewidth = 1.5)
    end

    return p
end

"""
    decay_solve_summary(result) -> DataFrame

Compact summary table: consensus time, std, and counts.
"""
function decay_solve_summary(result)
    t_c = result.consensus_time
    return DataFrame(
        consensus_time_s     = [t_c],
        consensus_time_h     = [t_c / 3600],
        consensus_time_label = [isfinite(t_c) ? _format_decay_time(t_c) : "N/A"],
        std_log10_t          = [result.std_log10_t],
        n_with_crossing      = [result.n_with_crossing],
        n_without_crossing   = [result.n_without_crossing],
    )
end

function _format_decay_time(seconds)
    seconds == 0.0 && return "0 s"
    minute, hour, day, year = 60.0, 3600.0, 86400.0, 31557600.0
    if seconds < hour
        return "$(round(seconds / minute; digits=2)) min"
    elseif seconds < day
        return "$(round(seconds / hour; digits=2)) h"
    elseif seconds < year
        return "$(round(seconds / day; digits=1)) d"
    end
    return "$(round(seconds / year; digits=2)) yr"
end

# ---------------------------------------------------------------------------
# Analytic solver
# ---------------------------------------------------------------------------

# Half-lives (seconds) for unstable isotopes in the nova network (Z ≤ 20).
# Values from NUBASE2020 (Kondev et al. 2021, Chinese Phys. C 45, 030001).
# Isotopes not listed here are assumed stable (no analytic solution possible).
const NOVA_HALF_LIVES_S = Dict{String,Float64}(
    "BE-7"  => 53.22   * 86400.0,            # 53.22 d
    "N-13"  => 9.965   * 60.0,               # 9.965 min
    "O-14"  => 70.62,                         # 70.62 s
    "O-15"  => 122.24,                        # 122.24 s
    "F-17"  => 64.49,                         # 64.49 s
    "NE-18" => 1.672,                         # 1.672 s
    "NE-19" => 17.26,                         # 17.26 s
    "F-18"  => 109.739 * 60.0,               # 109.739 min
    "NA-20" => 0.4479,                        # 0.4479 s
    "NA-21" => 22.49,                         # 22.49 s
    "NA-22" => 2.6018  * 365.25 * 86400.0,   # 2.6018 yr
    "NA-24" => 14.997  * 3600.0,             # 14.997 h
    "MG-23" => 11.317,                        # 11.317 s
    "MG-27" => 567.5,                         # 567.5 s
    "AL-25" => 7.183,                         # 7.183 s
    "AL-26" => 717300.0 * 365.25 * 86400.0,  # 0.7173 Myr (ground state)
    "AL-28" => 134.4,                         # 134.4 s
    "SI-27" => 4.16,                          # 4.16 s
    "SI-31" => 9440.0,                        # 2.622 h
    "P-30"  => 149.88,                        # 149.88 s
    "P-32"  => 14.268  * 86400.0,            # 14.268 d
    "P-33"  => 24.97   * 86400.0,            # 24.97 d
    "S-31"  => 2.5534,                        # 2.5534 s
    "S-35"  => 87.24   * 86400.0,            # 87.24 d
    "CL-33" => 150.7,                         # 2.511 min = 150.7 s
    "CL-34" => 1.5264,                        # 1.5264 s
    "AR-37" => 35.011  * 86400.0,            # 35.011 d
    "K-37"  => 1.2267,                        # 1.2267 s
    "K-38"  => 7.636   * 60.0,              # 7.636 min
    "CA-39" => 0.8605,                        # 0.8605 s
    "CA-41" => 1.02e5  * 365.25 * 86400.0,  # ~102,000 yr
)

"""
    solve_isotope_decay_times_analytic(paths; model, min_reference_abundance) -> NamedTuple

Analytic per-isotope decay time solver.  For each unstable isotope in the Iliadis
reference, applies the radioactive decay law directly:

    t_i = T½_i / ln(2) × ln(X_ppn₀_i / X_iliadis_i)

where X_ppn₀ is the PPN abundance at t=0 (end of nova, before any decay) read from
the t=0 entry in the decay_time_scan.

Returned NamedTuple fields:
  - `per_isotope`    — DataFrame, one row per Iliadis isotope (including stable)
  - `consensus_time` — geometric mean of t_i over isotopes with flag "ok"
  - `std_log10_t`    — std dev of log10(t_i) across those isotopes
  - `n_ok`           — isotopes with a valid solved t
  - `n_stable`       — isotopes with no solution because λ = 0
  - `n_underproduction` — isotopes where X_ppn₀ < X_iliadis (PPN underproduces)
  - `n_unknown`      — isotopes not in the half-life table and not obviously stable

Per-isotope DataFrame columns:
  isotope, X_iliadis, X_ppn_0, T_half_s, lambda_per_s,
  t_solved, log10_t_solved, flag
  (flag: "ok" | "stable" | "underproduction" | "zero_ppn" | "unknown_halflife")
"""
function solve_isotope_decay_times_analytic(
    paths;
    model = "JCH1",
    min_reference_abundance = 1.0e-30,
)
    isfile(paths.decay_time_scan_manifest) || error(
        "missing $(paths.decay_time_scan_manifest); run " *
        "`julia tools/decay_time_scan.jl --nova <case>` first",
    )

    manifest = CSV.read(paths.decay_time_scan_manifest, DataFrame)
    df_iliadis = read_iliadis_final_abundance(paths, model)

    # X_ppn_0 comes from the t=0 scan point (end of nova, no decay applied yet)
    t0_rows = filter(r -> r.decay_time_seconds == 0.0 && r.status == "ok", manifest)
    isempty(t0_rows) && error(
        "no t=0 entry in scan manifest; re-run `julia tools/decay_time_scan.jl` " *
        "with the default time grid (which includes t=0)",
    )
    df_t0 = read_iso_massf(String(t0_rows[1, :output]))

    rows = NamedTuple[]
    for ili_row in eachrow(df_iliadis)
        iso   = ili_row.isotope
        X_ili = ili_row.X
        X_ili >= min_reference_abundance || continue

        idx    = findfirst(==(iso), df_t0.isotope)
        X_ppn0 = idx === nothing ? 0.0 : df_t0.X[idx]

        T_half = get(NOVA_HALF_LIVES_S, iso, nothing)
        λ = if T_half === nothing
            nothing
        elseif T_half == Inf
            0.0
        else
            log(2.0) / T_half
        end

        t_solved, log10_t, flag = if T_half === nothing
            NaN, NaN, "unknown_halflife"
        elseif λ == 0.0
            NaN, NaN, "stable"
        elseif X_ppn0 <= 0.0
            NaN, NaN, "zero_ppn"
        elseif X_ppn0 < X_ili
            NaN, NaN, "underproduction"
        else
            t = log(X_ppn0 / X_ili) / λ
            t, log10(max(t, 1.0)), "ok"
        end

        push!(rows, (
            isotope      = iso,
            X_iliadis    = X_ili,
            X_ppn_0      = X_ppn0,
            T_half_s     = T_half === nothing ? NaN : T_half,
            lambda_per_s = λ === nothing ? NaN : λ,
            t_solved     = t_solved,
            log10_t_solved = log10_t,
            flag         = flag,
        ))
    end

    isempty(rows) && error("no isotopes found with sufficient abundance in Iliadis data")
    df = sort(DataFrame(rows), :isotope)

    ok_rows  = filter(r -> r.flag == "ok", df)
    log10_ok = filter(isfinite, ok_rows[!, :log10_t_solved])
    n_ok     = length(log10_ok)

    mean_log10 = n_ok > 0 ? sum(log10_ok) / n_ok : NaN
    std_log10  = n_ok >= 2 ?
        sqrt(sum((log10_ok .- mean_log10) .^ 2) / (n_ok - 1)) : NaN

    return (
        per_isotope       = df,
        consensus_time    = isfinite(mean_log10) ? 10^mean_log10 : NaN,
        std_log10_t       = std_log10,
        n_ok              = length(log10_ok),
        n_stable          = sum(r -> r.flag == "stable",          eachrow(df)),
        n_underproduction = sum(r -> r.flag == "underproduction", eachrow(df)),
        n_unknown         = sum(r -> r.flag == "unknown_halflife", eachrow(df)),
    )
end

"""
    compare_decay_solve_methods(paths; model, min_reference_abundance) -> NamedTuple

Run both the interpolation and analytic solvers and join their results on isotope.
Only isotopes where at least one method produces a finite t_i appear in the comparison.

Returned NamedTuple fields:
  - `merged`            — joined DataFrame (see columns below)
  - `interp_consensus`  — consensus time from interpolation solver
  - `analytic_consensus`— consensus time from analytic solver
  - `interp_std`        — std_log10_t from interpolation solver
  - `analytic_std`      — std_log10_t from analytic solver
  - `interp`            — full result from solve_isotope_decay_times
  - `analytic`          — full result from solve_isotope_decay_times_analytic

Merged DataFrame key columns:
  isotope, T_half_s, flag (analytic status),
  log10_t_interp, log10_bracket_width, direction,
  log10_t_analytic,
  log10_diff (|analytic − interp| in log10 space),
  agree_within_bracket (Bool: analytic within ±bracket of interp)
"""
function compare_decay_solve_methods(
    paths;
    model = "JCH1",
    min_reference_abundance = 1.0e-30,
)
    interp  = solve_isotope_decay_times(paths;          model = model, min_reference_abundance = min_reference_abundance)
    analytic = solve_isotope_decay_times_analytic(paths; model = model, min_reference_abundance = min_reference_abundance)

    df_i = select(
        interp.per_isotope,
        :isotope, :X_iliadis, :X_ppn_t0 => :X_ppn_0_interp,
        :log10_t_solved => :log10_t_interp,
        :log10_bracket_width,
        :has_crossing,
        :direction,
    )
    df_a = select(
        analytic.per_isotope,
        :isotope,
        :T_half_s,
        :log10_t_solved => :log10_t_analytic,
        :flag,
    )

    merged = innerjoin(df_i, df_a, on = :isotope)

    merged.log10_diff = [
        isfinite(r.log10_t_interp) && isfinite(r.log10_t_analytic) ?
            abs(r.log10_t_analytic - r.log10_t_interp) : NaN
        for r in eachrow(merged)
    ]

    merged.agree_within_bracket = [
        isfinite(r.log10_diff) && isfinite(r.log10_bracket_width) &&
            r.log10_diff <= r.log10_bracket_width
        for r in eachrow(merged)
    ]

    return (
        comparison         = sort(merged, :isotope),
        interp_consensus   = interp.consensus_time,
        analytic_consensus = analytic.consensus_time,
        interp_std         = interp.std_log10_t,
        analytic_std       = analytic.std_log10_t,
        interp_solve       = interp,
        analytic_solve     = analytic,
    )
end

"""
    plot_method_comparison(cmp; title) -> Plot

Side-by-side log10(t_solved) for the interpolation vs analytic method.
Only isotopes where both methods give a finite t_i are shown.
Agreement within the interpolation bracket is highlighted in green.
"""
function plot_method_comparison(cmp; title = "Decay time solve: interpolation vs analytic")
    df = filter(r -> isfinite(r.log10_t_interp) && isfinite(r.log10_t_analytic), cmp.comparison)
    isempty(df) && error("no isotopes with finite t_i from both methods")

    colors = [r.agree_within_bracket ? :seagreen : :red for r in eachrow(df)]
    lo = min(minimum(df.log10_t_interp), minimum(df.log10_t_analytic)) - 0.5
    hi = max(maximum(df.log10_t_interp), maximum(df.log10_t_analytic)) + 0.5

    p = scatter(
        df.log10_t_analytic,
        df.log10_t_interp;
        xlabel     = "log₁₀ t_analytic (s)",
        ylabel     = "log₁₀ t_interp (s)",
        title      = title,
        label      = "",
        markercolor = colors,
        markersize = 8,
        xlims      = (lo, hi),
        ylims      = (lo, hi),
        aspect_ratio = :equal,
    )

    plot!(p, [lo, hi], [lo, hi]; linestyle = :dash, color = :black, label = "exact agreement")

    for r in eachrow(df)
        isfinite(r.log10_bracket_width) || continue
        plot!(
            p,
            [r.log10_t_analytic, r.log10_t_analytic],
            [r.log10_t_interp - r.log10_bracket_width / 2,
             r.log10_t_interp + r.log10_bracket_width / 2];
            color = :gray, linewidth = 1, label = "",
        )
    end

    # Label each point with the isotope name
    for r in eachrow(df)
        annotate!(p, r.log10_t_analytic + 0.05, r.log10_t_interp, text(r.isotope, 7, :left))
    end

    return p
end
