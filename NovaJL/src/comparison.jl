using CSV
using JSON

const MODEL_TABLE8_FILES = Dict(
    "P1" => "table_05_sensitivity_P1.dat",
    "P2" => "table_06_sensitivity_P2.dat",
    "S1" => "table_07_sensitivity_S1.dat",
    "JCH1" => "table_08_sensitivity_JCH1.dat",
    "JCH2" => "table_09_sensitivity_JCH2.dat",
    "JH1" => "table_10_sensitivity_JH1.dat",
    "JH2" => "table_11_sensitivity_JH2.dat",
)

const TABLE8_FACTOR_COLUMNS = [:factor_100, :factor_10, :factor_2, Symbol("factor_0.5"), Symbol("factor_0.1"), Symbol("factor_0.01")]

function nova_case_paths(case_analysis_dir = pwd(); decay_scan_dir = "decay_time_scan")
    analysis_dir = abspath(case_analysis_dir)
    case_dir = dirname(analysis_dir)
    nova_cases_dir = dirname(case_dir)
    sensitivity_root = dirname(nova_cases_dir)
    data_root = joinpath(sensitivity_root, "data", "iliadis2002")
    return (
        analysis_dir = analysis_dir,
        case_dir = case_dir,
        runs_dir = joinpath(case_dir, "runs"),
        decay_time_scan_manifest = joinpath(case_dir, decay_scan_dir, "decay_time_scan_manifest.csv"),
        decay_time_scan_csv = joinpath(analysis_dir, "results", "decay_time_scan_jl", "decay_time_scan.csv"),
        config_path = joinpath(case_dir, "config", "reaction_plan.json"),
        data_root = data_root,
        results_dir = joinpath(analysis_dir, "results", "iliadis_comparison_jl"),
    )
end

function latest_validated_decay_run(paths)
    isfile(paths.decay_time_scan_csv) || error(
        "missing $(paths.decay_time_scan_csv); run `julia tools/decay_time_scan.jl --nova <case>` " *
        "then the 00 decay-time-calibration notebook first",
    )
    df = CSV.read(paths.decay_time_scan_csv, DataFrame)
    isempty(df) && error("no rows in $(paths.decay_time_scan_csv)")
    "output" in names(df) || error("no `output` column in $(paths.decay_time_scan_csv) — re-run the 00 notebook to regenerate it")
    return String(df.output[1])
end

function isotope_z(isotope)
    isotope == "n" && return 0
    parts = split(isotope, "-")
    length(parts) == 2 || error("cannot parse element symbol out of isotope label $(repr(isotope))")
    symbol = uppercase(parts[1])
    idx = findfirst(==(symbol), uppercase.(ELEMENT_SYMBOLS))
    idx === nothing && error("unknown element symbol $(repr(parts[1])) in isotope label $(repr(isotope))")
    return idx - 1
end

function read_iliadis_final_abundance(paths, model)
    path = joinpath(paths.data_root, "ppn_final_abundances", "iso_massf_iliadis2002_final_$(model).DAT")
    isfile(path) || error("missing Iliadis final abundance table: $path")
    return read_iso_massf(path)
end

function compare_baseline_to_iliadis(df_ppn::DataFrame, df_iliadis::DataFrame; max_z = 20, min_reference_abundance = 1.0e-30)
    ppn = select(df_ppn, :isotope, :X => :ppn)
    iliadis = select(df_iliadis, :isotope, :X => :iliadis)
    merged = innerjoin(ppn, iliadis, on = :isotope)
    merged.z = isotope_z.(merged.isotope)
    merged = filter(row -> row.z <= max_z && row.iliadis >= min_reference_abundance && row.ppn > 0, merged)

    merged.ratio = merged.ppn ./ merged.iliadis
    merged.log10ratio = log10.(merged.ratio)
    merged.absdev = abs.(merged.log10ratio)
    return sort(merged, [:z, :isotope])
end

function score_baseline(df_compare::DataFrame)
    return DataFrame([score_group("baseline", df_compare)])
end

function plot_baseline_comparison(df_compare::DataFrame;
        within_10_percent = 0.1,
        within_factor_2   = 2.0,
        title             = "Baseline final abundance",
        xlabel            = "Iliadis (2002) final abundance",
        ylabel            = "PPN decayed final abundance",
    )
    lo = minimum(vcat(df_compare.ppn, df_compare.iliadis)) / 1.5
    hi = maximum(vcat(df_compare.ppn, df_compare.iliadis)) * 1.5
    xs = [lo, hi]

    p = plot(
        xscale = :log10,
        yscale = :log10,
        xlims = (lo, hi),
        ylims = (lo, hi),
        xlabel = xlabel,
        ylabel = ylabel,
        title = title,
        legend = false,
        aspect_ratio = :equal,
    )
    plot!(p, xs, xs .* within_factor_2; fillrange = xs ./ within_factor_2, fillalpha = 0.18, linealpha = 0, color = :gray, label = "")
    plot!(p, xs, xs .* (1 + within_10_percent); fillrange = xs .* (1 - within_10_percent), fillalpha = 0.25, linealpha = 0, color = :seagreen, label = "")
    plot!(p, xs, xs; linestyle = :dash, color = :black, label = "")
    scatter!(p, df_compare.iliadis, df_compare.ppn; label = "", markersize = 4)
    return p
end

function decay_time_scan_results(paths; model = "JCH1", max_z = 20, min_reference_abundance = 1.0e-30, write_results = true)
    isfile(paths.decay_time_scan_manifest) || error(
        "missing $(paths.decay_time_scan_manifest); run `julia tools/decay_time_scan.jl --nova <case>` first",
    )
    manifest = CSV.read(paths.decay_time_scan_manifest, DataFrame)
    df_iliadis = read_iliadis_final_abundance(paths, model)

    rows = NamedTuple[]
    comparisons = Dict{Float64,DataFrame}()
    for entry in eachrow(manifest)
        entry.status != "ok" && continue
        df_ppn = read_iso_massf(entry.output)
        cmp = compare_baseline_to_iliadis(df_ppn, df_iliadis; max_z = max_z, min_reference_abundance = min_reference_abundance)
        comparisons[entry.decay_time_seconds] = cmp
        score = score_group("baseline", cmp)
        push!(
            rows,
            (
                decay_time_seconds = entry.decay_time_seconds,
                decay_time_label = entry.decay_time_label,
                output = String(entry.output),
                common_isotopes = score.common_isotopes,
                finite_ratios = score.finite_ratios,
                rms_log10 = score.rms_log10,
                median_abs_log10 = score.median_abs_log10,
                within_10_percent = score.within_10_percent,
                within_10_percent_fraction = score.within_10_percent_fraction,
                within_factor_2 = score.within_factor_2,
                within_factor_2_fraction = score.within_factor_2_fraction,
            ),
        )
    end
    isempty(rows) && error("no usable decay times found in $(paths.decay_time_scan_manifest)")
    scan = sort(DataFrame(rows), :rms_log10)

    if write_results
        mkpath(dirname(paths.decay_time_scan_csv))
        CSV.write(paths.decay_time_scan_csv, scan)
    end

    best_time = scan.decay_time_seconds[1]
    return (scan = scan, comparisons = comparisons, best_comparison = comparisons[best_time], best_time = best_time)
end

function plot_decay_time_scan(scan::DataFrame)
    sorted = filter(r -> r.decay_time_seconds > 0, sort(scan, :decay_time_seconds))
    p = plot(
        sorted.decay_time_seconds,
        sorted.rms_log10;
        seriestype = :scatter,
        xscale = :log10,
        xlabel = "Decay time (s)",
        ylabel = "RMS log10(PPN / Iliadis)",
        title = "Decay time scan vs Iliadis Table 4",
        label = "",
        markersize = 6,
    )
    plot!(p, sorted.decay_time_seconds, sorted.rms_log10; label = "", color = :gray, linealpha = 0.5)

    best = sort(scan, :rms_log10)[1, :]
    scatter!(
        p,
        [best.decay_time_seconds],
        [best.rms_log10];
        label = "best: $(best.decay_time_label)",
        markersize = 10,
        markershape = :star5,
        color = :red,
    )
    return p
end

function iliadis_isotope_to_io_style(iso)
    m = match(r"^(\d+)([A-Za-z]+)$", iso)
    m === nothing && return uppercase(iso)
    return uppercase(m.captures[2]) * "-" * m.captures[1]
end

function load_reaction_plan(path)
    isfile(path) || error("missing reaction plan: $path")
    data = JSON.parsefile(path)
    reactions = data["reactions"]
    entries = reactions isa AbstractDict ? values(reactions) : reactions

    rows = NamedTuple[]
    for reaction in entries
        get(reaction, "match_status", "") == "missing_from_network" && continue
        name = get(reaction, "name", get(reaction, "reaction_id", nothing))
        name === nothing && continue
        article = get(reaction, "article_reaction", name)
        isotopes = [iliadis_isotope_to_io_style(string(iso)) for iso in get(reaction, "affected_isotopes", String[])]
        factors = Float64.(get(reaction, "factors", Float64[]))
        notes = get(reaction, "notes", String[])
        push!(
            rows,
            (
                name = String(name),
                article_reaction = String(article),
                network_name = String(get(reaction, "network_name", name)),
                isotopes = isotopes,
                factors = factors,
                index = get(reaction, "index", missing),
                product_was_remapped = get(reaction, "product_was_remapped", false),
                notes = join(string.(notes), "; "),
            ),
        )
    end
    return rows
end

const _ELEMENT_SYMBOLS_UPPER = Set(uppercase(s) for s in ELEMENT_SYMBOLS)

function decode_isomer_label(label)
    s = string(label)
    m = match(r"^(\d+)([A-Za-z]+)$", s)
    m === nothing && return (isotope = s, isomer = "ground_or_unspecified")
    mass, element_raw = m.captures

    if uppercase(element_raw) in _ELEMENT_SYMBOLS_UPPER
        return (isotope = mass * uppercasefirst(lowercase(element_raw)), isomer = "ground_or_unspecified")
    end

    if length(element_raw) >= 2 && element_raw[end] in ('g', 'm', 'G', 'M') &&
       uppercase(element_raw[1:end-1]) in _ELEMENT_SYMBOLS_UPPER
        base = element_raw[1:end-1]
        state = lowercase(element_raw[end]) == "g" ? "ground" : "metastable"
        return (isotope = mass * uppercasefirst(lowercase(base)), isomer = state)
    end

    return (isotope = s, isomer = "ground_or_unspecified")
end

function isomer_reaction_table(reaction_plan)
    rows = NamedTuple[]
    for r in reaction_plan
        tokens = split(r.name, "_")
        for token in tokens
            decoded = decode_isomer_label(token)
            decoded.isomer == "ground_or_unspecified" && continue
            push!(
                rows,
                (
                    name = r.name,
                    article_reaction = r.article_reaction,
                    network_name = r.network_name,
                    tagged_token = String(token),
                    isotope = decoded.isotope,
                    isomer = decoded.isomer,
                    product_was_remapped = r.product_was_remapped,
                    notes = r.notes,
                ),
            )
        end
    end
    isempty(rows) && return DataFrame(
        name = String[], article_reaction = String[], network_name = String[], tagged_token = String[],
        isotope = String[], isomer = String[], product_was_remapped = Bool[], notes = String[],
    )
    return DataFrame(rows)
end

function final_cycle_label(xtime_path)
    isfile(xtime_path) || error("missing x-time.dat: $xtime_path")
    final_cycle = nothing
    for line in eachline(xtime_path)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue
        final_cycle = parse(Int, split(stripped)[1])
    end
    final_cycle === nothing && error("no cycle rows found in $xtime_path")
    return lpad(string(final_cycle), 5, '0')
end

function build_ppn_table8(run_dir::AbstractString, reaction_plan; verbose = false, iso_massf_filename = nothing)
    if iso_massf_filename === nothing
        cycle = final_cycle_label(joinpath(run_dir, "baseline", "x-time.dat"))
        iso_massf_filename = "iso_massf$(cycle).DAT"
    end
    base_path = joinpath(run_dir, "baseline", iso_massf_filename)

    reaction_isotopes = [(r.name, r.isotopes) for r in reaction_plan]
    all_factors = sort(unique(vcat([r.factors for r in reaction_plan]...)))

    return output_sensitivity_table(
        reaction_isotopes;
        reaction_run_path = run_dir,
        factors = all_factors,
        base_path = base_path,
        iso_massf_filename = iso_massf_filename,
        verbose = verbose,
    )
end

function _parse_factor_cell(x)
    ismissing(x) && return missing
    s = strip(string(x))
    s == "..." && return missing
    return parse(Float64, s)
end

function read_iliadis_table8(paths, model)
    haskey(MODEL_TABLE8_FILES, model) ||
        error("unknown Iliadis model $model; expected one of $(join(sort(collect(keys(MODEL_TABLE8_FILES))), ", "))")
    path = joinpath(paths.data_root, "tables", MODEL_TABLE8_FILES[model])
    isfile(path) || error("missing Iliadis table: $path")

    df = CSV.read(path, DataFrame; delim = '\t', comment = "#", types = String)
    rename!(df, :reaction_id => :reaction, :reaction => :reaction_label, :isotope_label => :isotope_io)
    for col in TABLE8_FACTOR_COLUMNS
        df[!, col] = _parse_factor_cell.(df[!, col])
    end
    df.isotope_io = uppercase.(df.isotope_io)
    return df[:, [:reaction, :reaction_label, :isotope, :isotope_io, TABLE8_FACTOR_COLUMNS...]]
end

function read_iliadis_table3(paths)
    path = joinpath(paths.data_root, "tables", "table_03_reaction_rate_uncertainties.dat")
    isfile(path) || error("missing Iliadis table: $path")

    df = CSV.read(path, DataFrame; delim = '\t', comment = "#", types = String)
    rename!(df, :reaction_id => :reaction, :reaction => :reaction_label)
    df.factor_up = _parse_factor_cell.(df.factor_up)
    df.factor_down = _parse_factor_cell.(df.factor_down)
    return df[:, [:reaction, :reaction_label, :factor_up, :factor_down]]
end

function classify_goodness(ratio, factor, table3_up, table3_down)
    ismissing(ratio) && return "missing"
    within_10 = 0.9 <= ratio <= 1.1
    if ismissing(table3_up) || ismissing(table3_down)
        return within_10 ? "within_10_unknown_table3" : "outside_10_unknown_table3"
    end
    within_table3 = table3_down <= factor <= table3_up
    within_10 && within_table3 && return "good_within_table3"
    within_10 && return "good_but_factor_outside_table3"
    within_table3 && return "bad_within_table3"
    return "bad_and_factor_outside_table3"
end

_factor_value_from_ppn_column(name) = parse(Float64, string(name))
_factor_value_from_iliadis_column(name) = parse(Float64, replace(string(name), "factor_" => ""))

function compare_to_iliadis(df_ppn::DataFrame, df_iliadis::DataFrame; df_table3 = nothing)
    ppn_factor_cols = [c for c in names(df_ppn) if !(c in ("reaction", "isotope"))]
    long_ppn = stack(df_ppn, ppn_factor_cols; variable_name = :factor_col, value_name = :ppn)
    long_ppn.factor = _factor_value_from_ppn_column.(long_ppn.factor_col)
    long_ppn.isotope = uppercase.(string.(long_ppn.isotope))
    long_ppn = select(long_ppn, :reaction, :isotope, :factor, :ppn)

    ili_factor_cols = string.(TABLE8_FACTOR_COLUMNS)
    ili_base = select(df_iliadis, :reaction, :reaction_label, :isotope_io => :isotope, ili_factor_cols...)
    long_ili = stack(ili_base, ili_factor_cols; variable_name = :factor_col, value_name = :iliadis)
    long_ili.factor = _factor_value_from_iliadis_column.(long_ili.factor_col)
    long_ili = select(long_ili, :reaction, :reaction_label, :isotope, :factor, :iliadis)

    merged = innerjoin(long_ppn, long_ili, on = [:reaction, :isotope, :factor])
    merged = filter(row -> !ismissing(row.ppn) && !ismissing(row.iliadis), merged)

    merged.ratio = merged.ppn ./ merged.iliadis
    merged.log10ratio = log10.(merged.ratio)
    merged.absdev = abs.(merged.log10ratio)

    if df_table3 !== nothing
        table3_lookup = Dict(row.reaction => (row.factor_up, row.factor_down) for row in eachrow(df_table3))
        merged.classification = [
            classify_goodness(row.ratio, row.factor, get(table3_lookup, row.reaction, (missing, missing))...)
            for row in eachrow(merged)
        ]
    end

    return sort(merged, [:reaction, :isotope, :factor])
end

function _median(xs)
    isempty(xs) && return missing
    s = sort(xs)
    n = length(s)
    isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

# OLS regression of log10(y) ~ slope * log10(x) + intercept.
# Returns (slope, intercept, r2); all missing if fewer than 2 valid points.
function regression_log10(x_vals::AbstractVector, y_vals::AbstractVector)
    mask = (x_vals .> 0) .& (y_vals .> 0) .& .!ismissing.(x_vals) .& .!ismissing.(y_vals)
    xs = log10.(x_vals[mask])
    ys = log10.(y_vals[mask])
    n  = length(xs)
    n < 2 && return (slope = missing, intercept = missing, r2 = missing)
    x_mean = sum(xs) / n
    y_mean = sum(ys) / n
    ss_xx  = sum((xs .- x_mean) .^ 2)
    ss_xy  = sum((xs .- x_mean) .* (ys .- y_mean))
    ss_yy  = sum((ys .- y_mean) .^ 2)
    ss_xx == 0.0 && return (slope = missing, intercept = missing, r2 = missing)
    slope     = ss_xy / ss_xx
    intercept = y_mean - slope * x_mean
    r2        = ss_yy == 0.0 ? missing : ss_xy ^ 2 / (ss_xx * ss_yy)
    return (slope = slope, intercept = intercept, r2 = r2)
end

function score_group(label, group)
    n = nrow(group)
    if n == 0
        return (
            factor = label,
            common_isotopes = 0,
            finite_ratios = 0,
            mean_log10 = missing,
            std_log10 = missing,
            rms_log10 = missing,
            median_abs_log10 = missing,
            within_10_percent = 0,
            within_10_percent_fraction = missing,
            within_factor_2 = 0,
            within_factor_2_fraction = missing,
        )
    end

    logs     = group.log10ratio
    mean_log = sum(logs) / n
    std_log  = n > 1 ? sqrt(sum((logs .- mean_log) .^ 2) / (n - 1)) : 0.0
    within10 = group.absdev .<= log10(1.1)
    within2  = group.absdev .<= log10(2.0)
    return (
        factor = label,
        common_isotopes = n,
        finite_ratios = n,
        mean_log10 = mean_log,
        std_log10 = std_log,
        rms_log10 = sqrt(sum(logs .^ 2) / n),
        median_abs_log10 = _median(group.absdev),
        within_10_percent = sum(within10),
        within_10_percent_fraction = sum(within10) / n,
        within_factor_2 = sum(within2),
        within_factor_2_fraction = sum(within2) / n,
    )
end

function score_comparison(df_long::DataFrame)
    rows = NamedTuple[score_group(string(sub.factor[1]), sub) for sub in groupby(df_long, :factor)]
    push!(rows, score_group("overall", df_long))
    return DataFrame(rows)
end

function combined_score_table(comparisons)
    tables = DataFrame[]
    for (pair_label, df_long) in comparisons
        s = score_comparison(df_long)
        s.pair = fill(String(pair_label), nrow(s))
        push!(tables, s)
    end
    return vcat(tables...)
end

overall_scores(combined::DataFrame) = filter(row -> row.factor == "overall", combined)

function _plot_ratio_parity(df_long::DataFrame, x_col::Symbol, y_col::Symbol, x_label, y_label; within_10_percent = 0.1, within_factor_2 = 2.0, show_legend = false, plot_title = "")
    factors = sort(unique(df_long.factor); rev = true)
    styles = reaction_styles(df_long.reaction)
    n_factors = length(factors)

    panels = map(enumerate(factors)) do (idx, factor)
        last_panel = show_legend && (idx == n_factors)
        sub = filter(row -> row.factor == factor && row[x_col] > 0 && row[y_col] > 0, df_long)
        x_vals = sub[!, x_col]
        y_vals = sub[!, y_col]
        lo = minimum(vcat(x_vals, y_vals)) / 1.5
        hi = maximum(vcat(x_vals, y_vals)) * 1.5
        xs = [lo, hi]

        p = plot(
            xscale = :log10,
            yscale = :log10,
            xlims = (lo, hi),
            ylims = (lo, hi),
            xlabel = x_label,
            ylabel = y_label,
            title = "Factor = $(factor)",
            legend = last_panel ? :outertopright : false,
            aspect_ratio = :equal,
        )

        plot!(p, xs, xs .* within_factor_2; fillrange = xs ./ within_factor_2, fillalpha = 0.18, linealpha = 0, color = :gray, label = "")
        plot!(p, xs, xs .* (1 + within_10_percent); fillrange = xs .* (1 - within_10_percent), fillalpha = 0.25, linealpha = 0, color = :seagreen, label = "")
        plot!(p, xs, xs; linestyle = :dash, color = :black, label = "")

        scatter_reactions!(p, sub, x_vals, y_vals; styles = styles, show_legend = last_panel)
        p
    end

    cols = min(n_factors, 3)
    rows_count = cld(n_factors, cols)
    if isempty(plot_title)
        return plot(panels...; layout = (rows_count, cols), size = (500 * cols, 500 * rows_count))
    else
        return plot(panels...; layout = (rows_count, cols), size = (500 * cols, 500 * rows_count),
                    plot_title = plot_title, plot_titlefontsize = 11)
    end
end

function plot_ppn_vs_iliadis(df_long::DataFrame; kwargs...)
    return _plot_ratio_parity(df_long, :iliadis, :ppn, "Iliadis (2002) ratio", "PPN ratio"; kwargs...)
end

function compare_ppn_rate_sets(df_a::DataFrame, df_b::DataFrame; label_a = "rate_set_a", label_b = "rate_set_b")
    factor_cols_a = [c for c in names(df_a) if !(c in ("reaction", "isotope"))]
    long_a = stack(df_a, factor_cols_a; variable_name = :factor_col, value_name = :value_a)
    long_a.factor = _factor_value_from_ppn_column.(long_a.factor_col)
    long_a.isotope = uppercase.(string.(long_a.isotope))
    long_a = select(long_a, :reaction, :isotope, :factor, :value_a)

    factor_cols_b = [c for c in names(df_b) if !(c in ("reaction", "isotope"))]
    long_b = stack(df_b, factor_cols_b; variable_name = :factor_col, value_name = :value_b)
    long_b.factor = _factor_value_from_ppn_column.(long_b.factor_col)
    long_b.isotope = uppercase.(string.(long_b.isotope))
    long_b = select(long_b, :reaction, :isotope, :factor, :value_b)

    merged = innerjoin(long_a, long_b, on = [:reaction, :isotope, :factor])
    merged = filter(row -> !ismissing(row.value_a) && !ismissing(row.value_b), merged)

    merged.ratio = merged.value_b ./ merged.value_a
    merged.log10ratio = log10.(merged.ratio)
    merged.absdev = abs.(merged.log10ratio)
    merged.label_a = fill(label_a, nrow(merged))
    merged.label_b = fill(label_b, nrow(merged))

    return sort(merged, [:reaction, :isotope, :factor])
end

function plot_rate_set_comparison(df_long::DataFrame; label_a = "rate_set_a", label_b = "rate_set_b", kwargs...)
    return _plot_ratio_parity(df_long, :value_a, :value_b, "$(label_a) ratio", "$(label_b) ratio"; kwargs...)
end

# 6-panel parity plot: x = baseline abundance, y = factored abundance (log-log, same format as run_rate_set_comparison).
# Requires df_ref (a DataFrame with :isotope and :X columns for the reference/baseline run) so that
# absolute abundances can be back-computed from the ratio table: X_factor = ratio × X_ref.
# When df_ref is omitted, falls back to the original log₁₀(ratio) vs isotope-index scatter.
function plot_sensitivity_table(df_wide::DataFrame; df_ref = nothing, show_legend = false, plot_title = "PPN sensitivity ratios")
    factor_cols = [c for c in names(df_wide) if !(c in ("reaction", "isotope"))]
    isempty(factor_cols) && error("no factor columns found in df_wide")

    df_long = stack(df_wide, factor_cols; variable_name = :factor_col, value_name = :ratio)
    df_long.factor = _factor_value_from_ppn_column.(df_long.factor_col)
    filter!(row -> !ismissing(row.ratio) && row.ratio > 0, df_long)

    if df_ref !== nothing
        # Parity format: back-compute absolute abundances from ratios and the reference file.
        df_ref_clean = select(df_ref, :isotope, :X => :X_ref)
        df_long = innerjoin(df_long, df_ref_clean, on = :isotope)
        filter!(row -> row.X_ref > 0, df_long)
        df_long.baseline = df_long.X_ref
        df_long.value    = df_long.ratio .* df_long.X_ref
        return _plot_ratio_parity(df_long, :baseline, :value, "baseline (decay)", "factored (decay)";
                                   show_legend = show_legend, plot_title = plot_title)
    end

    # Fallback: log₁₀(ratio) vs isotope-index scatter (original format, no df_ref needed).
    df_long.log10_ratio = log10.(df_long.ratio)
    isotopes = sort(unique(df_long.isotope); by = iso -> (isotope_z(iso), iso))
    iso_idx  = Dict(iso => Float64(i) for (i, iso) in enumerate(isotopes))
    df_long.x = [iso_idx[iso] for iso in df_long.isotope]
    styles = reaction_styles(df_long.reaction)

    panels = map(sort(unique(df_long.factor); rev = true)) do f
        sub = filter(row -> row.factor == f, df_long)
        isempty(sub) && return plot(title = "factor × $f")

        yvals = sub.log10_ratio
        p = plot(
            xticks    = (1:length(isotopes), isotopes),
            xrotation = 55,
            ylabel    = "log₁₀(X_factor / X_baseline)",
            title     = "factor × $f",
            ylims     = (min(-0.15, minimum(yvals) - 0.1), max(0.15, maximum(yvals) + 0.1)),
            legend    = false,
            bottom_margin = 8Plots.mm,
        )
        hline!(p, [0.0]; linestyle = :dash, color = :black, linewidth = 1, label = "")
        scatter_reactions!(p, sub, sub.x, yvals; styles = styles)
        p
    end

    cols = min(length(panels), 3)
    rows = cld(length(panels), cols)
    return plot(panels...; layout = (rows, cols), size = (560 * cols, 450 * rows),
                plot_title = plot_title, plot_titlefontsize = 11)
end

# Return long-format rows where |X_factor/X_baseline − 1| > threshold.
function sensitive_reactions(df_wide::DataFrame; threshold = 0.1)
    factor_cols = [c for c in names(df_wide) if !(c in ("reaction", "isotope"))]
    isempty(factor_cols) && return DataFrame(reaction = String[], isotope = String[], factor = Float64[], ratio = Float64[], deviation = Float64[])
    df_long = stack(df_wide, factor_cols; variable_name = :factor_col, value_name = :ratio)
    df_long.factor = _factor_value_from_ppn_column.(df_long.factor_col)
    filter!(row -> !ismissing(row.ratio) && row.ratio > 0, df_long)
    df_long.deviation = abs.(df_long.ratio .- 1.0)
    flagged = filter(row -> row.deviation > threshold, df_long)
    isempty(flagged) && return DataFrame(reaction = String[], isotope = String[], factor = Float64[], ratio = Float64[], deviation = Float64[])
    return sort(select(flagged, :reaction, :isotope, :factor, :ratio, :deviation), [:reaction, :isotope, :factor])
end

# Per-reaction summary.
# min_effective_factor: smallest ×N (N ≥ 1) at which the threshold is first exceeded.
# max_effective_ratio:  worst max(X_f/X_base, X_base/X_f) across all (isotope, factor) pairs.
#   Always ≥ 1; "×150" means the abundance changes by a factor of 150 in either direction.
#   Equivalent to what Iliadis reports directly in Table 8.
function sensitivity_summary(df_wide::DataFrame; threshold = 0.1)
    flagged = sensitive_reactions(df_wide; threshold = threshold)
    if isempty(flagged)
        return DataFrame(
            reaction = String[], n_isotopes = Int[],
            min_effective_factor = Float64[], max_effective_ratio = Float64[],
            worst_isotope = String[], affected_isotopes = String[],
        )
    end
    rows = NamedTuple[]
    for rxn_df in groupby(flagged, :reaction)
        rxn         = first(rxn_df.reaction)
        isotopes    = sort(unique(rxn_df.isotope))
        eff_ratios  = max.(rxn_df.ratio, 1.0 ./ rxn_df.ratio)
        eff_factors = max.(rxn_df.factor, 1.0 ./ rxn_df.factor)
        worst_idx   = argmax(eff_ratios)
        push!(rows, (
            reaction             = rxn,
            n_isotopes           = length(isotopes),
            min_effective_factor = round(minimum(eff_factors); digits = 1),
            max_effective_ratio  = round(maximum(eff_ratios); digits = 2),
            worst_isotope        = rxn_df.isotope[worst_idx],
            affected_isotopes    = join(isotopes, ", "),
        ))
    end
    return sort(DataFrame(rows), :min_effective_factor)
end

function _markdown_cell(value, digits)
    ismissing(value) && return ""
    value isa AbstractFloat && return string(round(value, digits = digits))
    return string(value)
end

function dataframe_to_markdown(df::DataFrame; digits = 4)
    cols = names(df)
    header = "| " * join(cols, " | ") * " |"
    separator = "| " * join(fill("---", length(cols)), " | ") * " |"
    body = ["| " * join((_markdown_cell(row[c], digits) for c in cols), " | ") * " |" for row in eachrow(df)]
    return join(vcat(header, separator, body), "\n")
end

struct RenderedHTML
    html::String
end

Base.show(io::IO, ::MIME"text/html", x::RenderedHTML) = print(io, x.html)
Base.show(io::IO, ::MIME"text/plain", x::RenderedHTML) = print(io, x.html)

function _html_escape(value)
    s = replace(string(value), "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    return replace(s, ">" => "&gt;")
end

function _html_cell(value, digits)
    ismissing(value) && return ""
    value isa AbstractFloat && return _html_escape(round(value, digits = digits))
    return _html_escape(value)
end

function dataframe_to_html(df::DataFrame; digits = 4)
    cols = names(df)
    header = "<tr>" * join(("<th>$(_html_escape(c))</th>" for c in cols)) * "</tr>"
    body = join(
        ("<tr>" * join(("<td>$(_html_cell(row[c], digits))</td>" for c in cols)) * "</tr>" for row in eachrow(df)),
    )
    return RenderedHTML("<table>$header$body</table>")
end

summary_markdown_table(df::DataFrame; digits = 4) = dataframe_to_markdown(df; digits = digits)

function run_iliadis_comparison(; nova = pwd(), model = "JCH1", decay_runs_dir = nothing, write_results = true)
    paths = nova_case_paths(nova)
    reaction_plan = load_reaction_plan(paths.config_path)
    if decay_runs_dir !== nothing
        df_ppn = build_ppn_table8(decay_runs_dir, reaction_plan; iso_massf_filename = "iso_massfdecay.DAT")
    else
        df_ppn = build_ppn_table8(paths.runs_dir, reaction_plan)
    end
    df_iliadis = read_iliadis_table8(paths, model)
    df_table3 = read_iliadis_table3(paths)

    comparison = compare_to_iliadis(df_ppn, df_iliadis; df_table3 = df_table3)
    score = score_comparison(comparison)
    figure = plot_ppn_vs_iliadis(comparison)
    markdown_table = summary_markdown_table(score)

    if write_results
        mkpath(paths.results_dir)
        CSV.write(joinpath(paths.results_dir, "ppn_vs_iliadis_$(model)_comparison.csv"), comparison)
        CSV.write(joinpath(paths.results_dir, "ppn_vs_iliadis_$(model)_score.csv"), score)
        write(joinpath(paths.results_dir, "ppn_vs_iliadis_$(model)_score.md"), markdown_table)
        savefig(figure, joinpath(paths.results_dir, "ppn_vs_iliadis_$(model)_factor_summary.png"))
    end

    return (comparison = comparison, score = score, figure = figure, markdown_table = markdown_table, paths = paths)
end

# Compute and plot X_factor / X_reference for every (reaction, isotope, factor).
# reference_path: any iso_massf*.DAT file used as the denominator across all factors.
#   Defaults to decay_runs_dir/baseline/iso_massfdecay.DAT when omitted.
# label: used as the plot title and output file stem.
function run_sensitivity(;
    nova            = pwd(),
    decay_runs_dir,
    reference_path  = nothing,
    label           = "sensitivity",
    write_results   = true,
)
    paths         = nova_case_paths(nova)
    reaction_plan = load_reaction_plan(paths.config_path)

    reaction_isotopes = [(r.name, r.isotopes) for r in reaction_plan]
    all_factors       = sort(unique(vcat([r.factors for r in reaction_plan]...)))

    ref = something(reference_path, joinpath(decay_runs_dir, "baseline", "iso_massfdecay.DAT"))
    isfile(ref) || error("reference file not found: $ref")
    df_ref = read_iso_massf(ref)

    df_table = output_sensitivity_table(
        reaction_isotopes;
        reaction_run_path  = decay_runs_dir,
        factors            = all_factors,
        df_base            = df_ref,
        iso_massf_filename = "iso_massfdecay.DAT",
        verbose            = false,
    )

    figure = plot_sensitivity_table(df_table; df_ref = df_ref, show_legend = false, plot_title = label)

    if write_results
        results_dir = joinpath(paths.analysis_dir, "results", "sensitivity_jl")
        mkpath(results_dir)
        CSV.write(joinpath(results_dir, "$(label)_table.csv"), df_table)
        savefig(figure, joinpath(results_dir, "$(label)_factor_summary.png"))
    end

    return (table = df_table, figure = figure, paths = paths)
end

# X_ppn/X_iliadis vs isotope index scatter, sorted by Z.
# Reference bands: gray = within factor 2, green = within 10%.
function plot_abundance_ratio(df_compare::DataFrame;
        title             = "Abundance ratio: PPN / Iliadis",
        within_factor_2   = 2.0,
        within_10_percent = 0.1,
    )
    cmp = sort(df_compare, :z)
    n   = nrow(cmp)

    p = plot(
        xticks        = (1:n, cmp.isotope),
        xrotation     = 55,
        ylabel        = "X_PPN / X_Iliadis",
        yscale        = :log10,
        title         = title,
        legend        = false,
        bottom_margin = 8Plots.mm,
    )
    hline!(p, [within_factor_2, 1.0 / within_factor_2]; linestyle = :dot, color = :gray,     linewidth = 1.5, label = "")
    hline!(p, [1 + within_10_percent, 1 - within_10_percent]; linestyle = :dot, color = :seagreen, linewidth = 1.5, label = "")
    hline!(p, [1.0]; linestyle = :dash, color = :black, linewidth = 1, label = "")
    scatter!(p, 1:n, cmp.ratio; markersize = 5, color = :steelblue, label = "")
    return p
end

# Per-reaction summary of what Iliadis found sensitive in Table 8, sorted by min_effective_factor.
# Only includes reactions where at least one (isotope, factor) entry deviates from 1 by > threshold.
function iliadis_sensitivity_summary(paths, model = "JCH1"; threshold = 0.1)
    df = read_iliadis_table8(paths, model)
    factor_str_vals = [
        ("factor_100", 100.0), ("factor_10", 10.0), ("factor_2", 2.0),
        ("factor_0.5", 0.5), ("factor_0.1", 0.1), ("factor_0.01", 0.01),
    ]

    rows = NamedTuple[]
    for rxn_df in groupby(df, :reaction)
        rxn       = first(rxn_df.reaction)
        rxn_label = first(rxn_df.reaction_label)

        flagged_isotopes = Set{String}()
        min_eff_factor   = Inf
        max_eff_ratio    = 0.0
        worst_iso        = ""

        for row in eachrow(rxn_df)
            for (fc, fv) in factor_str_vals
                val = row[Symbol(fc)]
                ismissing(val) && continue
                dev = abs(val - 1.0)
                dev <= threshold && continue
                push!(flagged_isotopes, row.isotope_io)
                eff_f = max(fv, 1.0 / fv)
                eff_r = max(val, 1.0 / val)
                if eff_f < min_eff_factor
                    min_eff_factor = eff_f
                end
                if eff_r > max_eff_ratio
                    max_eff_ratio = eff_r
                    worst_iso     = row.isotope_io
                end
            end
        end

        isempty(flagged_isotopes) && continue
        push!(rows, (
            reaction             = rxn,
            article_reaction     = rxn_label,
            n_isotopes           = length(flagged_isotopes),
            min_effective_factor = round(min_eff_factor; digits = 1),
            max_effective_ratio  = round(max_eff_ratio; digits = 2),
            worst_isotope        = worst_iso,
            affected_isotopes    = join(sort(collect(flagged_isotopes)), ", "),
        ))
    end

    isempty(rows) && return DataFrame(
        reaction = String[], article_reaction = String[], n_isotopes = Int[],
        min_effective_factor = Float64[], max_effective_ratio = Float64[],
        worst_isotope = String[], affected_isotopes = String[],
    )
    return sort(DataFrame(rows), :min_effective_factor)
end

# Side-by-side comparison of PPN and Iliadis sensitive reaction lists.
# Both inputs are summary DataFrames with at least :reaction, :min_effective_factor,
# :max_deviation, :affected_isotopes columns (as returned by sensitivity_summary /
# iliadis_sensitivity_summary).
function compare_sensitive_lists(df_ppn::DataFrame, df_ili::DataFrame)
    ppn_set  = Set(df_ppn.reaction)
    ili_set  = Set(df_ili.reaction)
    all_rxns = sort(collect(union(ppn_set, ili_set)))

    rows = NamedTuple[]
    for rxn in all_rxns
        in_ppn = rxn in ppn_set
        in_ili = rxn in ili_set
        status = in_ppn && in_ili ? "both" : (in_ppn ? "ppn_only" : "iliadis_only")

        ppn_r = in_ppn ? first(filter(r -> r.reaction == rxn, eachrow(df_ppn))) : nothing
        ili_r = in_ili ? first(filter(r -> r.reaction == rxn, eachrow(df_ili))) : nothing

        push!(rows, (
            reaction              = rxn,
            status                = status,
            ppn_min_factor        = in_ppn ? ppn_r.min_effective_factor  : missing,
            ili_min_factor        = in_ili ? ili_r.min_effective_factor  : missing,
            ppn_max_effective_ratio = in_ppn ? ppn_r.max_effective_ratio : missing,
            ili_max_effective_ratio = in_ili ? ili_r.max_effective_ratio : missing,
            ppn_isotopes          = in_ppn ? ppn_r.affected_isotopes     : missing,
            ili_isotopes          = in_ili ? ili_r.affected_isotopes     : missing,
        ))
    end

    isempty(rows) && return DataFrame(
        reaction = String[], status = String[],
        ppn_min_factor          = Union{Float64,Missing}[],
        ili_min_factor          = Union{Float64,Missing}[],
        ppn_max_effective_ratio = Union{Float64,Missing}[],
        ili_max_effective_ratio = Union{Float64,Missing}[],
        ppn_isotopes            = Union{String,Missing}[],
        ili_isotopes            = Union{String,Missing}[],
    )
    return sort(DataFrame(rows), [:status, :reaction])
end

# Jaccard index and F1 score for the sensitive-reaction overlap.
# Treats Iliadis as ground truth: both = TP, ppn_only = FP, iliadis_only = FN.
# Jaccard = TP / (TP + FP + FN); F1 = 2TP / (2TP + FP + FN).
function sensitivity_agreement_score(df_compare::DataFrame)
    both_n     = count(r -> r.status == "both",         eachrow(df_compare))
    ppn_only_n = count(r -> r.status == "ppn_only",     eachrow(df_compare))
    ili_only_n = count(r -> r.status == "iliadis_only", eachrow(df_compare))
    denom_j    = both_n + ppn_only_n + ili_only_n
    denom_f    = 2 * both_n + ppn_only_n + ili_only_n
    return (
        both         = both_n,
        ppn_only     = ppn_only_n,
        iliadis_only = ili_only_n,
        jaccard      = denom_j > 0 ? round(both_n / denom_j; digits = 4) : missing,
        f1           = denom_f > 0 ? round(2 * both_n / denom_f; digits = 4) : missing,
    )
end

function run_rate_set_comparison(; nova = pwd(), run_dir_a, run_dir_b, label_a = "rate_set_a", label_b = "rate_set_b", use_decay = false, write_results = true)
    paths = nova_case_paths(nova)
    reaction_plan = load_reaction_plan(paths.config_path)
    iso_fn = use_decay ? "iso_massfdecay.DAT" : nothing
    df_a = build_ppn_table8(run_dir_a, reaction_plan; iso_massf_filename = iso_fn)
    df_b = build_ppn_table8(run_dir_b, reaction_plan; iso_massf_filename = iso_fn)

    comparison = compare_ppn_rate_sets(df_a, df_b; label_a = label_a, label_b = label_b)
    score = score_comparison(comparison)
    figure = plot_rate_set_comparison(comparison; label_a = label_a, label_b = label_b)
    markdown_table = summary_markdown_table(score)

    if write_results
        results_dir = joinpath(paths.analysis_dir, "results", "rate_set_comparison_jl")
        mkpath(results_dir)
        tag = "$(label_a)_vs_$(label_b)"
        CSV.write(joinpath(results_dir, "$(tag)_comparison.csv"), comparison)
        CSV.write(joinpath(results_dir, "$(tag)_score.csv"), score)
        write(joinpath(results_dir, "$(tag)_score.md"), markdown_table)
        savefig(figure, joinpath(results_dir, "$(tag)_factor_summary.png"))
    end

    return (comparison = comparison, score = score, figure = figure, markdown_table = markdown_table, paths = paths)
end
