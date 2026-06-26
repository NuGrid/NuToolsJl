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

function nova_case_paths(case_analysis_dir = pwd())
    analysis_dir = abspath(case_analysis_dir)
    case_dir = dirname(analysis_dir)
    nova_cases_dir = dirname(case_dir)
    sensitivity_root = dirname(nova_cases_dir)
    data_root = joinpath(sensitivity_root, "data", "iliadis2002")
    return (
        analysis_dir = analysis_dir,
        case_dir = case_dir,
        runs_dir = joinpath(case_dir, "runs"),
        decay_time_scan_csv = joinpath(analysis_dir, "results", "baseline_decay_checker", "decay_time_scan.csv"),
        config_path = joinpath(case_dir, "config", "reaction_plan.json"),
        data_root = data_root,
        results_dir = joinpath(analysis_dir, "results", "iliadis_comparison_jl"),
    )
end

function latest_validated_decay_run(paths)
    isfile(paths.decay_time_scan_csv) || error(
        "missing $(paths.decay_time_scan_csv); run `julia tools/decay_time_scan.jl --nova <case>` " *
        "then `python3 baseline_decay_checker.py` first",
    )
    df = CSV.read(paths.decay_time_scan_csv, DataFrame)
    isempty(df) && error("no rows in $(paths.decay_time_scan_csv)")
    "output" in names(df) || error(
        "no `output` column in $(paths.decay_time_scan_csv) — this file predates the decay-run-path " *
        "fix in baseline_decay_checker.py; re-run `python3 baseline_decay_checker.py` to regenerate it",
    )
    return String(df.output[1])
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

function build_ppn_table8(run_dir::AbstractString, reaction_plan; verbose = false)
    cycle = final_cycle_label(joinpath(run_dir, "baseline", "x-time.dat"))
    iso_massf_filename = "iso_massf$(cycle).DAT"
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

function score_group(label, group)
    n = nrow(group)
    if n == 0
        return (
            factor = label,
            common_isotopes = 0,
            finite_ratios = 0,
            rms_log10 = missing,
            median_abs_log10 = missing,
            within_10_percent = 0,
            within_10_percent_fraction = missing,
            within_factor_2 = 0,
            within_factor_2_fraction = missing,
        )
    end

    within10 = group.absdev .<= log10(1.1)
    within2 = group.absdev .<= log10(2.0)
    return (
        factor = label,
        common_isotopes = n,
        finite_ratios = n,
        rms_log10 = sqrt(sum(group.log10ratio .^ 2) / n),
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

function _plot_ratio_parity(df_long::DataFrame, x_col::Symbol, y_col::Symbol, x_label, y_label; within_10_percent = 0.1, within_factor_2 = 2.0)
    factors = sort(unique(df_long.factor); rev = true)
    styles = reaction_styles(df_long.reaction)

    panels = map(factors) do factor
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
            legend = false,
            aspect_ratio = :equal,
        )

        plot!(p, xs, xs .* within_factor_2; fillrange = xs ./ within_factor_2, fillalpha = 0.18, linealpha = 0, color = :gray, label = "")
        plot!(p, xs, xs .* (1 + within_10_percent); fillrange = xs .* (1 - within_10_percent), fillalpha = 0.25, linealpha = 0, color = :seagreen, label = "")
        plot!(p, xs, xs; linestyle = :dash, color = :black, label = "")

        scatter_reactions!(p, sub, x_vals, y_vals; styles = styles)
        p
    end

    cols = min(length(panels), 3)
    rows = cld(length(panels), cols)
    return plot(panels...; layout = (rows, cols), size = (500 * cols, 500 * rows))
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

function run_iliadis_comparison(; nova = pwd(), model = "JCH1", write_results = true)
    paths = nova_case_paths(nova)
    reaction_plan = load_reaction_plan(paths.config_path)
    df_ppn = build_ppn_table8(paths.runs_dir, reaction_plan)
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

function run_rate_set_comparison(; nova = pwd(), run_dir_a, run_dir_b, label_a = "rate_set_a", label_b = "rate_set_b", write_results = true)
    paths = nova_case_paths(nova)
    reaction_plan = load_reaction_plan(paths.config_path)
    df_a = build_ppn_table8(run_dir_a, reaction_plan)
    df_b = build_ppn_table8(run_dir_b, reaction_plan)

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
