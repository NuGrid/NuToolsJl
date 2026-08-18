# tables.jl — rendering and saving DataFrames: markdown, HTML, CSV.

"""
    RenderedHTML

Wraps an HTML string so it displays as rendered markup in a notebook/REPL
(`text/html`) instead of as a literal string. Returned by
[`dataframe_to_html`](@ref).
"""
struct RenderedHTML
    html::String
end

Base.show(io::IO, ::MIME"text/html", x::RenderedHTML) = print(io, x.html)
Base.show(io::IO, ::MIME"text/plain", x::RenderedHTML) = print(io, x.html)

_table_cell(value, digits) = ismissing(value) ? "" :
    value isa AbstractFloat ? string(round(value; digits)) : string(value)

"""
    dataframe_to_markdown(df::DataFrame; digits = 4) -> String

Render `df` as a Markdown pipe table, rounding floats to `digits` and blanking
`missing` cells.
"""
function dataframe_to_markdown(df::DataFrame; digits = 4)
    cols = names(df)
    header = "| " * join(cols, " | ") * " |"
    separator = "| " * join(fill("---", length(cols)), " | ") * " |"
    body = ["| " * join((_table_cell(row[c], digits) for c in cols), " | ") * " |" for row in eachrow(df)]
    return join(vcat(header, separator, body), "\n")
end

_html_escape(value) = replace(replace(replace(string(value), "&" => "&amp;"), "<" => "&lt;"), ">" => "&gt;")
_html_cell(value, digits) = ismissing(value) ? "" : _html_escape(
    value isa AbstractFloat ? round(value; digits) : value)

"""
    dataframe_to_html(df::DataFrame; digits = 4) -> RenderedHTML

Render `df` as an HTML `<table>`, rounding floats to `digits` and blanking
`missing` cells.
"""
function dataframe_to_html(df::DataFrame; digits = 4)
    cols = names(df)
    header = "<tr>" * join(("<th>$(_html_escape(c))</th>" for c in cols)) * "</tr>"
    body = join(("<tr>" * join(("<td>$(_html_cell(row[c], digits))</td>" for c in cols)) * "</tr>" for row in eachrow(df)))
    return RenderedHTML("<table>$header$body</table>")
end

"""
    save_table(df::DataFrame, path) -> DataFrame

Write `df` to `path` as CSV (creating parent directories as needed) and return
`df` unchanged, so this can be chained: `save_table(df, "out.csv") |> dataframe_to_html`.
"""
function save_table(df::DataFrame, path::AbstractString)
    mkpath(dirname(abspath(path)))
    CSV.write(path, df)
    return df
end
