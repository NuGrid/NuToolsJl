# chart_primitives.jl — shared N/Z nuclear-chart drawing.
#
# Every chart in this package (abundance, flux, ratio) draws the same kind of
# picture: square tiles on a neutron-number/proton-number grid, optionally
# joined by flow arrows. These helpers draw that picture once so the chart
# functions only have to decide what data goes into it. Not exported — treat
# them as the package's internal drawing toolkit.

_tile_corners(n, z) = CM.Point2f[
    (n - 0.5, z - 0.5), (n + 0.5, z - 0.5),
    (n + 0.5, z + 0.5), (n - 0.5, z + 0.5),
]

"Color for a chart tile: `colormap` sampled at `value`'s fraction between `color_limits`."
function chart_colormap_color(value, color_limits, colormap)
    fraction = clamp((value - color_limits[1]) / (color_limits[2] - color_limits[1]), 0.0, 1.0)
    return CM.cgrad(colormap)[fraction]
end

"""
    draw_isotope_tiles!(ax, df; value_col, color_limits, colormap, mass_label_size = 8)

Draw one colored, mass-labeled tile per row of `df` (needs `:N`, `:Z`, `:A`
columns), colored by `log10(df[!, value_col])` scaled across `color_limits`.
"""
function draw_isotope_tiles!(ax, df::DataFrame; value_col, color_limits, colormap, mass_label_size = 8)
    for row in eachrow(df)
        value = log10(max(row[value_col], 10.0^color_limits[1]))
        CM.poly!(ax, _tile_corners(row.N, row.Z);
                  color = chart_colormap_color(value, color_limits, colormap),
                  strokecolor = :black, strokewidth = 1)
        CM.text!(ax, string(row.A); position = (row.N, row.Z), align = (:center, :center),
                  fontsize = mass_label_size, color = :black)
    end
end

"""
    draw_empty_isotope_tiles!(ax, df; mass_label_size = 8)

Draw white, mass-labeled tiles (no abundance color) for every `(N, Z, A)` row of `df`.
"""
function draw_empty_isotope_tiles!(ax, df::DataFrame; mass_label_size = 8)
    for row in eachrow(df)
        CM.poly!(ax, _tile_corners(row.N, row.Z); color = :white, strokecolor = :black, strokewidth = 1)
        CM.text!(ax, string(row.A); position = (row.N, row.Z), align = (:center, :center),
                  fontsize = mass_label_size, color = :black)
    end
end

"""
    add_element_labels!(ax, df, min_n, max_z; min_z = 1, element_label_size = 16)

Draw an element symbol for each proton number `min_z:max_z`, placed just left
of that row's leftmost tile in `df` (needs `:N`, `:Z` columns).
"""
function add_element_labels!(ax, df::DataFrame, min_n, max_z; min_z = 1, element_label_size = 16)
    for z in min_z:max_z
        rows = filter(:Z => ==(z), df)
        n_label = nrow(rows) == 0 ? min_n - 1.5 : minimum(rows.N) - 0.75
        CM.text!(ax, element_symbol(z); position = (n_label, z), align = (:right, :center),
                  fontsize = element_label_size, font = :bold, color = :black)
    end
end

"""
    draw_arrow!(ax, x1, y1, x2, y2; color, linewidth = 2, head_length = 0.22, head_width = 0.16)

A flow arrow from `(x1, y1)` to `(x2, y2)`, inset from both endpoints so it
doesn't overlap the tile borders it connects.
"""
function draw_arrow!(ax, x1, y1, x2, y2; color, linewidth = 2, head_length = 0.22, head_width = 0.16)
    dx, dy = x2 - x1, y2 - y1
    distance = hypot(dx, dy)
    distance == 0 && return
    ux, uy = dx / distance, dy / distance
    px, py = -uy, ux

    start_margin = min(0.35, 0.25 * distance)
    end_margin = min(0.45, 0.35 * distance)
    sx, sy = x1 + start_margin * ux, y1 + start_margin * uy
    tip_x, tip_y = x2 - end_margin * ux, y2 - end_margin * uy
    base_x, base_y = tip_x - head_length * ux, tip_y - head_length * uy

    CM.lines!(ax, [sx, base_x], [sy, base_y]; color, linewidth)
    CM.poly!(ax, CM.Point2f[
        (tip_x, tip_y),
        (base_x + 0.5 * head_width * px, base_y + 0.5 * head_width * py),
        (base_x - 0.5 * head_width * px, base_y - 0.5 * head_width * py),
    ]; color, strokecolor = color)
end

"""
    draw_arrow_label!(ax, x1, y1, x2, y2, label; fontsize = 8, offset = 0.18)

A small text label offset to the side of the arrow from `(x1, y1)` to `(x2,
y2)`, rotated to follow the arrow's direction.
"""
function draw_arrow_label!(ax, x1, y1, x2, y2, label; fontsize = 8, offset = 0.18)
    dx, dy = x2 - x1, y2 - y1
    distance = hypot(dx, dy)
    distance == 0 && return
    px, py = -dy / distance, dx / distance
    angle = atan(dy, dx)
    CM.text!(ax, label; position = ((x1 + x2) / 2 + offset * px, (y1 + y2) / 2 + offset * py),
              align = (:center, :center), rotation = angle, fontsize, color = :black)
end

"""
    flow_tile_data(fluxes) -> DataFrame

Every unique `(N, Z, A)` tile touched by `fluxes` (both reactant and product
endpoints of each row), for drawing the empty tile grid under a flux chart.
"""
function flow_tile_data(fluxes::DataFrame)
    tiles = Dict{Tuple{Int,Int},Int}()
    for row in eachrow(fluxes)
        tiles[(row.n_start, row.z_start)] = row.a_start
        tiles[(row.n_end, row.z_end)] = row.a_end
    end
    return DataFrame(
        N = [k[1] for k in keys(tiles)],
        Z = [k[2] for k in keys(tiles)],
        A = collect(values(tiles)),
    )
end
