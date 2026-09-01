# residual_chart.jl — N/Z chart of isotopes that crossed the tracking
# threshold between two abundance sets: net creation or net destruction,
# exactly the transitions ratio_chart can't show (a ratio against a
# near-zero side is undefined, so ratio_chart hatches those cells instead of
# coloring them).

const _RESIDUAL_DESTROYED_COLOR = :firebrick
const _RESIDUAL_CREATED_COLOR = :steelblue
const _RESIDUAL_BOTH_COLOR = (:gray, 0.55)

"""
    residual_chart(ab1::Abundances, ab2::Abundances; element_limit = "Ca", tolerance = 1e-10,
                   title = "Residual Chart", figure_size = (950, 650),
                   element_label_size = 16, mass_label_size = 8) -> CM.Figure

N/Z chart flagging every isotope that crossed `tolerance` between `ab1` and
`ab2`: **destroyed** (at or above `tolerance` in `ab1`, below it in `ab2` —
firebrick, matching [`ratio_chart`](@ref)'s color for `X1 > X2`) or
**created** (below `tolerance` in `ab1`, at or above it in `ab2` — steelblue,
matching `ratio_chart`'s color for `X1 < X2`). Isotopes tracked at or above
`tolerance` on both sides are shown neutral gray — `ratio_chart` already
covers their relative magnitude, a ratio chart and a residual chart are
meant to be read together, right after one another, not as substitutes.
Isotopes below `tolerance` on both sides are left blank.
"""
function residual_chart(ab1::Abundances, ab2::Abundances; element_limit = "Ca", tolerance = 1e-10,
                         title = "Residual Chart", figure_size = (950, 650),
                         element_label_size = 16, mass_label_size = 8)
    max_z = proton_number(element_limit)
    max_z === nothing && throw(ArgumentError("unknown element_limit \"$element_limit\""))

    isos = filter(iso -> 1 <= iso.Z <= max_z, unique(vcat(isotopes(ab1), isotopes(ab2))))
    isempty(isos) && throw(ArgumentError("no isotopes at or below element_limit=$element_limit"))

    created = NamedTuple[]
    destroyed = NamedTuple[]
    both = NamedTuple[]
    neither = NamedTuple[]
    for iso in isos
        present1, present2 = ab1[iso] >= tolerance, ab2[iso] >= tolerance
        row = (N = neutron_number(iso), Z = iso.Z, A = iso.A)
        if present1 && present2
            push!(both, row)
        elseif present1 && !present2
            push!(destroyed, row)
        elseif !present1 && present2
            push!(created, row)
        else
            push!(neither, row)
        end
    end
    all_n = vcat([r.N for r in created], [r.N for r in destroyed], [r.N for r in both], [r.N for r in neither])
    min_n, max_n = extrema(all_n)

    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1];
            xlabel = "neutron number (A-Z)", ylabel = "proton number (Z)", title,
            aspect = CM.DataAspect(), xgridvisible = false, ygridvisible = false,
            xticks = min_n:max_n, yticks = 0:2:max_z,
            limits = (min_n - 3.0, max_n + 1.0, -0.5, max_z + 1.0))

        _draw_residual_tiles!(ax, destroyed, _RESIDUAL_DESTROYED_COLOR, mass_label_size)
        _draw_residual_tiles!(ax, created, _RESIDUAL_CREATED_COLOR, mass_label_size)
        _draw_residual_tiles!(ax, both, _RESIDUAL_BOTH_COLOR, mass_label_size)
        _draw_residual_tiles!(ax, neither, (:gray, 0.15), mass_label_size)

        elem_df = DataFrame(N = all_n, Z = vcat([r.Z for r in created], [r.Z for r in destroyed],
                                                  [r.Z for r in both], [r.Z for r in neither]))
        add_element_labels!(ax, elem_df, min_n, max_z; element_label_size)

        CM.Legend(fig[1, 2],
            [CM.PolyElement(color = _RESIDUAL_DESTROYED_COLOR), CM.PolyElement(color = _RESIDUAL_CREATED_COLOR),
             CM.PolyElement(color = _RESIDUAL_BOTH_COLOR), CM.PolyElement(color = (:gray, 0.15))],
            ["destroyed (X1→below tol.)", "created (below tol.→X2)",
             "tracked in both", "tracked in neither"])
        fig
    end
end

function _draw_residual_tiles!(ax, rows, color, mass_label_size)
    for r in rows
        CM.poly!(ax, _tile_corners(r.N, r.Z); color, strokecolor = :black, strokewidth = 1)
        CM.text!(ax, string(r.A); position = (r.N, r.Z), align = (:center, :center),
                  fontsize = mass_label_size, color = :black)
    end
end
