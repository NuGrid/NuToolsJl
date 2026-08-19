# ratio_chart.jl — nuclear-chart comparison of two abundance sets.

"""
    ratio_chart(ab1::Abundances, ab2::Abundances; element_limit = "Ca", tolerance = 1e-10,
                color_range = 2.0, title = "Abundance ratio chart", figure_size = (950, 650),
                element_label_size = 16, mass_label_size = 8,
                colorbar_label = "log₁₀(ratio)") -> CM.Figure

N/Z chart of `X(ab1) / X(ab2)`, colored on the diverging
[`RATIO_COLORMAP`](@ref) (steelblue below 1, firebrick above 1), symmetrically
clipped at `±color_range` dex. Every isotope tracked by `ab1` or `ab2` within
`element_limit` gets a tile and a mass-number label; isotopes where the ratio
isn't well-defined (below `tolerance` in either side) are drawn hatched
(gray) instead of colored, rather than dropped or pinned to an extreme
color, since a ratio against a near-zero denominator isn't meaningful — a
blank patch always means "tracked by neither," never "below threshold." See
[`changed_isotopes`](@ref) for the ranked-list companion to this chart.
"""
function ratio_chart(ab1::Abundances, ab2::Abundances; element_limit = "Ca", tolerance = 1e-10,
                      color_range = 2.0, title = "Abundance ratio chart", figure_size = (950, 650),
                      element_label_size = 16, mass_label_size = 8,
                      colorbar_label = "log₁₀(ratio)")
    max_z = proton_number(element_limit)
    max_z === nothing && throw(ArgumentError("unknown element_limit \"$element_limit\""))

    isos = filter(iso -> 1 <= iso.Z <= max_z, unique(vcat(isotopes(ab1), isotopes(ab2))))
    isempty(isos) && throw(ArgumentError("no isotopes at or below element_limit=$element_limit"))

    rows = NamedTuple[]
    hatched = NamedTuple[]
    for iso in isos
        x1, x2 = ab1[iso], ab2[iso]
        if x1 >= tolerance && x2 >= tolerance
            push!(rows, (N = neutron_number(iso), Z = iso.Z, A = iso.A, log_ratio = log10(x1 / x2)))
        else
            push!(hatched, (N = neutron_number(iso), Z = iso.Z, A = iso.A))
        end
    end
    all_n = vcat([r.N for r in rows], [r.N for r in hatched])
    min_n, max_n = extrema(all_n)

    return with_nugrid_theme() do
        fig = CM.Figure(size = figure_size)
        ax = CM.Axis(fig[1, 1];
            xlabel = "neutron number (A-Z)", ylabel = "proton number (Z)", title,
            aspect = CM.DataAspect(), xgridvisible = false, ygridvisible = false,
            xticks = min_n:max_n, yticks = 0:2:max_z,
            limits = (min_n - 3.0, max_n + 1.0, -0.5, max_z + 1.0))

        for r in rows
            CM.poly!(ax, _tile_corners(r.N, r.Z);
                      color = chart_colormap_color(r.log_ratio, (-color_range, color_range), RATIO_COLORMAP),
                      strokecolor = :black, strokewidth = 1)
            CM.text!(ax, string(r.A); position = (r.N, r.Z), align = (:center, :center),
                      fontsize = mass_label_size, color = :black)
        end
        for r in hatched
            CM.poly!(ax, _tile_corners(r.N, r.Z); color = (:gray, 0.3), strokecolor = :black, strokewidth = 1)
            CM.text!(ax, string(r.A); position = (r.N, r.Z), align = (:center, :center),
                      fontsize = mass_label_size, color = :black)
        end

        elem_df = DataFrame(N = all_n, Z = vcat([r.Z for r in rows], [r.Z for r in hatched]))
        add_element_labels!(ax, elem_df, min_n, max_z; element_label_size)

        CM.Colorbar(fig[1, 2]; colormap = RATIO_COLORMAP, colorrange = (-color_range, color_range), label = colorbar_label)
        fig
    end
end

"""
    changed_isotopes(ab1::Abundances, ab2::Abundances; threshold = 0.1, tolerance = 1e-10) -> DataFrame

Isotopes present in both `ab1` and `ab2` (each above `tolerance`) whose ratio
`X(ab1)/X(ab2)` differs from 1 by at least `threshold` in log10, ranked by the
size of that change — the "what moved the most" companion to
[`ratio_chart`](@ref), for spotting an isotope worth investigating with
[`reactions_for_isotope`](@ref).
"""
function changed_isotopes(ab1::Abundances, ab2::Abundances; threshold = 0.1, tolerance = 1e-10)
    rows = NamedTuple[]
    for iso in intersect(isotopes(ab1), isotopes(ab2))
        x1, x2 = ab1[iso], ab2[iso]
        (x1 < tolerance || x2 < tolerance) && continue
        log_ratio = log10(x1 / x2)
        abs(log_ratio) >= threshold || continue
        push!(rows, (isotope = isotope_name(iso), Z = iso.Z, A = iso.A, X1 = x1, X2 = x2,
                      ratio = x1 / x2, log_ratio = log_ratio))
    end
    isempty(rows) && return DataFrame(isotope = String[], Z = Int[], A = Int[], X1 = Float64[],
                                        X2 = Float64[], ratio = Float64[], log_ratio = Float64[])
    return sort(DataFrame(rows), :log_ratio; by = abs, rev = true)
end
