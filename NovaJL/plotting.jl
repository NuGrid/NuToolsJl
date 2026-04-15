module Plotting
    using ..NovaPPN: AbundanceSet
    using DataFrames
    using Plots
    using CairoMakie
    using LaTeXStrings

    export analyze_factor, plot_trajectory, abundance_chart, reaction_stylepoint_table

    const REACTION_COLORS = [
        :blue, :orange, :green, :red, :purple, :brown, :pink, :gray,
        :olive, :cyan, :magenta, :gold, :navy, :teal, :coral, :black
    ]

    const REACTION_MARKERS = [
        :circle, :rect, :diamond, :utriangle, :dtriangle, :star5,
        :hexagon, :cross, :xcross
    ]

    function _read_reaction_plan_entries(config_path::AbstractString)
        if !isfile(config_path)
            candidate_paths = [
                joinpath("..", config_path),
                joinpath("..", "..", config_path),
                joinpath("config", "reaction_plan.json"),
                joinpath("..", "config", "reaction_plan.json"),
            ]

            found = findfirst(isfile, candidate_paths)
            if found !== nothing
                config_path = candidate_paths[found]
            else
                error("Reaction plan not found: $config_path")
            end
        end

        text = read(config_path, String)
        entries = Vector{NamedTuple{(:reaction, :index), Tuple{String, Union{Missing, Int}}}}()

        for m in eachmatch(r"\{[^{}]*\"name\"\s*:\s*\"([^\"]+)\"[^{}]*\}"s, text)
            block = m.match
            name = m.captures[1]
            index_match = match(r"\"index\"\s*:\s*([0-9]+)", block)
            index = index_match === nothing ? missing : parse(Int, index_match.captures[1])
            push!(entries, (reaction=name, index=index))
        end

        return entries
    end

    function _stylepoint_table_from_reactions(reactions; indices=nothing)
        markers = Symbol[]
        colors = Symbol[]
        style_ids = Int[]

        for i in eachindex(reactions)
            color_i = REACTION_COLORS[mod1(i, length(REACTION_COLORS))]
            marker_i = REACTION_MARKERS[mod1(i, length(REACTION_MARKERS))]
            push!(colors, color_i)
            push!(markers, marker_i)
            push!(style_ids, i)
        end

        index_col = indices === nothing ? fill(missing, length(reactions)) : indices

        return DataFrame(
            style_id = style_ids,
            reaction = collect(reactions),
            index = index_col,
            color = colors,
            marker = markers
        )
    end

    function reaction_stylepoint_table(config_path::AbstractString="novae/nova_test/config/reaction_plan.json")
        entries = _read_reaction_plan_entries(config_path)
        reactions = [entry.reaction for entry in entries]
        indices = [entry.index for entry in entries]
        return _stylepoint_table_from_reactions(reactions; indices=indices)
    end

    function _style_lookup(style_table)
        lookup = Dict{String, NamedTuple{(:color, :marker), Tuple{Symbol, Symbol}}}()

        for row in eachrow(style_table)
            lookup[String(row.reaction)] = (color=Symbol(row.color), marker=Symbol(row.marker))
        end

        return lookup
    end

    function _offset_by_isotope_and_reaction(df, isotope_index)
        offsets = zeros(Float64, nrow(df))

        for iso in unique(df.isotope)
            rows = findall(==(iso), df.isotope)
            reactions = hasproperty(df, :reaction) ? unique(df.reaction[rows]) : [""]
            reaction_rank = Dict(reaction => i for (i, reaction) in enumerate(reactions))
            n = length(reactions)
            step = n <= 1 ? 0.0 : min(0.12, 0.45 / max(n - 1, 1))

            for row_i in rows
                rank = reaction_rank[hasproperty(df, :reaction) ? df.reaction[row_i] : ""]
                offsets[row_i] = (rank - (n + 1) / 2) * step
            end
        end

        return [isotope_index[iso] + offsets[i] for (i, iso) in enumerate(df.isotope)]
    end

    function _scatter_reactions!(plot_obj, df, xcol, ycol, style_table; show_reaction_legend=false)
        lookup = _style_lookup(style_table)

        for reaction in unique(df.reaction)
            rows = findall(==(reaction), df.reaction)
            style = get(lookup, String(reaction), (color=:black, marker=:circle))

            Plots.scatter!(
                plot_obj,
                df[rows, xcol],
                df[rows, ycol],
                color = style.color,
                marker = style.marker,
                label = show_reaction_legend ? String(reaction) : "",
            )
        end

        return plot_obj
    end

    function analyze_factor(df_compare, f; style_table=nothing, config_path=nothing, show_reaction_legend=false, return_style_table=false)
        
        # -------------------------
        # COLUMN NAMES
        # -------------------------
        col_ppn = Symbol(f * "_ppn")
        col_iliadis = Symbol(f * "_iliadis")
        col_ratio = Symbol(f * "_ratio")
    
        # -------------------------
        # FILTER VALID ROWS
        # -------------------------
        df_valid = filter(row ->
            !ismissing(row[col_ppn]) &&
            !ismissing(row[col_iliadis]),
            df_compare
        )

        # -------------------------
        # COMPUTE METRICS
        # -------------------------
        df_valid.ratio = df_valid[!, col_ppn] ./ df_valid[!, col_iliadis]
    
        df_valid.logratio = log10.(df_valid.ratio)
        df_valid.dev = abs.(df_valid.logratio)
        has_reaction = hasproperty(df_valid, :reaction)

        if has_reaction
            df_valid.reaction = String.(df_valid.reaction)

            if style_table === nothing
                if config_path !== nothing
                    style_table = reaction_stylepoint_table(config_path)
                else
                    reactions = unique(df_valid.reaction)
                    style_table = _stylepoint_table_from_reactions(reactions)
                end
            end
        end

        # -------------------------
        # SCATTER: PPN vs Iliadis
        # -------------------------
        p1 = Plots.plot(
            xlabel = "Iliadis",
            ylabel = "PPN",
            title = "PPN vs Iliadis (factor = $f)",
            label = "",
            xscale = :log10,
            yscale = :log10
        )

        if has_reaction
            _scatter_reactions!(p1, df_valid, col_iliadis, col_ppn, style_table; show_reaction_legend=show_reaction_legend)
        else
            Plots.scatter!(p1, df_valid[!, col_iliadis], df_valid[!, col_ppn], label = "")
        end

        Plots.plot!(p1, [1e-3, 10], [1e-3, 10], linestyle=:dash, label="y = x")

        # -------------------------
        # GROUP BY ISOTOPE
        # -------------------------
        isotopes = unique(df_valid.isotope)
        iso_index = Dict(iso => i for (i, iso) in enumerate(isotopes))
        x = has_reaction ? _offset_by_isotope_and_reaction(df_valid, iso_index) : [iso_index[iso] for iso in df_valid.isotope]
        df_valid.plot_x = x

        # -------------------------
        # LOG RATIO PLOT
        # -------------------------
        p2 = Plots.plot(
            xticks = (1:length(isotopes), isotopes),
            xlabel = "Isotope",
            ylabel = "log10(PPN / Iliadis)",
            title = "Log Ratio (factor = $f)",
            label = "",
            xrotation = 60,
            size = (900, 400),
            ylims = (-1, 3)
        )

        if has_reaction
            _scatter_reactions!(p2, df_valid, :plot_x, :logratio, style_table; show_reaction_legend=show_reaction_legend)
        else
            Plots.scatter!(p2, x, df_valid.logratio, label = "")
        end

        Plots.hline!(p2, [0.0], linestyle=:dash)
    
        # -------------------------
        # DEVIATION PLOT
        # -------------------------
        p3 = Plots.plot(
            xticks = (1:length(isotopes), isotopes),
            xlabel = "Isotope",
            ylabel = "|log10(PPN / Iliadis)|",
            title = "Deviation (factor = $f)",
            label = "",
            xrotation = 60,
            size = (900, 400)
        )

        if has_reaction
            _scatter_reactions!(p3, df_valid, :plot_x, :dev, style_table; show_reaction_legend=show_reaction_legend)
        else
            Plots.scatter!(p3, x, df_valid.dev, label = "")
        end

        return return_style_table && has_reaction ? (p1, p2, p3, df_valid, style_table) : (p1, p2, p3, df_valid)
    end

    # Dependent on the ELEMENT_Z found in utils, if you want to see more isotopes
    # please go and extend the list of isotopes you would like to see. 
    # Since this is a nova analysis tool, we have only implemented isotopes occurent in nova
    # which is usually up to Ca
    function abundance_chart(ab::AbundanceSet;
                logscale=true,
                zmax=nothing,
                title="Abundance Chart")
    
        # Extract arrays
        Z = [iso.Z for iso in ab.isotopes]
        N = [iso.N for iso in ab.isotopes]
        X = [iso.X for iso in ab.isotopes]
        names = [iso.name for iso in ab.isotopes]

        valid = [(Z[i] != 0 || N[i] != 0) && X[i] > 0 for i in eachindex(Z)]
        Z = Z[valid]
        N = N[valid]
        X = X[valid]
        names = names[valid]

        # Optional Z cutoff
        if zmax !== nothing
            mask = Z .<= zmax
            Z, N, X, names = Z[mask], N[mask], X[mask], names[mask]
        end

        # turn into log scale
        values = logscale ? log10.(X .+ 1e-30) : X

        fig = CairoMakie.Figure(size=(1000, 800))
        ax = CairoMakie.Axis(fig[1,1],
                  xlabel="Neutron number (N)",
                  ylabel="Atomic number (Z)",
                  title=title,
                  xgridvisible=false,
                  ygridvisible=false
                )
    
        rects = Rect[]
        colors = Float64[]
        size = 1.0
        for i in eachindex(names)
            push!(rects, Rect(N[i]-size/2, Z[i]-size/2, size, size))
            push!(colors, values[i])
        end

        p = CairoMakie.poly!(ax, rects,
                  color=colors,
                  colormap=:linear_wyor_100_45_c55_n256,
                  strokecolor=:black,
                  strokewidth=1.0
                )

        # add isotope labels
        for i in eachindex(names)
            CairoMakie.text!(ax, N[i], Z[i],
                  text=names[i],
                  align=(:center, :center),
                  fontsize=8,
                  font=:bold,
                  color=:black)
        end

        CairoMakie.Colorbar(fig[1,2], p, label="log10(X)")
        return fig
    end

    function plot_trajectory(df)

        fig = CairoMakie.Figure(size=(1000, 600))

        ax1 = CairoMakie.Axis(fig[1,1],
            xlabel = "Time (minutes)",
            ylabel = "Temperature (GK)",
            yticklabelcolor = :red,
            ylabelcolor = :red
        )

        ax2 = CairoMakie.Axis(fig[1,1],
            yaxisposition = :right,
            ylabel = "Density (g/cm⁻³)",
            yticklabelcolor = :blue,
            ylabelcolor = :blue,
            yscale = log10
        )

        # Temperature (left axis)
        CairoMakie.lines!(ax1, df.time, df.T, color=:red, label="Temperature")

        # Density (right axis)
        CairoMakie.lines!(ax2, df.time, df.rho, color=:blue, label="Density")

        # link x-axis
        CairoMakie.linkxaxes!(ax1, ax2)

        imax = argmax(df.T)
    
        CairoMakie.scatter!(ax1,
            [df.time[imax]],
            [df.T[imax]],
            color = :black,
            markersize = 10
        )

        return fig
    end

end # module
