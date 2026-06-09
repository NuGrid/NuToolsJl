using DataFrames

const DEFAULT_REACTION_ISOTOPES = [
    ("3He_ag_7Be", ["BE-7"]),
    ("7Be_pg_8B", ["BE-7"]),
    ("8B_pg_9C", ["BE-7"]),

    ("13N_pg_14O", ["C-12", "C-13", "O-18"]),
    ("14N_pg_15O", ["C-12", "C-13", "N-14", "N-15"]),
    ("15N_pa_12C", ["N-15"]),

    ("16O_pg_17F", ["C-13", "O-16", "O-17", "O-18", "F-18", "F-19"]),
    ("17O_pg_18F", ["O-17", "O-18", "F-18", "F-19"]),
    ("17O_pa_14N", ["O-17", "O-18", "F-18", "F-19"]),
    ("17F_pg_18Ne", ["O-17", "O-18", "F-18", "F-19"]),

    ("18F_pg_19Ne", ["F-19"]),
    ("18F_pa_15O", ["O-18", "F-18", "F-19"]),
    ("19F_pa_16O", ["F-19"]),

    ("20Ne_pg_21Na", ["NE-21", "NA-22", "NA-23", "MG-24", "MG-25", "MG-26", "AL-26", "AL-27", "SI-28", "SI-29"]),
    ("21Ne_pg_22Na", ["NE-21"]),
    ("22Ne_pg_23Na", ["NE-22"]),

    ("21Na_pg_22Mg", ["NE-21", "NA-22", "MG-25", "MG-26", "AL-26", "AL-27"]),
    ("22Na_pg_23Mg", ["NA-22"]),

    ("23Na_pg_24Mg", ["NE-20", "NE-21", "NA-22", "NA-23", "MG-24", "MG-25", "MG-26", "AL-26", "AL-27", "SI-28", "SI-29", "SI-30", "P-31", "S-32", "S-33"]),
    ("23Na_pa_20Ne", ["NA-23", "MG-24", "MG-25", "MG-26", "AL-26", "AL-27", "SI-28", "SI-29"]),

    ("25Mg_pg_26Alg", ["MG-25", "MG-26", "AL-26"]),
    ("25Mg_pg_26Alm", ["MG-26"]),
    ("26Mg_pg_27Al", ["MG-26"]),
    ("26Alg_pg_27Si", ["AL-26", "AL-27"]),
    ("26Alm_pg_27Si", ["MG-26"]),
    ("27Al_pg_28Si", ["AL-27"]),

    ("28Si_pg_29P", ["SI-28", "SI-29", "SI-30", "P-31", "S-32", "S-33"]),
    ("29Si_pg_30P", ["SI-29", "SI-30", "P-31", "S-32", "S-33"]),
    ("30Si_pg_31P", ["SI-30", "P-31"]),
    ("30P_pg_31S", ["SI-30", "P-31", "S-32", "S-33", "S-34"]),

    ("31P_pa_28Si", ["P-31", "S-32", "S-33"]),
    ("32S_pg_33Cl", ["S-33", "S-34"]),
    ("33S_pg_34Cl", ["S-33", "S-34", "CL-35"]),
    ("34S_pg_35Cl", ["S-34", "CL-35"]),

    ("36Ar_pg_37K", ["AR-36", "AR-37", "AR-38"]),
    ("37Ar_pg_38K", ["AR-37", "AR-38", "CL-37", "K-39"]),
    ("38Ar_pg_39K", ["AR-38", "K-39"]),
    ("38K_pg_39Ca", ["AR-38", "K-39"]),
    ("39K_pg_40Ca", ["K-39"])
]

const REACTION_ISOTOPES = DEFAULT_REACTION_ISOTOPES

function prepare_base_massf(df_base::DataFrame)
    prepared = copy(df_base)
    if :X_base in propertynames(prepared)
        return prepared[:, [:isotope, :X_base]]
    elseif :X in propertynames(prepared)
        rename!(prepared, :X => :X_base)
        return prepared[:, [:isotope, :X_base]]
    end

    throw(ArgumentError("df_base must contain either an :X or :X_base column"))
end

function candidate_reaction_run_roots(reaction_run_path, base_path)
    roots = [string(reaction_run_path)]

    if base_path !== nothing
        base_run_root = dirname(dirname(string(base_path)))
        push!(roots, base_run_root)
    end

    return unique(roots)
end

function factor_file_candidates(reaction_run_path, reaction_name, factor, iso_massf_filename, base_path)
    roots = candidate_reaction_run_roots(reaction_run_path, base_path)
    labels = factor_directory_labels(factor)
    candidates = String[]

    for root in roots
        for label in labels
            for prefix in ("factor_", "fact_")
                push!(candidates, joinpath(root, reaction_name, "$(prefix)$(label)", iso_massf_filename))
            end
        end
    end

    return unique(candidates)
end

function first_existing_factor_file(reaction_run_path, reaction_name, factor, iso_massf_filename, base_path)
    candidates = factor_file_candidates(reaction_run_path, reaction_name, factor, iso_massf_filename, base_path)
    for candidate in candidates
        isfile(candidate) && return candidate, candidates
    end
    return nothing, candidates
end

function output_sensitivity_table(
    reaction_isotopes = DEFAULT_REACTION_ISOTOPES;
    reaction_run_path = "../runs",
    factors = [100, 10, 2, 0.5, 0.1, 0.01],
    df_base = nothing,
    base_path = joinpath(reaction_run_path, "baseline", "iso_massf00804.DAT"),
    iso_massf_filename = "iso_massf00804.DAT",
    verbose = true,
)
    base = df_base === nothing ? read_iso_massf(base_path) : df_base
    base = prepare_base_massf(base)
    factor_labels = factor_to_folder.(factors)
    factor_cols = Symbol.(factor_labels)

    df_final = DataFrame()

    for (reaction_name, isotopes_keep) in reaction_isotopes
        verbose && println("Processing: ", reaction_name)

        isotopes_keep = string.(isotopes_keep)
        base_subset = filter(row -> row.isotope in isotopes_keep, base)
        allowmissing!(base_subset, :X_base)

        missing_isotopes = setdiff(isotopes_keep, base_subset.isotope)
        if !isempty(missing_isotopes)
            append!(
                base_subset,
                DataFrame(
                    isotope = missing_isotopes,
                    X_base = Union{Missing, Float64}[missing for _ in missing_isotopes],
                ),
                cols = :union,
            )
        end

        df_long = DataFrame(
            reaction = String[],
            isotope = String[],
            factor = String[],
            ratio = Union{Missing, Float64}[],
        )

        for (factor, factor_label) in zip(factors, factor_labels)
            filepath, candidates = first_existing_factor_file(
                reaction_run_path,
                reaction_name,
                factor,
                iso_massf_filename,
                base_path,
            )

            if filepath === nothing
                verbose && println("Missing: ", first(candidates))
                append!(
                    df_long,
                    DataFrame(
                        reaction = fill(reaction_name, nrow(base_subset)),
                        isotope = base_subset.isotope,
                        factor = fill(factor_label, nrow(base_subset)),
                        ratio = Union{Missing, Float64}[missing for _ in 1:nrow(base_subset)],
                    ),
                )
                continue
            end

            df_run = rename(read_iso_massf(filepath), :X => :X_run)
            merged = leftjoin(base_subset, df_run, on = :isotope)
            ratios = [
                (ismissing(x_run) || ismissing(x_base)) ? missing : x_run / x_base
                for (x_run, x_base) in zip(merged.X_run, merged.X_base)
            ]

            append!(
                df_long,
                DataFrame(
                    reaction = fill(reaction_name, nrow(merged)),
                    isotope = merged.isotope,
                    factor = fill(factor_label, nrow(merged)),
                    ratio = ratios,
                ),
            )
        end

        df_table = unstack(df_long, [:reaction, :isotope], :factor, :ratio)
        df_table = df_table[:, [:reaction, :isotope, factor_cols...]]
        append!(df_final, df_table, cols = :union)
    end

    return df_final
end
