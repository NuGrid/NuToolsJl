module MesaJL

using CSV
using DataFrames

export read_history, read_profile, read_profiles_index, profile_path_for_model

"""
    read_history(path) -> DataFrame

Read a MESA `LOGS/history.data` file. MESA's history/profile files share a
fixed preamble: a row of column indices, a row of names, one row of values
(run metadata, ignored here), a blank separator, then the same
indices/names/data pattern again for the per-model table that this function
returns. Column dtypes are inferred by CSV.jl; MESA prints logical columns
as `T`/`F`, which CSV.jl reads as `Bool`.
"""
function read_history(path)
    return CSV.read(path, DataFrame; delim = ' ', ignorerepeated = true,
                     header = 6, skipto = 7, truestrings = ["T"], falsestrings = ["F"])
end

"""
    read_profile(path) -> (info::DataFrameRow, zones::DataFrame)

Read a MESA `LOGS/profileN.data` file. Returns the single-row run-info table
(model_number, star_age, Teff, ...) as `info`, and the zone-by-zone structure
(mass, logR, logT, logRho, ...) as `zones`, one row per zone. Zone 1 is the
surface; the last zone is the center (MESA's usual convention).
"""
function read_profile(path)
    info_df = CSV.read(path, DataFrame; delim = ' ', ignorerepeated = true,
                        header = 2, skipto = 3, limit = 1)
    zones = CSV.read(path, DataFrame; delim = ' ', ignorerepeated = true,
                      header = 6, skipto = 7, truestrings = ["T"], falsestrings = ["F"])
    return (info = info_df[1, :], zones = zones)
end

"""
    read_profiles_index(path) -> DataFrame

Read a MESA `LOGS/profiles.index` file, mapping each saved profile number to
its model number (columns `model_number`, `priority`, `profile_number`).
"""
function read_profiles_index(path)
    return CSV.read(path, DataFrame; delim = ' ', ignorerepeated = true,
                     header = ["model_number", "priority", "profile_number"],
                     skipto = 2)
end

"""
    profile_path_for_model(logs_dir, model_number)

Find the `profileN.data` file in `logs_dir` corresponding to a given MESA
`model_number`, using `profiles.index`. Returns `nothing` if no profile was
saved for that model.
"""
function profile_path_for_model(logs_dir, model_number)
    idx = read_profiles_index(joinpath(logs_dir, "profiles.index"))
    row = findfirst(==(model_number), idx.model_number)
    row === nothing && return nothing
    return joinpath(logs_dir, "profile$(idx.profile_number[row]).data")
end

end # module MesaJL
