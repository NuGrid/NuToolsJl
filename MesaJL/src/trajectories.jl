using Printf

"""
    extract_zone_trajectory(logs_dir; target_mass = nothing) -> DataFrame

Extract the (time, temperature, density) history of a single Lagrangian
mass zone across every saved profile in a MESA `LOGS/` directory, for
feeding a single-zone nucleosynthesis post-processor (nuppn's `ppn`).

Every `profileN.data` listed in `profiles.index` is read via
[`read_profile`](@ref); within each one, the zone whose `mass` column (in
Msun) is closest to `target_mass` is picked, and its `temperature` and
`logRho` are recorded against that profile's `star_age`. MESA remeshes
between profiles, so zone *index* isn't a stable tracer — matching on the
`mass` coordinate itself is.

If `target_mass` is left as `nothing`, it is chosen automatically as the
hottest-reaching zone *among those that actually end up ejected* (mass
coordinate above the final profile's surface mass, i.e. zones that are no
longer part of the star once the wind is done) — not simply the single
hottest zone over the whole run. The unconditional hottest zone is almost
always just inside the retained envelope (deepest layers are hottest, and
the deepest layers are the last to be stripped, if they're stripped at
all), so picking it would track material that settles back into the WD
rather than material that is ejected — not useful for ejecta
nucleosynthesis, which is the point of this extraction. If nothing was
ejected over the run at all, falls back to the unconditional hottest zone
with a warning printed.

Once `target_mass` exceeds the current profile's surface mass coordinate
(i.e. the wind has stripped away everything out to and including that
zone), the tracked parcel is no longer part of MESA's Lagrangian grid at
all — it has actually been ejected. At that point extraction **stops**:
later profiles are not included. The alternative (nearest-available-zone
matching) would silently glue the real zone's pre-ejection history onto the
*new, post-ejection surface's* history instead, producing a spurious
kink/discontinuity in the trajectory that has nothing to do with the
tracked parcel's actual thermodynamics.

Returns a `DataFrame` with columns `model_number`, `age_yr`, `mass_Msun`
(the zone's actual matched mass coordinate, for sanity-checking against
`target_mass`), `temperature_K`, `density_cgs`, sorted by increasing model
number (so increasing age). Also prints a note if the trajectory was cut
short by ejection.
"""
function extract_zone_trajectory(logs_dir; target_mass = nothing)
    idx = read_profiles_index(joinpath(logs_dir, "profiles.index"))
    order = sortperm(idx.model_number)
    profile_paths = [joinpath(logs_dir, "profile$(p).data")
                      for p in idx.profile_number[order]]

    profiles = [read_profile(p) for p in profile_paths]

    if target_mass === nothing
        final_surface_mass = maximum(profiles[end].zones.mass)
        peak_T = -Inf
        for pr in profiles
            ejected = pr.zones.mass .> final_surface_mass
            any(ejected) || continue
            local_zones = findall(ejected)
            k = local_zones[argmax(pr.zones.temperature[local_zones])]
            if pr.zones.temperature[k] > peak_T
                peak_T = pr.zones.temperature[k]
                target_mass = pr.zones.mass[k]
            end
        end
        if target_mass === nothing
            println("extract_zone_trajectory: nothing was ejected over this run; ",
                    "falling back to the unconditional hottest zone (will track ",
                    "retained, not ejected, material).")
            peak_T = -Inf
            for pr in profiles
                k = argmax(pr.zones.temperature)
                if pr.zones.temperature[k] > peak_T
                    peak_T = pr.zones.temperature[k]
                    target_mass = pr.zones.mass[k]
                end
            end
        end
    end

    model_number = Int[]
    age_yr = Float64[]
    mass_Msun = Float64[]
    temperature_K = Float64[]
    density_cgs = Float64[]

    ever_present = false
    for pr in profiles
        zones = pr.zones
        if target_mass > maximum(zones.mass)
            if ever_present
                # was in the grid before, isn't now: genuinely ejected, stop for good.
                println("extract_zone_trajectory: target_mass = ", target_mass,
                        " Msun left the grid (was ejected) before model ",
                        pr.info.model_number, " (age = ", pr.info.star_age,
                        " yr) — trajectory truncated there.")
                break
            else
                # not accreted yet at this early profile: skip, don't stop.
                continue
            end
        end
        ever_present = true
        k = argmin(abs.(zones.mass .- target_mass))
        push!(model_number, pr.info.model_number)
        push!(age_yr, pr.info.star_age)
        push!(mass_Msun, zones.mass[k])
        push!(temperature_K, zones.temperature[k])
        push!(density_cgs, 10.0^zones.logRho[k])
    end

    return DataFrame(model_number = model_number, age_yr = age_yr,
                      mass_Msun = mass_Msun, temperature_K = temperature_K,
                      density_cgs = density_cgs)
end

"""
    write_ppn_trajectory(traj, path; time_col = :age_yr,
                          temp_col = :temperature_K, dens_col = :density_cgs,
                          id = "0123456789")

Write `traj` (as returned by [`extract_zone_trajectory`](@ref)) to `path` in
the single-zone `ppn` `trajectory.input` format that nuppn's `ppn` frame
reads (`nuppn/utils/source/trajectories.F90`): a fixed 7-line header — 3
comment lines, then `AGEUNIT`/`TUNIT`/`RHOUNIT` tags read with the Fortran
format `(10x,A3)`, then an `ID` line — followed by one `time  T9  rho[cgs]`
row per timestep. The header spacing is load-bearing (the unit code must
land in columns 11-13); don't reflow the `println` lines below. Ages are
written in years, temperatures in GK (T9 = T/1e9 K), densities in g/cc,
matching the tags.
"""
function write_ppn_trajectory(traj, path; time_col = :age_yr,
                               temp_col = :temperature_K,
                               dens_col = :density_cgs,
                               id = "0123456789")
    open(path, "w") do io
        println(io, "# time         T       rho")
        println(io, "# YRS/SEC; T8K/T9K; CGS/LOG")
        println(io, "# FORMAT: '(10x,A3)'")
        println(io, "AGEUNIT = YRS")
        println(io, "TUNIT   = T9K")
        println(io, "RHOUNIT = CGS")
        println(io, "ID = ", id)
        for row in eachrow(traj)
            t9 = row[temp_col] / 1.0e9
            @printf(io, " %.5e  %.4f  %.5e\n", row[time_col], t9, row[dens_col])
        end
    end
end
