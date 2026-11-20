#!/usr/bin/env julia
# run_nova_grid.jl — build a grid of MESA nova cases from a template
# directory, varying accretion/envelope/output knobs, and run them all in
# parallel via MesaJL.run_mesa_parallel.
#
# Usage:
#   julia run_nova_grid.jl [mesa_root] [jobs]
#
# mesa_root (default: PostTNR/mesa/) is where new case directories get
# created. jobs (default: 4) is how many cases run concurrently.
#
# Edit the `cases` vector near the bottom of this file to define your grid,
# then run it. Each NovaCase only has to set the knobs it wants to change
# from the template; everything else is left as the template already has it.
#
# IMPORTANT LIMITATIONS, read before reaching for a knob that isn't here:
#   - WD mass is NOT a free numeric knob. It comes from whichever starting
#     .mod you point `start_model` at, and only `1.1M_lgTc_7.7.mod` exists
#     locally today. To get another mass, run $MESA_DIR/star/test_suite/
#     make_co_wd (or make_o_ne_wd for ONe) for the mass/Tc you want, then
#     point `start_model` at the .mod it produces -- this script builds
#     nothing new on its own.
#   - nova_type = "one" (ONe WD) is NOT wired up: there's no ONe starting
#     model or ONe-appropriate network set up locally yet (only "co" is).
#     Requesting it is a hard error, not a silent fallback.
#   - Only inlist_wd_nova_burst is touched by the named knobs below (that's
#     where accretion/composition/network/output-cadence actually live in
#     this test case's inlist split -- inlist_setup only relaxes the
#     starting model and sets the network/v_flag once). Use
#     `extra_controls` for anything else in that same file; if you need to
#     edit inlist_setup itself, do it by hand on the built case directory
#     before running.

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using MesaJL

const POSTTNR = joinpath(@__DIR__, "..", "..", "..")
const DEFAULT_TEMPLATE = joinpath(POSTTNR, "mesa", "co_nova_1.1_fiducial")

"""
    NovaCase

One grid point for [`build_case`](@ref). `label` names the new directory
(`<mesa_root>/<label>/`). `template` is the case directory copied as the
starting point for it (must already be a working, compiled MESA case: has
`make/`, a starting `.mod`, etc. — see `co_nova_1.1_fiducial`). Every other
field is a knob into `inlist_wd_nova_burst`; leave a field at its default
(`nothing` / empty) to keep whatever the template already has.
"""
Base.@kwdef struct NovaCase
    label::String
    nova_type::String = "co"                        # only "co" is wired up locally — see limitations above
    template::String = DEFAULT_TEMPLATE
    start_model::Union{Nothing,String} = nothing     # path to a .mod (sets WD mass/Tc); nothing = keep template's
    accretion_rate::Union{Nothing,Float64} = nothing # Msun/yr -> mass_change
    accretion_h1::Union{Nothing,Float64} = nothing
    accretion_he3::Union{Nothing,Float64} = nothing
    accretion_he4::Union{Nothing,Float64} = nothing
    net_name::Union{Nothing,String} = nothing        # reaction network, e.g. "cno_extras_o18_to_mg26_plus_fe56.net"
    profile_interval::Union{Nothing,Int} = nothing
    extra_controls::Dict{String,String} = Dict{String,String}()  # raw &controls key => already-Fortran-formatted value
end

"""
    set_inlist_value!(path, key, value_str)

Replace the value of an existing `key = ...` line in a MESA inlist file
with `value_str` (a literal Fortran value as a string — quote it yourself
for strings, e.g. `"'some.net'"`). Errors if `key` isn't already present:
this only edits existing namelist entries, it never appends new ones, so a
typo'd key fails loudly instead of silently landing outside any namelist.
"""
function set_inlist_value!(path::AbstractString, key::AbstractString, value_str::AbstractString)
    lines = readlines(path)
    pattern = Regex("^(\\s*" * key * "\\s*=\\s*).*\$")
    for i in eachindex(lines)
        m = match(pattern, lines[i])
        if m !== nothing
            lines[i] = m.captures[1] * value_str
            write(path, join(lines, "\n") * "\n")
            return nothing
        end
    end
    error("set_inlist_value!: key '$key' not found in $path")
end

"""
    build_case(case::NovaCase, mesa_root) -> String

Copy `case.template` to `<mesa_root>/<case.label>/`, strip its prior run
output (`LOGS/`, `photos/`, `.mesa_temp_cache/`, saved `.mod`s, `run.log`)
so it starts clean, then apply `case`'s knobs to `inlist_wd_nova_burst` (and
`inlist_setup` for `start_model`/`net_name`). Returns the new directory.
Errors if the target directory already exists, or if `case.nova_type`
isn't `"co"` (see the module-level limitations note).
"""
function build_case(case::NovaCase, mesa_root::AbstractString)
    case.nova_type == "co" ||
        error("build_case: nova_type = \"$(case.nova_type)\" isn't wired up locally " *
              "(no starting model/network set up for it) — see the limitations note " *
              "at the top of run_nova_grid.jl.")

    dir = joinpath(mesa_root, case.label)
    isdir(dir) && error("build_case: $dir already exists; remove it or pick a different label")
    cp(case.template, dir)

    for sub in ("LOGS", "photos", ".mesa_temp_cache")
        rm(joinpath(dir, sub); recursive = true, force = true)
    end
    for f in ("final.mod", "ready.mod", "run.log")
        rm(joinpath(dir, f); force = true)
    end

    burst_inlist = joinpath(dir, "inlist_wd_nova_burst")
    setup_inlist = joinpath(dir, "inlist_setup")

    case.start_model === nothing || set_inlist_value!(setup_inlist, "load_model_filename", "'$(case.start_model)'")
    case.net_name === nothing || set_inlist_value!(setup_inlist, "new_net_name", "'$(case.net_name)'")

    case.accretion_rate === nothing || set_inlist_value!(burst_inlist, "mass_change", string(case.accretion_rate))
    case.accretion_h1 === nothing || set_inlist_value!(burst_inlist, "accretion_h1", string(case.accretion_h1))
    case.accretion_he3 === nothing || set_inlist_value!(burst_inlist, "accretion_he3", string(case.accretion_he3))
    case.accretion_he4 === nothing || set_inlist_value!(burst_inlist, "accretion_he4", string(case.accretion_he4))
    case.profile_interval === nothing || set_inlist_value!(burst_inlist, "profile_interval", string(case.profile_interval))

    for (key, value_str) in case.extra_controls
        set_inlist_value!(burst_inlist, key, value_str)
    end

    return dir
end

# -----------------------------------------------------------------------
# Define your grid here, then run this file.
# -----------------------------------------------------------------------
cases = [
    NovaCase(label = "co_nova_1.1_fiducial_A", accretion_rate = 1e-10),
    NovaCase(label = "co_nova_1.1_fiducial_C", accretion_rate = 1e-9),
]

mesa_root = length(ARGS) >= 1 ? ARGS[1] : joinpath(POSTTNR, "mesa")
jobs = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 4

dirs = [build_case(case, mesa_root) for case in cases]
println("Built $(length(dirs)) case(s): ", join(basename.(dirs), ", "))

results = run_mesa_parallel(dirs; jobs = jobs)
for (dir, ok) in results
    println(basename(dir), ": ", ok ? "ok" : "FAILED (see run.log)")
end
