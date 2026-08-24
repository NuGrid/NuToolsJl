# build_decay_sweep.jl — post-process a finished trajectory run through
# nuppn's built-in decay-only mode (physics_knobs.F90's `decay`/`decay_time`),
# at a range of decay durations, to see how far a Table-4 comparison closes up
# once short-lived species (Na-22, Al-26, C-14, Be-7, ...) are allowed to
# finish decaying instead of being read off mid-decay at the trajectory's own
# last cycle.
#
# Ported from the proven working pattern in
# NovaSensitivityStudy/single-zone/tools/decay_ppn_sweep.jl and
# decay_time_scan.jl (per user pointer) after an `iabuini = 5` version of this
# script (reading the raw copied iso_massf#####.DAT directly via
# abundances.F90's `load_ppn`) crashed instantly and unsymbolized on every
# decay_time value, with no diagnosable difference in the input files. The
# working pattern instead uses `iabuini = 11` (`load_urs_frischknecht_xin`),
# which expects a *reformatted* abundance file (fixed-width `frisch_fmt =
# '(4x,a5,9x,d16.10)'` in abundances.F90), not nuppn's own iso_massf output
# format -- `write_post_abundance` below builds exactly that.
#
# nuppn's decay mode itself (confirmed by reading physics_knobs.F90, ppn.F90,
# decays.F90) is a single, non-iterative step: ppn.F90 loads the initial
# abundance, and if `decay = .true.` it calls `do_decay` (one call to
# `integrate_network` spanning all of `decay_time` seconds at a hardcoded
# T9=0.01/rho=10 -- cold enough that every charged-particle reaction is
# negligible and only weak/alpha decay matter) then writes a single output
# file and stops. The output is *not* a numbered `iso_massf#####.DAT` --
# ppn.F90's printonecycle hardcodes the cycle tag to the literal string
# "decay" whenever `decay` is set, i.e. `iso_massfdecay.DAT`.
#
# Since `decay`/`decay_time` are runtime ppn_physics.input knobs (not
# compile-time array-sizing parameters), every decay-time variant can share
# one already-compiled ppn.exe -- no recompilation needed, same as
# build_sweep.jl's fact_0.5/fact_2.0 variants (here: symlinked, not copied,
# matching decay_ppn_sweep.jl -- no reason to duplicate a ~300MB binary
# per decay-time directory).
#
# NuGridJl's `PPNRun`/`abundances` already know how to read `iso_massfdecay.DAT`
# directly via `abundances(run, :decay)` (see run.jl) -- no renaming needed.

include(joinpath(@__DIR__, "run_parallel.jl"))

using Printf

const _STALE_OUTPUT_RE = r"^(iso_massf.*\.DAT|flux_.*\.DAT|x-time\.dat|OUT|fort\.6|run\.log)$"i

function _copy_template!(template_dir::AbstractString, dest_dir::AbstractString)
    mkpath(dest_dir)
    npdata_target = nothing
    for candidate in (joinpath(template_dir, "NPDATA"), joinpath(template_dir, "..", "NPDATA"))
        if islink(candidate) || isdir(candidate)
            npdata_target = realpath(candidate)
            break
        end
    end
    for name in readdir(template_dir)
        name in ("NPDATA", "ppn.exe") && continue
        occursin(_STALE_OUTPUT_RE, name) && continue
        src = joinpath(template_dir, name)
        dst = joinpath(dest_dir, name)
        ispath(dst) && rm(dst; recursive = true, force = true)
        if islink(src)
            symlink(readlink(src), dst)
        else
            cp(src, dst; force = true)
        end
    end
    # ppn.exe is symlinked, not copied -- every decay-time variant shares one binary
    exe_dst = joinpath(dest_dir, "ppn.exe")
    ispath(exe_dst) || symlink(relpath(realpath(joinpath(template_dir, "ppn.exe")), dest_dir), exe_dst)
    if npdata_target !== nothing
        # ppn_physics.input references data files via '../NPDATA/...' --
        # relative to the run directory's *parent* -- so both dest_dir/NPDATA
        # and dirname(dest_dir)/NPDATA need to exist (decay_<label>/ sits one
        # level deeper than baseline/ did, under decay_analysis/, so the
        # parent-level link lands somewhere build_sweep.jl's own two-level
        # case never needed to but the identical reasoning still applies).
        for npdata_link in (joinpath(dest_dir, "NPDATA"), joinpath(dirname(dest_dir), "NPDATA"))
            ispath(npdata_link) || symlink(npdata_target, npdata_link)
        end
    end
end

"""
    last_cycle_file(baseline_dir) -> String

The highest-numbered `iso_massf#####.DAT` in `baseline_dir` -- the trajectory
run's final abundance, and the seed for every decay-time variant.
"""
function last_cycle_file(baseline_dir::AbstractString)
    candidates = filter(f -> occursin(r"^iso_massf\d+\.DAT$", f), readdir(baseline_dir))
    isempty(candidates) && throw(ArgumentError("no iso_massf#####.DAT files found in $baseline_dir"))
    return joinpath(baseline_dir, sort(candidates)[end])
end

"""
    parse_iso_massf_rows(path) -> Vector{NamedTuple}

Parse a raw `iso_massf#####.DAT`: each data row is `index Z A isom abundance
name`; header/comment lines (non-numeric first token) are skipped.
"""
function parse_iso_massf_rows(path::AbstractString)
    isfile(path) || throw(ArgumentError("missing iso_massf input: $path"))
    rows = NamedTuple[]
    for line in eachline(path)
        parts = split(line)
        isempty(parts) && continue
        all(isdigit, parts[1]) || continue
        push!(rows, (
            z = round(Int, parse(Float64, parts[2])),
            a = round(Int, parse(Float64, parts[3])),
            abundance = parse(Float64, parts[5]),
            label_parts = parts[6:end],
        ))
    end
    isempty(rows) && throw(ArgumentError("no abundance rows found in $path"))
    return rows
end

function _format_iso_label(a::Integer, label_parts)
    if length(label_parts) >= 2
        return @sprintf("%-2s%3d", lowercase(label_parts[1]), a)
    end
    # single token: normally a heavy species whose 3-digit mass number fills
    # zis's "(A2,I3)" field with no internal space (e.g. "FE254") -- split
    # and reformat that case for consistency with the length>=2 branch above.
    # Anything else single-token (NEUT, PROT, OOOOO, or one of this project's
    # isomer-tagged names like "AL26M" -- ppn_physics.F90's "(A2,I2,A1)"
    # naming for a build-time isomer species, exactly 5 characters already)
    # is already a complete, correctly-formatted zis-style name -- pass it
    # through verbatim rather than erroring on a pattern it was never meant
    # to match.
    m = match(r"^([A-Za-z]+)(\d+)$", label_parts[1])
    m === nothing && return label_parts[1]
    return @sprintf("%-2s%3d", lowercase(m.captures[1]), a)
end

"""
    write_post_abundance(source_iso_massf, out_path)

Reformat a raw `iso_massf#####.DAT` into the fixed-width format
`abundances.F90`'s `iabuini = 11` path (`load_urs_frischknecht_xin`,
`frisch_fmt = '(4x,a5,9x,d16.10)'`) actually expects: not the same format
nuppn writes its own per-cycle output in.
"""
function write_post_abundance(source_iso_massf::AbstractString, out_path::AbstractString)
    rows = parse_iso_massf_rows(source_iso_massf)
    open(out_path, "w") do io
        for row in rows
            iso = _format_iso_label(row.a, row.label_parts)
            @printf(io, "%3d %-5s         %16.10E\n", row.z, iso, row.abundance)
        end
    end
end

"""
    update_namelist(text, replacements) -> String

Set each `name => value` pair in a Fortran namelist's text: overwrite the
value in place if `name` already has a line, otherwise insert it just before
the terminating `/`. Case-insensitive on `name` (Fortran namelists are).
"""
function update_namelist(text::AbstractString, replacements::Vector{<:Pair})
    lines = split(text, '\n')
    !isempty(lines) && lines[end] == "" && pop!(lines)

    replacement_keys = Dict(lowercase(k) => v for (k, v) in replacements)
    found = Set{String}()
    output = String[]

    for line in lines
        if strip(line) == "/"
            for (key, value) in replacements
                lowercase(key) in found || push!(output, "        $key = $value")
            end
            push!(output, line)
            continue
        end
        m = match(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
        if m !== nothing && haskey(replacement_keys, lowercase(m.captures[2]))
            key = m.captures[2]
            push!(found, lowercase(key))
            push!(output, "$(m.captures[1])$key = $(replacement_keys[lowercase(key)])")
        else
            push!(output, line)
        end
    end
    return join(output, '\n') * "\n"
end

fortran_double(value::Real) = replace(@sprintf("%.10E", value), "E" => "d")

"""
    patch_decay_inputs!(decay_run, decay_time)

Edit `decay_run`'s already-copied `ppn_frame.input`/`ppn_physics.input` in
place for decay mode: `iabuini = 11` reading `post_abundance.DAT`
(`write_post_abundance` must have already written it), `decay = .true.` /
`decay_time`, `detailed_balance = .false.` (physics_knobs.F90 forces this
whenever `decay` is set anyway, printing a warning if it wasn't already
false -- setting it here just avoids that warning).
"""
function patch_decay_inputs!(decay_run::AbstractString, decay_time::Real)
    frame = joinpath(decay_run, "ppn_frame.input")
    physics = joinpath(decay_run, "ppn_physics.input")
    write(frame, update_namelist(read(frame, String), [
        "nsource" => "0",
        "iabuini" => "11",
        "ini_filename" => "'post_abundance.DAT'",
        "iplot_flux_option" => "0",
        "i_flux_integrated" => "0",
    ]))
    write(physics, update_namelist(read(physics, String), [
        "decay" => ".true.",
        "decay_time" => fortran_double(decay_time),
        "detailed_balance" => ".false.",
    ]))
end

"""
    build_decay_sweep(baseline_dir, out_dir, decay_times;
                       jobs = 4, dry_run = false) -> Vector{Pair{String,Float64}}

For each `(label, seconds)` pair in `decay_times`, build `out_dir/decay_<label>/`
by copying `baseline_dir` (its compiled `ppn.exe` is symlinked, not copied --
decay mode is a runtime knob, no recompilation needed), seeding
`post_abundance.DAT` from `baseline_dir`'s last trajectory cycle
([`write_post_abundance`](@ref)), and patching `ppn_physics.input`/
`ppn_frame.input` for decay mode ([`patch_decay_inputs!`](@ref)). Unless
`dry_run`, runs all of them via [`run_parallel`](@ref). Read a result back
with `abundances(PPNRun(dir), :decay)`. Returns the `(directory, seconds)`
pairs built.
"""
function build_decay_sweep(baseline_dir::AbstractString, out_dir::AbstractString,
                            decay_times::Vector{<:Pair}; jobs::Integer = 4, dry_run::Bool = false)
    seed_file = last_cycle_file(baseline_dir)

    built = Pair{String,Float64}[]
    for (label, seconds) in decay_times
        dest = joinpath(out_dir, "decay_$(label)")
        _copy_template!(baseline_dir, dest)
        write_post_abundance(seed_file, joinpath(dest, "post_abundance.DAT"))
        patch_decay_inputs!(dest, seconds)
        built = push!(built, dest => Float64(seconds))
    end

    if !dry_run
        dirs = first.(built)
        results = run_parallel(dirs; jobs)
        for (dir, ok) in results
            ok || @warn "decay run failed" dir
        end
    end
    return built
end
