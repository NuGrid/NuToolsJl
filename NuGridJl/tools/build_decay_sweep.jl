# build_decay_sweep.jl — post-process finished trajectory run(s) through
# nuppn's built-in decay-only mode (physics_knobs.F90's `decay`/`decay_time`),
# once unstable species (Na-22, Al-26, C-14, Be-7, ...) are allowed to finish
# decaying instead of being read off mid-decay at the trajectory's own last
# cycle. Two orchestration layers over the same per-run building block
# ([`build_decay_run!`](@ref)): [`build_decay_sweep`](@ref) (one source run,
# many decay times -- e.g. to find which decay time best matches a reference
# baseline) and [`decay_run_tree`](@ref) (one decay time, many source runs --
# e.g. to mirror an entire `build_sweep.jl` sensitivity sweep post-decay).
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
include(joinpath(@__DIR__, "namelist_utils.jl"))

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
    build_decay_run!(source_dir, dest_dir, decay_time)

The shared, one-run building block behind both [`build_decay_sweep`](@ref)
(one source run, many decay times) and [`decay_run_tree`](@ref) (one decay
time, many source runs): copy `source_dir` to `dest_dir` (its compiled
`ppn.exe` is symlinked, not copied -- decay mode is a runtime knob, no
recompilation needed), seed `post_abundance.DAT` from `source_dir`'s last
trajectory cycle ([`write_post_abundance`](@ref)), and patch
`ppn_physics.input`/`ppn_frame.input` for decay mode
([`patch_decay_inputs!`](@ref)).
"""
function build_decay_run!(source_dir::AbstractString, dest_dir::AbstractString, decay_time::Real)
    _copy_template!(source_dir, dest_dir)
    write_post_abundance(last_cycle_file(source_dir), joinpath(dest_dir, "post_abundance.DAT"))
    patch_decay_inputs!(dest_dir, decay_time)
end

"""
    build_decay_sweep(baseline_dir, out_dir, decay_times;
                       jobs = 4, dry_run = false) -> Vector{Pair{String,Float64}}

For each `(label, seconds)` pair in `decay_times`, build `out_dir/decay_<label>/`
from `baseline_dir` via [`build_decay_run!`](@ref). Unless `dry_run`, runs
all of them via [`run_parallel`](@ref). Read a result back with
`abundances(PPNRun(dir), :decay)`. Returns the `(directory, seconds)` pairs
built.
"""
function build_decay_sweep(baseline_dir::AbstractString, out_dir::AbstractString,
                            decay_times::Vector{<:Pair}; jobs::Integer = 4, dry_run::Bool = false)
    built = Pair{String,Float64}[]
    for (label, seconds) in decay_times
        dest = joinpath(out_dir, "decay_$(label)")
        build_decay_run!(baseline_dir, dest, seconds)
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

"""
    discover_ppn_runs(tree_dir) -> Vector{String}

Every directory under `tree_dir` (`tree_dir` included) that's a genuine
trajectory-integrated source run -- has `ppn.exe` *and* at least one
numbered `iso_massf#####.DAT` cycle -- e.g. a `build_sweep.jl` output's
`baseline/` plus every `<reaction>/fact_<factor>/`. The numbered-cycle
requirement is what excludes a previous [`build_decay_sweep`](@ref)/
[`decay_run_tree`](@ref) output sitting nested in the same tree (e.g.
`decay_analysis/decay_3hr/`): decay-mode runs only ever write the single
`iso_massfdecay.DAT`, never a numbered one, so re-decaying a tree that
already contains decayed output doesn't pick those back up as if they were
new sources. Sorted for deterministic ordering.
"""
function discover_ppn_runs(tree_dir::AbstractString)
    runs = String[]
    for (root, _, files) in walkdir(tree_dir)
        "ppn.exe" in files && any(f -> occursin(r"^iso_massf\d+\.DAT$", f), files) && push!(runs, root)
    end
    return sort!(runs)
end

"""
    decay_run_tree(source_tree_dir, out_dir, decay_time;
                    jobs = 4, dry_run = false) -> Vector{String}

Mirror an entire `build_sweep.jl` output (`source_tree_dir`, discovered via
[`discover_ppn_runs`](@ref) -- `baseline/` plus every `<reaction>/
fact_<factor>/`) into `out_dir`, decaying every single run by the same
`decay_time` (seconds) via [`build_decay_run!`](@ref), preserving the
relative directory structure so a sensitivity comparison can be redone
post-decay the same way it was done pre-decay. Complementary to
[`build_decay_sweep`](@ref) (one run, many times) -- this is one time, many
runs. Unless `dry_run`, runs all of them via [`run_parallel`](@ref). Returns
the built directories (mirroring `discover_ppn_runs(source_tree_dir)`
1:1, same relative paths under `out_dir`).
"""
function decay_run_tree(source_tree_dir::AbstractString, out_dir::AbstractString,
                         decay_time::Real; jobs::Integer = 4, dry_run::Bool = false)
    sources = discover_ppn_runs(source_tree_dir)
    isempty(sources) && throw(ArgumentError("no ppn.exe found anywhere under $source_tree_dir"))

    built = String[]
    for source in sources
        dest = joinpath(out_dir, relpath(source, source_tree_dir))
        build_decay_run!(source, dest, decay_time)
        push!(built, dest)
    end

    if !dry_run
        results = run_parallel(built; jobs)
        for (dir, ok) in results
            ok || @warn "decay run failed" dir
        end
    end
    return built
end

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

"""
    auto_label(seconds) -> String

A human-readable directory-name label for a decay time in seconds (`60.0` ->
`"1min"`, `10800.0` -> `"3hr"`) -- what `--times`/`--decay-time` use to name
`decay_<label>/` when the CLI derives labels from bare numbers instead of
the `"label" => seconds` pairs the Julia API takes directly.
"""
function auto_label(seconds::Real)
    seconds == 0 && return "0s"
    minute, hour, day, year = 60.0, 3600.0, 86400.0, 31557600.0
    seconds < minute && return @sprintf("%gs", seconds)
    seconds < hour && return @sprintf("%gmin", seconds / minute)
    seconds < day && return @sprintf("%ghr", seconds / hour)
    seconds < year && return @sprintf("%gd", seconds / day)
    return @sprintf("%gyr", seconds / year)
end

"""
    parse_times(text) -> Vector{Pair{String,Float64}}

Parse a `--times` value: comma-separated entries, each either a bare number
of seconds (auto-labeled via [`auto_label`](@ref)) or an explicit
`label=seconds` pair (e.g. `"1min=60,custom=12345"`).
"""
function parse_times(text::AbstractString)
    pairs = Pair{String,Float64}[]
    for chunk in split(text, ',')
        chunk = strip(chunk)
        isempty(chunk) && continue
        if occursin('=', chunk)
            label, secs = split(chunk, '='; limit = 2)
            push!(pairs, String(strip(label)) => parse(Float64, strip(secs)))
        else
            secs = parse(Float64, chunk)
            push!(pairs, auto_label(secs) => secs)
        end
    end
    isempty(pairs) && throw(ArgumentError("--times must contain at least one value"))
    return pairs
end

function usage()
    println("""
Usage:
  julia build_decay_sweep.jl times <source_dir> <out_dir> --times s1,s2,... [options]
  julia build_decay_sweep.jl tree <source_tree_dir> <out_dir> --decay-time SECONDS [options]

Post-process finished nuppn trajectory run(s) through decay-only mode
(physics_knobs.F90's decay/decay_time), once unstable species (Na-22, Al-26,
C-14, Be-7, ...) are allowed to keep decaying past the trajectory's last cycle.

Modes:
  times   One source run, many decay times -- e.g. to find which decay time
          best matches a reference baseline.
  tree    One decay time, many source runs -- mirrors an entire
          build_sweep.jl output (baseline/ + every <reaction>/fact_<factor>/),
          decaying each run by the same time, so a sensitivity comparison can
          be redone post-decay.

Arguments (times mode):
  source_dir       A single finished trajectory run (has iso_massf#####.DAT files)
  out_dir          Where to build decay_<label>/ for each time

Arguments (tree mode):
  source_tree_dir  A build_sweep.jl output (or any tree of trajectory runs)
  out_dir          Mirrors source_tree_dir's structure into out_dir, decayed

Options:
  --times s1,s2,...       (times mode, required) comma-separated decay times
                          in seconds; each entry is either a bare number
                          (auto-labeled, e.g. "3600" -> decay_1hr/) or an
                          explicit "label=seconds" pair (e.g. "custom=12345")
  --decay-time SECONDS     (tree mode, required) single decay time in seconds
  --jobs N, -j N           Number of ppn.exe runs in parallel (default: 4)
  --dry-run                Build directories without launching ppn.exe
  -h, --help                Show this help

Examples:
  julia build_decay_sweep.jl times run/baseline analysis \\
      --times 60,300,600,3600,1yr=3.156e7 --jobs 8
  julia build_decay_sweep.jl tree run_sweep decay_run_sweep \\
      --decay-time 10800 --jobs 8
""")
end

if abspath(PROGRAM_FILE) == @__FILE__
    let jobs = 4, dry_run = false, times_arg = nothing, decay_time_arg = nothing,
        positional = String[], i = 1

        if isempty(ARGS) || ARGS[1] in ("-h", "--help")
            usage(); exit(0)
        end
        mode = ARGS[1]
        mode in ("times", "tree") || (println(stderr, "First argument must be 'times' or 'tree', got '$mode'.\n"); usage(); exit(1))
        i = 2
        while i <= length(ARGS)
            a = ARGS[i]
            if a in ("-h", "--help")
                usage(); exit(0)
            elseif a == "--times"
                i == length(ARGS) && error("--times requires a value")
                times_arg = ARGS[i + 1]; i += 2
            elseif a == "--decay-time"
                i == length(ARGS) && error("--decay-time requires a value")
                decay_time_arg = parse(Float64, ARGS[i + 1]); i += 2
            elseif a in ("--jobs", "-j")
                i == length(ARGS) && error("--jobs requires a value")
                jobs = parse(Int, ARGS[i + 1]); i += 2
            elseif a == "--dry-run"
                dry_run = true; i += 1
            elseif startswith(a, "-")
                println(stderr, "Unknown option: $a\n")
                usage(); exit(1)
            else
                push!(positional, a); i += 1
            end
        end
        if length(positional) != 2
            println(stderr, "Expected 2 positional arguments (<source> <out_dir>), got $(length(positional)).\n")
            usage(); exit(1)
        end
        source, out_dir = positional

        if mode == "times"
            times_arg === nothing && (println(stderr, "times mode requires --times.\n"); usage(); exit(1))
            built = build_decay_sweep(source, out_dir, parse_times(times_arg); jobs, dry_run)
            println(length(built), " decay run(s) at ", source, " -> ", out_dir)
        else
            decay_time_arg === nothing && (println(stderr, "tree mode requires --decay-time.\n"); usage(); exit(1))
            built = decay_run_tree(source, out_dir, decay_time_arg; jobs, dry_run)
            println(length(built), " decayed run(s) mirrored from ", source, " -> ", out_dir)
        end
    end
end
