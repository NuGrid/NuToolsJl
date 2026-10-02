#!/usr/bin/env julia

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools
using Dates

Base.@kwdef mutable struct SweepOptions
    nova::String = "nova_test"
    jobs::Int = get(ENV, "MAX_JOBS", nothing) === nothing ? 4 : parse(Int, ENV["MAX_JOBS"])
    baseline_only::Bool = false
    build_only::Bool = false
    apply_network_edits::Bool = false
    starlib_option::Int = 1
    runs_name::String = "runs_star"
end

function usage()
    println("""
Usage:
  julia tools/run_ppn_sweep_starlib.jl [options]

Options:
  --nova NAME             Nova directory under nova_cases/ (default: nova_test)
  --jobs N               Maximum concurrent ppn jobs (default: MAX_JOBS or 4)
  --baseline-only        Build and run only the baseline; do not build factored runs
  --build-only           Build all run directories but do not execute ppn.exe for factored runs
  --apply-network-edits  Apply config/network_edits.json before building
  --starlib-option N     STARLIB option to write into ppn_physics.input:
                         1 = MC10+MC13 (STARLIB 2010), 2 = ETR25 2025 (default: 1)
  --runs-name NAME       Output runs directory under the nova case (default: runs_star)
  -h, --help             Show this help

Note:
  The baseline is always run first with the requested starlib_option so that its
  networksetup.txt (which reflects the actual PPN reaction ordering for that STARLIB
  option) can be used to resolve reaction indices for the factored runs. This is
  essential for options that insert extra reactions (e.g. STARLIB v6.10 adds ~42
  reactions that shift indices relative to old REACLIB).
""")
end

function parse_args(args)
    opts = SweepOptions()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--nova"
            i += 1
            i <= length(args) || error("--nova requires a value")
            opts.nova = args[i]
        elseif arg == "--jobs"
            i += 1
            i <= length(args) || error("--jobs requires a value")
            opts.jobs = parse(Int, args[i])
            opts.jobs >= 1 || error("--jobs must be >= 1")
        elseif arg == "--baseline-only"
            opts.baseline_only = true
        elseif arg == "--build-only"
            opts.build_only = true
        elseif arg == "--apply-network-edits"
            opts.apply_network_edits = true
        elseif arg == "--starlib-option"
            i += 1
            i <= length(args) || error("--starlib-option requires a value")
            opts.starlib_option = parse(Int, args[i])
            opts.starlib_option in (1, 2) || error("--starlib-option must be 1 or 2")
        elseif arg == "--runs-name"
            i += 1
            i <= length(args) || error("--runs-name requires a value")
            opts.runs_name = args[i]
            isempty(strip(opts.runs_name)) && error("--runs-name cannot be empty")
        elseif arg in ("-h", "--help")
            usage()
            exit(0)
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end
    return opts
end

function relative_run_name(runs_dir, run_dir)
    rel = relpath(run_dir, runs_dir)
    replace(rel, '/' => '_')
end

function list_ppn_executables(runs_dir; baseline_only=false)
    if baseline_only
        exe = joinpath(runs_dir, "baseline", "ppn.exe")
        isfile(exe) || error("Missing baseline executable: $exe")
        return [exe]
    end

    exes = String[]
    for (root, _, files) in walkdir(runs_dir)
        if "ppn.exe" in files
            push!(exes, joinpath(root, "ppn.exe"))
        end
    end
    sort!(exes)
    isempty(exes) && error("No ppn.exe files found under $runs_dir")
    return exes
end

function set_starlib_option!(path, option)
    lines = readlines(path, keep=true)
    changed = false
    for i in eachindex(lines)
        if occursin(r"(?i)^\s*starlib_option\s*=", lines[i])
            lines[i] = replace(lines[i], r"(?i)^(\s*starlib_option\s*=\s*)\d+" => SubstitutionString("\\g<1>$option"); count=1)
            changed = true
            break
        end
    end
    if !changed
        insert_at = findfirst(line -> strip(line) == "/", lines)
        insert_at === nothing && error("Could not find ppn_physics namelist terminator in $path")
        insert!(lines, insert_at, "        starlib_option = $option\n")
    end
    write(path, join(lines))
end

function enable_starlib_in_runs!(runs_dir, option)
    count = 0
    for (root, _, files) in walkdir(runs_dir)
        if "ppn_physics.input" in files
            set_starlib_option!(joinpath(root, "ppn_physics.input"), option)
            count += 1
        end
    end
    count > 0 || error("No ppn_physics.input files found under $runs_dir")
    println("Enabled starlib_option=$option in $count ppn_physics.input file(s).")
end

function run_and_log(cmd::Cmd, logfile)
    open(logfile, "a") do io
        process = run(pipeline(cmd, stdout=io, stderr=io), wait=false)
        wait(process)
        success(process) || error("Command failed: $cmd")
    end
end

function run_one_ppn(exe, runs_dir, logs_dir)
    run_dir = dirname(exe)
    name = relative_run_name(runs_dir, run_dir)
    logfile = joinpath(logs_dir, "$name.log")
    start_time = time()

    open(logfile, "w") do io
        println(io, "===========================================")
        println(io, "Running $run_dir")
        println(io, "Started $(Dates.now())")
        println(io, "===========================================")
    end

    println("Running $name")
    run_and_log(Cmd(`./ppn.exe`, dir=run_dir), logfile)

    elapsed = round(Int, time() - start_time)
    open(logfile, "a") do io
        println(io, "Finished $run_dir in $elapsed seconds")
    end
    println("Finished $name in $(elapsed)s")
end

function run_parallel(exes, runs_dir, logs_dir, jobs)
    queue = Channel{String}(length(exes))
    for exe in exes
        put!(queue, exe)
    end
    close(queue)

    failures = Channel{Tuple{String,Any}}(length(exes))
    workers = Task[]

    for _ in 1:min(jobs, length(exes))
        push!(workers, @async begin
            for exe in queue
                try
                    run_one_ppn(exe, runs_dir, logs_dir)
                catch err
                    put!(failures, (exe, err))
                end
            end
        end)
    end

    foreach(wait, workers)
    close(failures)

    failed = collect(failures)
    if !isempty(failed)
        println("Failures:")
        for (exe, err) in failed
            println("  $exe")
            println("    $err")
        end
        error("$(length(failed)) ppn jobs failed")
    end
end

function main()
    opts = parse_args(ARGS)

    ENV["OMP_NUM_THREADS"]        = get(ENV, "OMP_NUM_THREADS", "1")
    ENV["OPENBLAS_NUM_THREADS"]   = get(ENV, "OPENBLAS_NUM_THREADS", "1")
    ENV["MKL_NUM_THREADS"]        = get(ENV, "MKL_NUM_THREADS", "1")
    ENV["BLIS_NUM_THREADS"]       = get(ENV, "BLIS_NUM_THREADS", "1")
    ENV["VECLIB_MAXIMUM_THREADS"] = get(ENV, "VECLIB_MAXIMUM_THREADS", "1")

    base = NovaRunTools.nova_dir(opts.nova)
    logs_dir = joinpath(base, "logs_$(opts.runs_name)")
    mkpath(logs_dir)
    runs_dir = joinpath(base, opts.runs_name)

    global_start = time()

    if opts.apply_network_edits
        println("Applying network edits before building runs...")
        setup_network(opts.nova; check=false)
    else
        println("Skipping network setup; using current networksetup.txt unchanged.")
    end

    if ispath(runs_dir)
        println("Deleting existing STARLIB runs directory: $runs_dir")
        rm(runs_dir; recursive=true)
    else
        println("No existing STARLIB runs directory found; building fresh runs.")
    end

    # ── Phase 1: build baseline, enable STARLIB, run baseline ────────────────
    println("\n─── Phase 1: Baseline (starlib_option=$(opts.starlib_option)) ───")
    create_factored_runs(opts.nova; baseline_only=true, runs_name=opts.runs_name)
    enable_starlib_in_runs!(runs_dir, opts.starlib_option)

    println("Running baseline to establish reaction network ordering...")
    run_one_ppn(joinpath(runs_dir, "baseline", "ppn.exe"), runs_dir, logs_dir)

    opts.baseline_only && return

    # ── Phase 2: resolve indices from baseline networksetup.txt ──────────────
    # The baseline run generates a networksetup.txt that reflects the actual
    # reaction ordering for this STARLIB option (e.g. STL02 inserts ~42 extra
    # reactions that shift indices relative to old REACLIB).
    baseline_network = joinpath(runs_dir, "baseline", "networksetup.txt")
    println("\n─── Phase 2: Resolving reaction indices (starlib_option=$(opts.starlib_option)) ───")
    create_factored_runs(opts.nova;
        skip_baseline   = true,
        network_path    = baseline_network,
        runs_name       = opts.runs_name,
        prefer_sources  = ["STL02", "STL01"])
    enable_starlib_in_runs!(runs_dir, opts.starlib_option)

    if opts.build_only
        println("Build-only requested; not running factored ppn.exe runs.")
        return
    end

    # ── Phase 3: run factored runs in parallel ────────────────────────────────
    exes = list_ppn_executables(runs_dir)
    filter!(e -> basename(dirname(e)) != "baseline", exes)
    println("\nStarting $(length(exes)) STARLIB ppn run(s) with up to $(opts.jobs) concurrent job(s).")
    run_parallel(exes, runs_dir, logs_dir, opts.jobs)

    elapsed = round(Int, time() - global_start)
    println("==================================================")
    println("ALL REQUESTED STARLIB RUNS COMPLETE")
    println("Total time: $elapsed seconds")
    println("==================================================")
end

main()
