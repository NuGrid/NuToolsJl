#!/usr/bin/env julia
# collect_vital_rates.jl
#
# Run PPN with print_rate = .true. for each of the 117 VITAL reactions
# (and their preferred non-VITAL counterparts) to collect rate(T) curves.
#
# Output:
#   <nova_case>/VITAL_rate_comp/VITAL_rates/rate_<name>.txt
#   <nova_case>/VITAL_rate_comp/non_VITAL_rates/rate_<name>.txt
#
# Usage:
#   julia tools/collect_vital_rates.jl --nova ne_nova_1.15_12_X_weiss_mixed

include(joinpath(@__DIR__, "NovaRunTools.jl"))
using .NovaRunTools
using Printf

# ─── Options ──────────────────────────────────────────────────────────────────

Base.@kwdef mutable struct CollectOptions
    nova::String             = "nova_test"
    jobs::Int                = get(ENV, "MAX_JOBS", nothing) === nothing ? 4 : parse(Int, ENV["MAX_JOBS"])
    vital_run::String        = "runs_VITAL"
    baseline_run::String     = "runs"
    vital_only::Bool         = false
end

function usage()
    println("""
Usage:
  julia tools/collect_vital_rates.jl [options]

Options:
  --nova NAME            Nova directory under nova_cases/ (default: nova_test)
  --jobs N               Parallel jobs (default: MAX_JOBS or 4)
  --vital-run DIR        Sweep dir containing VITAL baseline/ (default: runs_VITAL)
  --baseline-run DIR     Sweep dir containing non-VITAL baseline/ (default: runs)
  --vital-only           Collect only VITAL rates, skip non-VITAL
  -h, --help             Show this help

Output:
  <nova_case>/VITAL_rate_comp/VITAL_rates/rate_<name>.txt
  <nova_case>/VITAL_rate_comp/non_VITAL_rates/rate_<name>.txt

  Both directories sit alongside VITAL_rate_comp/non_VITAL_rates/.
  Reactions with no non-VITAL alternative (pp-chain, some compound reactions)
  are skipped with a warning for the non-VITAL set.
""")
end

function parse_args(args)
    opts = CollectOptions()
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--nova"
            i += 1; i <= length(args) || error("--nova requires a value")
            opts.nova = args[i]
        elseif arg == "--jobs"
            i += 1; i <= length(args) || error("--jobs requires a value")
            opts.jobs = parse(Int, args[i])
            opts.jobs >= 1 || error("--jobs must be >= 1")
        elseif arg == "--vital-run"
            i += 1; i <= length(args) || error("--vital-run requires a value")
            opts.vital_run = args[i]
        elseif arg == "--baseline-run"
            i += 1; i <= length(args) || error("--baseline-run requires a value")
            opts.baseline_run = args[i]
        elseif arg == "--vital-only"
            opts.vital_only = true
        elseif arg in ("-h", "--help")
            usage(); exit(0)
        else
            error("Unknown argument: $arg")
        end
        i += 1
    end
    return opts
end

# ─── Reaction naming ──────────────────────────────────────────────────────────

const RTYPE_CODES = Dict(
    "(p,g)" => "pg",  "(p,a)" => "pa",  "(a,g)" => "ag",  "(a,n)" => "an",
    "(a,p)" => "ap",  "(g,p)" => "gp",  "(g,a)" => "ga",  "(g,n)" => "gn",
    "(n,g)" => "ng",  "(n,p)" => "np",  "(n,a)" => "na",  "(p,n)" => "pn",
    "(v,v)" => "vv",  "(+,g)" => "bpg", "(-,g)" => "bmg",
    "(b,a)" => "ba",  "(b,n)" => "bn",  "(e,v)" => "ev",
)

function fmt_sp(sp)
    (sp === nothing || sp == (0, "G")) && return nothing
    "$(sp[1])$(titlecase(lowercase(sp[2])))"
end

function rate_name(row::NovaRunTools.Row, idx::Int)
    # Use the standard row_name for the five supported channels (pg, pa, ag, an, ap).
    n = NovaRunTools.row_name(row)
    n !== nothing && return n

    # Fallback: reactant + type-code + first non-empty product
    type_code = get(RTYPE_CODES, row.rtype,
        replace(replace(replace(row.rtype, "(" => ""), ")" => ""), "," => ""))

    r1 = fmt_sp(row.reactant)
    p1 = fmt_sp(row.product_1)
    p1 === nothing && (p1 = fmt_sp(row.product_2))

    parts = filter(!isnothing, [r1, type_code, p1])
    isempty(parts) && return @sprintf("rxn%03d", idx)
    return join(parts, "_")
end

# Zero-padded index prefix ensures lexicographic order matches network order.
rate_file_stem(row::NovaRunTools.Row, idx::Int) = @sprintf("%03d_%s", idx, rate_name(row, idx))

# ─── ppn_frame.input editing ──────────────────────────────────────────────────

function write_frame_input(src_path, dst_path, which_rate_idx)
    lines = readlines(src_path, keep=true)
    buf = IOBuffer()
    in_namelist = false
    for line in lines
        stripped = strip(line)
        stripped == "&ppn_frame" && (in_namelist = true)
        # Drop any pre-existing print_rate / which_rate entries
        if in_namelist && (occursin(r"^\s*print_rate\s*=", line) || occursin(r"^\s*which_rate\s*=", line))
            continue
        end
        # Inject before the closing slash
        if in_namelist && stripped == "/"
            write(buf, "   print_rate = .true.\n")
            write(buf, "   which_rate = $which_rate_idx\n")
            in_namelist = false
        end
        write(buf, line)
    end
    write(dst_path, String(take!(buf)))
end

# ─── Single rate collection run ───────────────────────────────────────────────

function collect_one_rate(source_dir, which_rate_idx, out_path, work_dir, run_label)
    run_dir = joinpath(work_dir, run_label)
    ispath(run_dir) && rm(run_dir; recursive=true)
    mkpath(run_dir)

    try
        # Copy all source files; symlinks are reproduced as symlinks (Julia default).
        for entry in readdir(source_dir; join=false)
            src = joinpath(source_dir, entry)
            dst = joinpath(run_dir, entry)
            entry == "ppn_frame.input" && continue   # written below
            entry == "rate.txt"        && continue   # output file
            if islink(src)
                symlink(realpath(src), dst)
            elseif isdir(src)
                cp(src, dst; force=true)
            else
                cp(src, dst; force=true)
            end
        end

        # Write modified ppn_frame.input.
        write_frame_input(
            joinpath(source_dir, "ppn_frame.input"),
            joinpath(run_dir,    "ppn_frame.input"),
            which_rate_idx,
        )

        # Ensure ../NPDATA exists from within run_dir (vital.F90 opens "../NPDATA/…").
        # run_dir parent = work_dir, so work_dir/NPDATA must exist.
        npdata_target = realpath(joinpath(source_dir, "NPDATA"))
        parent_npdata = joinpath(work_dir, "NPDATA")
        ispath(parent_npdata) || symlink(npdata_target, parent_npdata)

        # Run PPN — it writes rate.txt and stops immediately.
        logfile = joinpath(work_dir, "$(run_label).log")
        open(logfile, "w") do io
            proc = run(pipeline(Cmd(`./ppn.exe`, dir=run_dir), stdout=io, stderr=io), wait=false)
            wait(proc)
            success(proc) || error("ppn.exe failed; see $logfile")
        end

        # Collect output.
        rate_src = joinpath(run_dir, "rate.txt")
        isfile(rate_src) || error("rate.txt was not produced (check $logfile)")
        cp(rate_src, out_path; force=true)
        println("  ✓ $(basename(out_path))")

    finally
        ispath(run_dir) && rm(run_dir; recursive=true)
    end
end

# ─── Parallel runner ──────────────────────────────────────────────────────────

# Each task: (source_dir, which_rate_idx, out_path, run_label)
function run_all_parallel(tasks, work_dir, jobs)
    isempty(tasks) && return

    queue    = Channel{Any}(length(tasks))
    for t in tasks; put!(queue, t); end
    close(queue)

    failures = Channel{Any}(length(tasks))
    workers  = Task[]
    for _ in 1:min(jobs, length(tasks))
        push!(workers, @async begin
            for (src_dir, idx, out_path, label) in queue
                try
                    collect_one_rate(src_dir, idx, out_path, work_dir, label)
                catch err
                    put!(failures, (label, err))
                end
            end
        end)
    end

    foreach(wait, workers)
    close(failures)

    failed = collect(failures)
    if !isempty(failed)
        println("\nFailures:")
        for (label, err) in failed
            println("  $label: $err")
        end
        error("$(length(failed)) rate collection job(s) failed")
    end
end

# ─── Main ─────────────────────────────────────────────────────────────────────

function main()
    opts = parse_args(ARGS)

    base             = NovaRunTools.nova_dir(opts.nova)
    vital_basedir    = joinpath(base, opts.vital_run,    "baseline")
    baseline_basedir = joinpath(base, opts.baseline_run, "baseline")
    output_dir       = joinpath(base, "VITAL_rate_comp")
    vital_out        = joinpath(output_dir, "VITAL_rates")
    nonvital_out     = joinpath(output_dir, "non_VITAL_rates")
    work_dir         = joinpath(output_dir, "_work")

    isdir(vital_basedir) || error(
        "VITAL baseline not found: $vital_basedir\n" *
        "  Run the VITAL sweep first or use --vital-run to specify its directory.")
    isdir(baseline_basedir) || error(
        "Non-VITAL baseline not found: $baseline_basedir\n" *
        "  Run the regular sweep first or use --baseline-run to specify its directory.")

    mkpath(vital_out)
    mkpath(nonvital_out)
    mkpath(work_dir)

    # ── Parse networks ────────────────────────────────────────────────────────
    println("Parsing VITAL networksetup from $vital_basedir...")
    vital_net  = NovaRunTools.parse_networksetup(joinpath(vital_basedir,    "networksetup.txt"))
    vital_rows = sort(filter(r -> r.active && r.index <= 117, vital_net), by = r -> r.index)
    length(vital_rows) == 117 || @warn "Expected 117 active VITAL rows, got $(length(vital_rows))"

    println("Parsing baseline networksetup from $baseline_basedir...")
    base_net = NovaRunTools.parse_networksetup(joinpath(baseline_basedir, "networksetup.txt"))

    # Active rows at index > 117 in the regular baseline are the non-VITAL alternatives.
    # First encountered active row wins for each reaction_key (prefer lower index = REACLIB default ordering).
    nonvital_by_key = Dict{Any, NovaRunTools.Row}()
    for row in base_net
        row.active       || continue
        row.index <= 117 && continue     # VITAL block — all F/VITAL, skip
        key = NovaRunTools.reaction_key(row)
        haskey(nonvital_by_key, key) || (nonvital_by_key[key] = row)
    end

    # ── Build task lists ──────────────────────────────────────────────────────
    vital_tasks    = Any[]
    nonvital_tasks = Any[]
    n_missing = 0

    println("\nReactions:")
    for row in vital_rows
        stem = rate_file_stem(row, row.index)
        key  = NovaRunTools.reaction_key(row)

        push!(vital_tasks, (
            vital_basedir,
            row.index,
            joinpath(vital_out, "rate_$(stem).txt"),
            "v_$(stem)",
        ))

        if !opts.vital_only
            if haskey(nonvital_by_key, key)
                nv = nonvital_by_key[key]
                push!(nonvital_tasks, (
                    baseline_basedir,
                    nv.index,
                    joinpath(nonvital_out, "rate_$(stem).txt"),
                    "n_$(stem)",
                ))
                println("  [$(row.index)] $(stem)  VITAL→idx $(row.index)  non-VITAL→idx $(nv.index) ($(nv.source))")
            else
                println("  [$(row.index)] $(stem)  VITAL only (no non-VITAL active row in baseline)")
                n_missing += 1
            end
        else
            println("  [$(row.index)] $(stem)")
        end
    end

    # ── Run ──────────────────────────────────────────────────────────────────
    println("\n─── Collecting VITAL rates ($(length(vital_tasks)) runs, $(opts.jobs) parallel) ───")
    run_all_parallel(vital_tasks, work_dir, opts.jobs)

    if !opts.vital_only && !isempty(nonvital_tasks)
        println("\n─── Collecting non-VITAL rates ($(length(nonvital_tasks)) runs, $(opts.jobs) parallel) ───")
        run_all_parallel(nonvital_tasks, work_dir, opts.jobs)
    end

    # Cleanup work dir (NPDATA symlink + empty _work/).
    try; rm(work_dir; recursive=true); catch; end

    println("\n══════════════════════════════════════════")
    println("DONE")
    println("  VITAL rates   : $(length(vital_tasks)) files in $vital_out")
    opts.vital_only || println("  non-VITAL rates: $(length(nonvital_tasks)) files in $nonvital_out")
    n_missing > 0   && println("  Note: $n_missing reactions had no non-VITAL counterpart (expected for pp-chain, compound (v,v), some decays).")
    println("══════════════════════════════════════════")
end

main()
