"""
    pin_single_threaded_openmp!()

Force OMP/OpenBLAS/MKL to one thread each, in this process's environment
(inherited by every MESA child process it spawns). MESA's `star` binary is
built against OpenMP (matrix solves, EOS/opacity table interpolation);
without this, running several cases at once oversubscribes the CPU and each
individual run gets slower, defeating the point of [`run_mesa_parallel`](@ref).
"""
function pin_single_threaded_openmp!()
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
        ENV[var] = get(ENV, var, "1")
    end
end

"""
    run_mesa_parallel(dirs; jobs = 4, exe = "./rn", verbose = true) -> Vector{Tuple{String,Bool}}

Run `exe` (MESA's test_suite launcher script by default) in every directory
in `dirs`, up to `jobs` at a time. A single MESA evolution is inherently
sequential in time — you can't parallelize *one* run's timestepping — so the
speedup here comes from running several independent cases (different WD
mass, accretion rate, composition, ...) concurrently instead, bounded by
however many cores you can spare.

Concurrency is via `@async` tasks, not OS threads (appropriate since the
work is waiting on external processes, not CPU-bound Julia code — this
doesn't need `julia -t N`). Pins OpenMP to single-threaded first (see
[`pin_single_threaded_openmp!`](@ref)) so `jobs` concurrent MESA processes
don't each also try to spread across every core. Each run's combined
stdout/stderr is written to `<dir>/run.log` (overwriting any existing one).
Returns `(dir, success)` pairs in `dirs` order.

Each directory in `dirs` must already be a ready-to-run MESA case (compiled
via `./mk`, inlists already set the way you want them) — this only launches
`exe`, it doesn't build anything or edit inlists for you.

With `verbose = true` (the default), prints a line each time a run starts
and each time one finishes (`[k/n]`, elapsed time, ok/FAILED), plus a final
summary with total wall-clock time — useful since individual MESA cases can
take tens of minutes and the CLI would otherwise give no sign of progress.
"""
function run_mesa_parallel(dirs; jobs::Integer = 4, exe::AbstractString = "./rn",
                            verbose::Bool = true)
    pin_single_threaded_openmp!()
    n = length(dirs)
    results = fill(false, n)
    queue = Channel{Int}(n)
    foreach(i -> put!(queue, i), eachindex(dirs))
    close(queue)

    sweep_start = time()
    progress_lock = ReentrantLock()
    started = 0
    finished = 0

    @sync for _ in 1:max(1, min(jobs, n))
        @async for i in queue
            dir = dirs[i]
            if verbose
                lock(progress_lock) do
                    started += 1
                    println("[$started/$n started, $finished/$n done] running $(basename(dir))")
                end
            end
            run_start = time()
            log_path = joinpath(dir, "run.log")
            cmd = Cmd(`$exe`; dir = dir)
            ok = try
                success(pipeline(cmd; stdout = log_path, stderr = log_path))
            catch
                false
            end
            results[i] = ok
            if verbose
                run_elapsed = round(time() - run_start; digits = 1)
                lock(progress_lock) do
                    finished += 1
                    status = ok ? "ok" : "FAILED"
                    println("[$finished/$n done] $(basename(dir)) $status ($(run_elapsed)s)")
                end
            end
        end
    end

    if verbose
        total_elapsed = round(time() - sweep_start; digits = 1)
        n_failed = count(!, results)
        println("sweep finished: $n runs, $(n - n_failed) ok, $n_failed failed, " *
                "$(total_elapsed)s total ($(jobs) jobs)")
    end

    return collect(zip(dirs, results))
end
