# run_parallel.jl — shared worker-pool helper for launching many ppn.exe runs.
#
# Not part of the NuGridJl package: this is orchestration (it launches an
# external process), which NuGridJl itself never does. `include` this from
# another tools/ script.

"""
    pin_single_threaded_blas!()

Force OMP/OpenBLAS/MKL/BLIS/vecLib to one thread each, in this process's
environment (inherited by every `ppn.exe` child it spawns). `ppn.exe` links
against a threaded BLAS; without this, running several instances at once
oversubscribes the CPU and each individual run gets slower, defeating the
point of `run_parallel`.
"""
function pin_single_threaded_blas!()
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
        ENV[var] = get(ENV, var, "1")
    end
end

"""
    run_parallel(dirs; jobs = 4, exe = "./ppn.exe", verbose = true) -> Vector{Tuple{String,Bool}}

Run `exe` in every directory in `dirs`, up to `jobs` at a time (concurrency is
via `@async` tasks, not OS threads — appropriate here since the work is
waiting on external processes, not CPU-bound Julia code, so this doesn't need
`julia -t N`). Pins BLAS to single-threaded first (see
[`pin_single_threaded_blas!`](@ref)). Each run's combined stdout/stderr is
written to `<dir>/run.log`. Returns `(dir, success)` pairs in `dirs` order.

With `verbose = true` (the default), prints a line to stdout each time a run
starts and each time one finishes (`[k/n]`, elapsed time, ok/FAILED), plus a
final summary line with total wall-clock time once every run is done —
useful for a sweep of this size (minutes to hours) where otherwise the CLI
gives no sign of progress until everything finishes at once.
"""
function run_parallel(dirs; jobs::Integer = 4, exe::AbstractString = "./ppn.exe",
                       verbose::Bool = true)
    pin_single_threaded_blas!()
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
