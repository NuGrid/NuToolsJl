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
    run_parallel(dirs; jobs = 4, exe = "./ppn.exe") -> Vector{Tuple{String,Bool}}

Run `exe` in every directory in `dirs`, up to `jobs` at a time (concurrency is
via `@async` tasks, not OS threads — appropriate here since the work is
waiting on external processes, not CPU-bound Julia code, so this doesn't need
`julia -t N`). Pins BLAS to single-threaded first (see
[`pin_single_threaded_blas!`](@ref)). Each run's combined stdout/stderr is
written to `<dir>/run.log`. Returns `(dir, success)` pairs in `dirs` order.
"""
function run_parallel(dirs; jobs::Integer = 4, exe::AbstractString = "./ppn.exe")
    pin_single_threaded_blas!()
    results = fill(false, length(dirs))
    queue = Channel{Int}(length(dirs))
    foreach(i -> put!(queue, i), eachindex(dirs))
    close(queue)

    @sync for _ in 1:max(1, min(jobs, length(dirs)))
        @async for i in queue
            dir = dirs[i]
            log_path = joinpath(dir, "run.log")
            cmd = Cmd(`$exe`; dir = dir)
            results[i] = try
                success(pipeline(cmd; stdout = log_path, stderr = log_path))
            catch
                false
            end
        end
    end
    return collect(zip(dirs, results))
end
