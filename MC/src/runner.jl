# runner.jl: run a ppn reaction-rate Monte Carlo, one ppn.exe per seed, many at a time.

const PPN_INPUTS = ["ppn_frame.input", "ppn_physics.input", "ppn_solver.input"]

"""
    parse_seeds(s) -> Vector{Int}

`"1-100"`, `"1,5,9"` or a mix (`"1-100,200-250"`).
"""
function parse_seeds(s::AbstractString)
    out = Int[]
    for part in split(s, ',')
        part = strip(part)
        isempty(part) && continue
        if occursin('-', part)
            a, b = parse.(Int, split(part, '-'))
            append!(out, a:b)
        else
            push!(out, parse(Int, part))
        end
    end
    out
end

"""
    physics_with_seed(text, seed) -> String

A `ppn_physics.input` text with `factor_uncertainty_mc = .true.` and `factor_uncertainty_mc_seed = seed`
(existing lines for either are replaced), inserted before the namelist's closing `/`.
"""
function physics_with_seed(text::AbstractString, seed::Integer)
    lines = filter(l -> !occursin(r"^\s*factor_uncertainty_mc(_seed)?\s*="i, l), split(text, '\n'))
    i = findfirst(l -> strip(l) == "/", lines)
    i === nothing && error("ppn_physics.input: no closing '/' of the &ppn_physics namelist")
    splice!(lines, i:i-1, ["\tfactor_uncertainty_mc = .true.", "\tfactor_uncertainty_mc_seed = $seed"])
    join(lines, '\n')
end

function _setup_worker(template, wdir)
    run_dir = joinpath(wdir, "ppn")
    isdir(run_dir) && return run_dir
    mkpath(run_dir)
    symlink(realpath(joinpath(template, "..", "NPDATA")), joinpath(wdir, "NPDATA"))
    symlink(realpath(joinpath(template, "ppn.exe")), joinpath(run_dir, "ppn.exe"))
    frame = read(joinpath(template, "ppn_frame.input"), String)
    need = vcat(PPN_INPUTS, [namelist_value(frame, "trajectory_fn")],
                filter(f -> startswith(f, "isotopedatabase"), readdir(template)))
    ini = namelist_value(frame, "ini_filename")
    ini === nothing || push!(need, ini)
    for f in need
        f === nothing && continue
        cp(joinpath(template, f), joinpath(run_dir, f); force=true)
    end
    run_dir
end

"write a file atomically (other array tasks may read it at the same time)"
function _write_once(path, content)
    isfile(path) && return
    tmp = path * ".$(gethostname()).$(getpid())"
    write(tmp, content)
    mv(tmp, path; force=true)
end

"""
    run_seed(seed, template, out, wdir; timeout=1800) -> (reason, seconds)

Run one Monte Carlo sample in the worker directory `wdir` and store its result under `out`
(see [`run_mc`](@ref)). `reason` is `nothing` on success, otherwise why the seed failed.
"""
function run_seed(seed::Integer, template::AbstractString, out::AbstractString, wdir::AbstractString;
                  timeout::Real=1800)
    t0 = time()
    run_dir = _setup_worker(template, wdir)
    write(joinpath(run_dir, "ppn_physics.input"),
          physics_with_seed(read(joinpath(template, "ppn_physics.input"), String), seed))
    for f in readdir(run_dir)
        any(p -> startswith(f, p), ("iso_massf", "flux_", "x-time", "rate_uncertainty_")) && rm(joinpath(run_dir, f))
    end
    logf = joinpath(run_dir, "run.log")
    reason = nothing
    open(logf, "w") do lf
        p = run(pipeline(Cmd(`./ppn.exe`; dir=run_dir); stdout=lf, stderr=lf); wait=false)
        timer = Timer(timeout) do _
            process_running(p) && kill(p)
        end
        wait(p)
        close(timer)
        if p.termsignal != 0
            reason = "killed (signal $(p.termsignal); timeout $(timeout) s?)"
        elseif p.exitcode != 0
            reason = "exit $(p.exitcode)"
        end
    end
    if reason === nothing && !occursin("simulation complete", read(logf, String))
        reason = "no \"simulation complete\" in log"
    end
    local names, vals
    if reason === nothing
        names, vals = last_xtime_row(joinpath(run_dir, "x-time.dat"))
        length(vals) == length(names) || (reason = "incomplete last x-time.dat row")
    end
    if reason !== nothing
        mkpath(joinpath(out, "logs"))
        cp(logf, joinpath(out, "logs", @sprintf("s%06d.log", seed)); force=true)
        return reason, time() - t0
    end
    xs = read_rate_x(joinpath(run_dir, "rate_uncertainty_x.txt"))
    _write_once(joinpath(out, "species.txt"), join(names, '\n') * "\n")
    _write_once(joinpath(out, "rates.txt"),
                "# row\tlabel\treaction\thow (rate_uncertainty_x.txt of the first seed run)\n" *
                join(["$(r.row)\t$(r.label)\t$(r.reaction)\t$(strip(replace(r.how, "drawn" => "")))" for r in xs],
                     '\n') * "\n")
    write(joinpath(out, "x", @sprintf("s%06d.f64", seed)), Float64[r.x for r in xs])
    write(joinpath(out, "final", @sprintf("s%06d.f64", seed)), fnum.(vals))   # written last: marks the seed done
    nothing, time() - t0
end

"""
    run_mc(template, out, seeds; jobs=Sys.CPU_THREADS, timeout=1800, io=stdout) -> (ndone, nfailed)

Run a ppn reaction-rate Monte Carlo. `template` is a finished-setup ppn run directory (ppn.exe, the
`*.input` files, isotopedatabase, iniab, trajectory) whose parent holds the `NPDATA` link. Each seed is
one ppn run of that setup with `factor_uncertainty_mc = .true.` and `factor_uncertainty_mc_seed = seed`
added to its `ppn_physics.input`; `jobs` of them run at once, each killed after `timeout` seconds.

Per seed it writes, under `out`:
- `final/s<seed>.f64`: mass fractions of the last x-time.dat row, raw Float64;
- `x/s<seed>.f64`: the x drawn for every varied rate, in rate_uncertainty_x.txt row order, raw Float64;

and once `species.txt` (x-time.dat column names) and `rates.txt` (row, label, reaction, source of the
uncertainty). A failed seed is listed in `failed.txt` with its reason, and its log is kept in `logs/`.
Seeds whose `final/` file exists are skipped, so an interrupted run restarts with the same call, and
several processes (e.g. SLURM array tasks) can share one `out`. Gather the result with
[`collect_mc`](@ref).
"""
function run_mc(template::AbstractString, out::AbstractString, seeds::AbstractVector{<:Integer};
                jobs::Integer=Sys.CPU_THREADS, timeout::Real=1800, io::IO=stdout)
    template, out = abspath(template), abspath(out)
    for f in vcat(PPN_INPUTS, "ppn.exe")
        isfile(joinpath(template, f)) || error("$template: missing $f")
    end
    ispath(joinpath(template, "..", "NPDATA")) || error("$template/../NPDATA not found")
    mkpath(joinpath(out, "final")); mkpath(joinpath(out, "x"))

    todo = filter(s -> !isfile(joinpath(out, "final", @sprintf("s%06d.f64", s))), collect(seeds))
    tag = "$(gethostname())_$(getpid())"                  # unique across nodes and array tasks
    println(io, "$(length(todo)) seeds to run with $jobs parallel jobs -> $out"); flush(io)
    slots = Channel{Int}(jobs)
    foreach(i -> put!(slots, i), 1:jobs)
    t0, done, nfail = time(), Ref(0), Ref(0)
    lk = ReentrantLock()
    asyncmap(todo; ntasks=jobs) do seed
        slot = take!(slots)
        reason, _ = try
            run_seed(seed, template, out, joinpath(out, "work", tag, "w$slot"); timeout)
        catch err
            sprint(showerror, err), 0.0
        finally
            put!(slots, slot)
        end
        lock(lk) do
            done[] += 1
            if reason !== nothing
                nfail[] += 1
                open(f -> println(f, "$seed\t$reason"), joinpath(out, "failed.txt"), "a")
            end
            if done[] % max(1, length(todo) ÷ 50) == 0 || done[] == length(todo) || reason !== nothing
                el = time() - t0
                @printf(io, "%d/%d done, %d failed, %.1f min, ~%.1f min left%s\n", done[], length(todo), nfail[],
                        el / 60, el / done[] * (length(todo) - done[]) / 60,
                        reason === nothing ? "" : "  [seed $seed: $reason]")
                flush(io)
            end
        end
    end
    rm(joinpath(out, "work", tag); recursive=true, force=true)   # only this invocation's worker directories
    try                                                           # and work/ itself once no other task uses it
        isdir(joinpath(out, "work")) && isempty(readdir(joinpath(out, "work"))) && rm(joinpath(out, "work"))
    catch
    end
    done[] - nfail[], nfail[]
end
