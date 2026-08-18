# runset.jl — collections of PPNRuns: factored-rate sweeps and Monte Carlo
# ensembles.
#
# Both are read from directory conventions on disk, built by tools outside
# this package (NuGridJl never launches ppn.exe itself): a `baseline/` run,
# plus either `<reaction>/fact_<factor>/` sweep runs or `sample_NNNNNN/`
# Monte Carlo runs each carrying an `mc_manifest.json`. Runs are discovered by
# directory convention alone, so no external manifest is required to read a
# sweep back.

"One factored run in a [`PPNSweep`](@ref): the swept reaction, the rate multiplier `factor`, and the run itself."
struct FactoredRun
    reaction::String
    factor::Float64
    run::PPNRun
end

"""
    PPNSweep

A directory of single-reaction-at-a-time rate sweeps: a `baseline` run plus one
[`FactoredRun`](@ref) per `<reaction>/fact_<factor>/` subdirectory. Use
[`reactions`](@ref), [`factors`](@ref) and [`sweep_run`](@ref) to navigate it.
"""
struct PPNSweep
    dir::String
    baseline::PPNRun
    runs::Vector{FactoredRun}
end

const _FACT_DIR_RE = r"^fact_([0-9.eE+-]+)$"

"""
    PPNSweep(dir) -> PPNSweep

Discover a sweep at `dir`: a `baseline/` run and every `<reaction>/fact_<factor>/`
subdirectory beneath it.
"""
function PPNSweep(dir::AbstractString)
    baseline_dir = joinpath(dir, "baseline")
    isdir(baseline_dir) || throw(ArgumentError("no baseline/ run in sweep directory $dir"))
    baseline = PPNRun(baseline_dir)

    runs = FactoredRun[]
    for reaction in sort(readdir(dir))
        reaction in ("baseline", "NPDATA") && continue
        reaction_dir = joinpath(dir, reaction)
        isdir(reaction_dir) || continue
        for entry in sort(readdir(reaction_dir))
            m = match(_FACT_DIR_RE, entry)
            m === nothing && continue
            factor = tryparse(Float64, m.captures[1])
            factor === nothing && continue
            push!(runs, FactoredRun(reaction, factor, PPNRun(joinpath(reaction_dir, entry))))
        end
    end
    return PPNSweep(dir, baseline, runs)
end

function Base.show(io::IO, sweep::PPNSweep)
    print(io, "PPNSweep(\"", sweep.dir, "\", ", length(reactions(sweep)), " reactions, ",
          length(sweep.runs), " factored runs)")
end

"Reaction names present in `sweep`, sorted and de-duplicated."
reactions(sweep::PPNSweep) = sort(unique(r.reaction for r in sweep.runs))

"Factors tested for `reaction` in `sweep`, sorted."
factors(sweep::PPNSweep, reaction::AbstractString) =
    sort([r.factor for r in sweep.runs if r.reaction == reaction])

"""
    sweep_run(sweep, reaction, factor) -> PPNRun

The factored run for `reaction` at `factor` (exact match — see [`factors`](@ref)
for what's available).
"""
function sweep_run(sweep::PPNSweep, reaction::AbstractString, factor::Real)
    idx = findfirst(r -> r.reaction == reaction && r.factor == factor, sweep.runs)
    idx === nothing && throw(ArgumentError(
        "no run for reaction=$reaction factor=$factor in $(sweep.dir); " *
        "available factors: $(factors(sweep, reaction))"))
    return sweep.runs[idx].run
end

# ---------------------------------------------------------------------------
# Monte Carlo ensembles
# ---------------------------------------------------------------------------

"One Monte Carlo draw in a [`PPNEnsemble`](@ref): its `id`, parsed `mc_manifest.json`, and the run itself."
struct MCSample
    id::Int
    manifest::Dict{String,Any}
    run::PPNRun
end

"""
    PPNEnsemble

A Monte Carlo ensemble of nuppn runs: a `baseline` run plus every
`sample_NNNNNN/` subdirectory that carries an `mc_manifest.json`.
"""
struct PPNEnsemble
    dir::String
    baseline::PPNRun
    samples::Vector{MCSample}
end

const _SAMPLE_DIR_RE = r"^sample_(\d+)$"

"""
    PPNEnsemble(dir) -> PPNEnsemble

Discover a Monte Carlo ensemble at `dir`: a `baseline/` run and every
`sample_NNNNNN/` subdirectory that carries an `mc_manifest.json`.
"""
function PPNEnsemble(dir::AbstractString)
    baseline_dir = joinpath(dir, "baseline")
    isdir(baseline_dir) || throw(ArgumentError("no baseline/ run in ensemble directory $dir"))
    baseline = PPNRun(baseline_dir)

    samples = MCSample[]
    for entry in sort(readdir(dir))
        m = match(_SAMPLE_DIR_RE, entry)
        m === nothing && continue
        sample_dir = joinpath(dir, entry)
        manifest_path = joinpath(sample_dir, "mc_manifest.json")
        isfile(manifest_path) || continue
        manifest = JSON.parsefile(manifest_path)
        push!(samples, MCSample(parse(Int, m.captures[1]), manifest, PPNRun(sample_dir)))
    end
    return PPNEnsemble(dir, baseline, samples)
end

function Base.show(io::IO, ens::PPNEnsemble)
    print(io, "PPNEnsemble(\"", ens.dir, "\", ", length(ens.samples), " samples)")
end

"The sample with id `id` in `ens`."
function sample(ens::PPNEnsemble, id::Integer)
    idx = findfirst(s -> s.id == id, ens.samples)
    idx === nothing && throw(ArgumentError("no sample $id in $(ens.dir)"))
    return ens.samples[idx]
end
