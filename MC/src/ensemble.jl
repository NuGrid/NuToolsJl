# ensemble.jl: gather the per-seed output of run_mc and load it as one ensemble.

"""
    collect_mc(out; io=stdout) -> number of samples

Gather the per-seed files of a [`run_mc`](@ref) output directory into
- `X.f64`: raw Float64, nsp × n column-major; column j = final mass fractions of sample j, rows as in
  `species.txt` (the first columns are cycle, time, t9, rho, ...);
- `xr.f64`: raw Float64, nrate × n; column j = x of every varied rate of sample j, rows as in `rates.txt`;
- `seeds.txt`: the seed of each sample, in column order.

In Python: `numpy.fromfile("X.f64").reshape(n, nsp)` (row j = sample j).
"""
function collect_mc(out::AbstractString; io::IO=stdout)
    files = sort(filter(f -> occursin(r"^s\d+\.f64$", f), readdir(joinpath(out, "final"))))
    isempty(files) && error("no samples in $out/final")
    nsp = length(readlines(joinpath(out, "species.txt")))
    nrate = count(l -> !startswith(l, '#'), readlines(joinpath(out, "rates.txt")))
    seeds = Int[]
    X, xr = Float64[], Float64[]
    for f in files
        xf = joinpath(out, "x", f)
        isfile(xf) || continue
        a = reinterpret(Float64, read(joinpath(out, "final", f)))
        b = reinterpret(Float64, read(xf))
        if length(a) != nsp || length(b) != nrate
            println(io, "$f: shape mismatch ($(length(a)) species, $(length(b)) rates), skipped")
            continue
        end
        push!(seeds, parse(Int, f[2:end-4])); append!(X, a); append!(xr, b)
    end
    write(joinpath(out, "X.f64"), X)
    write(joinpath(out, "xr.f64"), xr)
    write(joinpath(out, "seeds.txt"), join(seeds, '\n') * "\n")
    ff = joinpath(out, "failed.txt")
    nf = isfile(ff) ? length(readlines(ff)) : 0
    println(io, "$(length(seeds)) samples, $nsp columns, $nrate varied rates -> $out/{X.f64, xr.f64, seeds.txt} ",
            "($nf failed seeds listed in failed.txt)")
    length(seeds)
end

"One varied rate of an ensemble: its row in networksetup, the 5-character label, the reaction, and how its uncertainty was set (`\"\"` = its own, `\"default_fu\"`, `\"shared ...\"`)."
const RateInfo = NamedTuple{(:row, :label, :reaction, :how),Tuple{Int,String,String,String}}

"""
    MCEnsemble(out)

A collected Monte Carlo (see [`collect_mc`](@ref)):
- `species`: x-time.dat column names; `X`: nsp × n final values;
- `rates`: one [`RateInfo`](@ref) per varied rate; `xr`: nrate × n drawn x;
- `seeds`: seed of each sample.

Use [`abundance`](@ref), [`ratio`](@ref) and [`rank_reactions`](@ref) on it.
"""
struct MCEnsemble
    dir::String
    species::Vector{String}
    X::Matrix{Float64}
    rates::Vector{RateInfo}
    xr::Matrix{Float64}
    seeds::Vector{Int}
    columns::Dict{String,Vector{Int}}     # isotope name -> X rows
end

function MCEnsemble(out::AbstractString)
    isfile(joinpath(out, "X.f64")) || error("$out/X.f64 not found: run collect_mc first")
    species = readlines(joinpath(out, "species.txt"))
    rates = RateInfo[]
    for l in readlines(joinpath(out, "rates.txt"))
        startswith(l, '#') && continue
        t = split(l, '\t')
        push!(rates, (row=parse(Int, t[1]), label=String(t[2]), reaction=String(t[3]),
                      how=length(t) > 3 ? String(t[4]) : ""))
    end
    X = reshape(collect(reinterpret(Float64, read(joinpath(out, "X.f64")))), length(species), :)
    xr = reshape(collect(reinterpret(Float64, read(joinpath(out, "xr.f64")))), length(rates), :)
    sf = joinpath(out, "seeds.txt")
    seeds = isfile(sf) ? parse.(Int, filter(!isempty, readlines(sf))) : collect(1:size(X, 2))
    columns = Dict{String,Vector{Int}}()
    for (i, c) in enumerate(species)
        nm = iso_name(c)
        nm === nothing || push!(get!(columns, nm, Int[]), i)
    end
    MCEnsemble(abspath(out), species, X, rates, xr, seeds, columns)
end

nsamples(ens::MCEnsemble) = size(ens.X, 2)

function Base.show(io::IO, ens::MCEnsemble)
    print(io, "MCEnsemble(\"", ens.dir, "\", ", nsamples(ens), " samples, ", length(ens.rates), " varied rates)")
end

"""
    abundance(ens, iso; isomers=true) -> Vector{Float64}

Final mass fraction of `iso` (`"Al26"`, `"C12"`, ...) in every sample; with `isomers` the isomer
species is added (`"Al26"` = ground + `"Al26m"`). An isotope not in the network gives zeros.
"""
function abundance(ens::MCEnsemble, iso::AbstractString; isomers::Bool=true)
    idx = vcat(get(ens.columns, iso, Int[]), isomers ? get(ens.columns, iso * "m", Int[]) : Int[])
    isempty(idx) ? zeros(nsamples(ens)) : vec(sum(ens.X[idx, :]; dims=1))
end

"""
    ratio(ens, a, b; X0=nothing, xi=0) -> Vector{Float64}

Mass-fraction ratio `a/b` per sample (`NaN` where `b` is 0). With initial abundances `X0` and `xi > 0`
it is the mixed ratio (X_a + ξ X0_a) / (X_b + ξ X0_b) of Ward et al. (2025) Eq. 4.
"""
function ratio(ens::MCEnsemble, a::AbstractString, b::AbstractString; X0=nothing, xi::Real=0)
    xa, xb = abundance(ens, a), abundance(ens, b)
    if X0 !== nothing && xi > 0
        return (xa .+ xi * get(X0, a, 0.0)) ./ (xb .+ xi * get(X0, b, 0.0))
    end
    xa ./ ifelse.(xb .> 0, xb, NaN)
end

"`\"own\"`, `\"default f.u.\"` or `\"shared x\"` for a rate's `how` field"
uncertainty_source(how::AbstractString) =
    occursin("default_fu", how) ? "default f.u." : occursin("shared", how) ? "shared x" : "own"
