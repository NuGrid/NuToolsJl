# run.jl — PPNRun, the single entry point to a nuppn output directory.
#
# One lazy handle over a run folder: it discovers the available cycles up front
# but only reads (and caches) an abundance vector, flux file, the time series,
# the network or the input decks when you ask for them.

"""
    PPNRun(dir)

A handle to a single-zone `nuppn` run.  `dir` may be the run folder itself or a
parent that contains a `ppn/` subfolder; whichever holds the `iso_massf*.DAT`
files is used.  Reads are lazy and cached.

Accessors: [`abundances`](@ref)`(run, cycle | :initial | :final | :decay)`,
[`fluxes`](@ref)`(run, cycle)`, [`xtime`](@ref)`(run)`, [`network`](@ref)`(run)`,
[`inputs`](@ref)`(run)`.
"""
struct PPNRun
    dir::String
    cycles::Vector{Int}
    _abund::Dict{Int,Abundances}
    _flux::Dict{Int,DataFrame}
    _cache::Dict{Symbol,Any}
end

const _ISO_RE = r"^iso_massf(\d+)\.DAT$"i

function _resolve_dir(dir)
    has_cycles(d) = isdir(d) && any(occursin(_ISO_RE, f) for f in readdir(d))
    has_cycles(dir) && return dir
    sub = joinpath(dir, "ppn")
    has_cycles(sub) && return sub
    return dir
end

function _discover_cycles(dir)
    cycles = Int[]
    isdir(dir) || return cycles
    for f in readdir(dir)
        m = match(_ISO_RE, f)
        m === nothing || push!(cycles, parse(Int, m.captures[1]))
    end
    return sort!(cycles)
end

function PPNRun(dir::AbstractString)
    resolved = _resolve_dir(String(dir))
    cycles = _discover_cycles(resolved)
    return PPNRun(resolved, cycles, Dict{Int,Abundances}(), Dict{Int,DataFrame}(), Dict{Symbol,Any}())
end

function Base.show(io::IO, run::PPNRun)
    print(io, "PPNRun(\"", run.dir, "\", ", length(run.cycles), " cycles")
    isempty(run.cycles) || print(io, " ", first(run.cycles), "–", last(run.cycles))
    print(io, ")")
end

_cycle_path(run, cycle) = joinpath(run.dir, string("iso_massf", lpad(cycle, 5, '0'), ".DAT"))
_flux_path(run, cycle) = joinpath(run.dir, string("flux_", lpad(cycle, 5, '0'), ".DAT"))

function _resolve_cycle(run::PPNRun, cycle::Symbol)
    isempty(run.cycles) && throw(ArgumentError("run has no cycles: $(run.dir)"))
    cycle === :final && return last(run.cycles)
    cycle === :initial && return first(run.cycles)
    throw(ArgumentError("unknown cycle selector :$cycle (use :initial, :final or :decay)"))
end

"""
    abundances(run, cycle) -> Abundances

`cycle` is a cycle number, `:initial`/`:final` for the first/last available
cycle, or `:decay` to read `iso_massfdecay.DAT` from the run directory.
"""
function abundances(run::PPNRun, cycle::Integer)
    get!(run._abund, Int(cycle)) do
        path = _cycle_path(run, cycle)
        isfile(path) || throw(ArgumentError(
            "cycle $cycle not found ($path); available: $(run.cycles)"))
        read_abundances(path)
    end
end

function abundances(run::PPNRun, cycle::Symbol)
    if cycle === :decay
        return get!(run._cache, :decay) do
            path = joinpath(run.dir, "iso_massfdecay.DAT")
            isfile(path) || throw(ArgumentError("no iso_massfdecay.DAT in $(run.dir)"))
            read_abundances(path)
        end
    end
    return abundances(run, _resolve_cycle(run, cycle))
end

"""
    fluxes(run, cycle) -> DataFrame

Reaction fluxes for `cycle` (a number or `:final`/`:initial`).
"""
function fluxes(run::PPNRun, cycle::Integer)
    get!(run._flux, Int(cycle)) do
        path = _flux_path(run, cycle)
        isfile(path) || throw(ArgumentError(
            "flux file for cycle $cycle not found ($path)"))
        read_fluxes(path)
    end
end
fluxes(run::PPNRun, cycle::Symbol) = fluxes(run, _resolve_cycle(run, cycle))

"The `x-time.dat` time series for the run (cached)."
function xtime(run::PPNRun)
    get!(run._cache, :xtime) do
        path = joinpath(run.dir, "x-time.dat")
        isfile(path) || throw(ArgumentError("no x-time.dat in $(run.dir)"))
        read_xtime(path)
    end
end

"The reaction `Network` for the run, from `networksetup.txt` (cached)."
function network(run::PPNRun)
    get!(run._cache, :network) do
        path = joinpath(run.dir, "networksetup.txt")
        isfile(path) || throw(ArgumentError("no networksetup.txt in $(run.dir)"))
        read_network(path)
    end
end

"The isotope database for the run, from `isotopedatabase.txt` (cached)."
function isotopedatabase(run::PPNRun)
    get!(run._cache, :isotopedatabase) do
        path = joinpath(run.dir, "isotopedatabase.txt")
        isfile(path) || throw(ArgumentError("no isotopedatabase.txt in $(run.dir)"))
        read_isotopedatabase(path)
    end
end

"""
    inputs(run) -> NamedTuple

The parsed ppn input decks as `(frame, physics, solver)`; a group whose file is
absent is an empty `Dict`.
"""
function inputs(run::PPNRun)
    get!(run._cache, :inputs) do
        read_group(name) = begin
            path = joinpath(run.dir, name)
            isfile(path) ? read_namelist(path) : Dict{Symbol,Any}()
        end
        (frame = read_group("ppn_frame.input"),
         physics = read_group("ppn_physics.input"),
         solver = read_group("ppn_solver.input"))
    end
end
