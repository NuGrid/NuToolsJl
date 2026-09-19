# namelist_utils.jl — small helpers for patching Fortran namelist text
# (ppn_physics.input / ppn_frame.input), shared by build_decay_sweep.jl and
# build_sigma_sweep.jl. Not part of the NuGridJl package: like run_parallel.jl,
# this is orchestration/file-munging, `include`d directly by whichever tool
# needs it.

using Printf

"""
    update_namelist(text, replacements) -> String

Set each `name => value` pair in a Fortran namelist's text: overwrite the
value in place if `name` already has a line, otherwise insert it just before
the terminating `/`. Case-insensitive on `name` (Fortran namelists are).
"""
function update_namelist(text::AbstractString, replacements::Vector{<:Pair})
    lines = split(text, '\n')
    !isempty(lines) && lines[end] == "" && pop!(lines)

    replacement_keys = Dict(lowercase(k) => v for (k, v) in replacements)
    found = Set{String}()
    output = String[]

    for line in lines
        if strip(line) == "/"
            for (key, value) in replacements
                lowercase(key) in found || push!(output, "        $key = $value")
            end
            push!(output, line)
            continue
        end
        m = match(r"^(\s*)([A-Za-z_][A-Za-z0-9_]*)\s*=", line)
        if m !== nothing && haskey(replacement_keys, lowercase(m.captures[2]))
            key = m.captures[2]
            push!(found, lowercase(key))
            push!(output, "$(m.captures[1])$key = $(replacement_keys[lowercase(key)])")
        else
            push!(output, line)
        end
    end
    return join(output, '\n') * "\n"
end

fortran_double(value::Real) = replace(@sprintf("%.10E", value), "E" => "d")

"""
    write_rate_factors!(ppn_physics_input_path, index_factor_pairs)

Insert `rate_index(i) = <index>` / `rate_factor(i) = <factor>` for each
`(index, factor)` pair, just before the `&ppn_physics` namelist's closing
`/`. Slots are numbered from 1 (NuPPN supports up to `num_rate_factors = 10`).
This is ppn's own runtime rate-factor mechanism (`apply_rate_factors` in
`evaluate_rates.F90`): a single flat, T9-independent multiplicative scalar
applied to whatever the reaction evaluates to at each timestep, addressed by
reaction index -- works for any rate source (STARLIB, NACRE, JINA, weak-rate
tables, ...), unlike `build_sigma_sweep.jl`'s STARLIB-file-rewriting sweep,
which only works for STARLIB-sourced reactions.
"""
function write_rate_factors!(ppn_physics_input_path::AbstractString, index_factor_pairs)
    lines = readlines(ppn_physics_input_path)
    terminator = findfirst(l -> strip(l) == "/", lines)
    terminator === nothing && throw(ArgumentError(
        "no namelist terminator '/' found in $ppn_physics_input_path"))
    new_lines = String[]
    for (i, (index, factor)) in enumerate(index_factor_pairs)
        push!(new_lines, "        rate_index($i) = $index")
        push!(new_lines, "        rate_factor($i) = $(factor)")
    end
    splice!(lines, terminator:(terminator - 1), new_lines)
    write(ppn_physics_input_path, join(lines, "\n") * "\n")
end
