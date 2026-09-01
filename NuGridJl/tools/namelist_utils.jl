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
