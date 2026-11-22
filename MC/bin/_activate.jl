# Shared by the bin/ scripts: activate the MC package (instantiating it on first use) and load it.
import Pkg
let dir = normpath(joinpath(@__DIR__, ".."))
    Pkg.activate(dir; io=devnull)
    isfile(joinpath(dir, "Manifest.toml")) || Pkg.instantiate(; io=devnull)
end
using MC

"(options Dict for `--key value`, positional arguments)"
function cli_args(args, defaults=Dict{String,String}())
    o, pos = copy(defaults), String[]
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            i < length(args) || error("$(args[i]) needs a value")
            o[args[i][3:end]] = args[i+1]; i += 2
        else
            push!(pos, args[i]); i += 1
        end
    end
    o, pos
end
