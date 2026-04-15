const ELEMENT_Z = Dict(
    "H"=>1, "HE"=>2, "LI"=>3, "BE"=>4, "B"=>5,
    "C"=>6, "N"=>7, "O"=>8, "F"=>9, "NE"=>10,
    "NA"=>11, "MG"=>12, "AL"=>13, "SI"=>14,
    "P"=>15, "S"=>16, "CL"=>17, "AR"=>18,
    "K"=>19, "CA"=>20
)

function parse_isotope(name::String)

    name = uppercase(strip(name))

    if name == "N" || name == "NEUT"
        return 0, 1, 1
    elseif name == "P" || name == "PROT"
        return 1, 0, 1
    end

    parts = split(name, "-")
    length(parts) < 2 && return nothing

    symbol = parts[1]
    A = parse(Int, parts[2])

    if !haskey(ELEMENT_Z, symbol)
        @warn "Unknown element: $symbol - skipping"
        return nothing
    end

    Z = ELEMENT_Z[symbol]
    N = A - Z

    return Z, N, A
end


function factor_to_folder(f)
    if isapprox(f, round(Int, f))
        return string(Int(round(f)))
    else
        return string(f)
    end
end


function dots_to_missing!(df)
    for col in names(df)[3:end]
        df[!, col] = [
            ismissing(x) ? missing :
            x == "..." ? missing :
            x isa Number ? x :
            parse(Float64, x)
            for x in df[!, col]
        ]
    end
end
