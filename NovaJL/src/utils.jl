using Printf

function factor_to_folder(f::Real)
    if isapprox(f, round(Int, f))
        return string(Int(round(f)))
    else
        return string(f)
    end
end

function factor_to_folder(f::AbstractString)
    label = strip(f)
    for prefix in ("factor_", "fact_")
        if startswith(label, prefix)
            return label[(lastindex(prefix) + 1):end]
        end
    end
    return label
end

function factor_to_scientific_folder(f::Real)
    return replace(@sprintf("%.3E", Float64(f)), "+" => "p", "-" => "m")
end

function parse_factor_folder_label(label)
    normalized = replace(strip(string(label)), "p" => "+", "m" => "-")
    return tryparse(Float64, normalized)
end

function factor_directory_labels(f::Real)
    return unique([factor_to_folder(f), factor_to_scientific_folder(f)])
end

function factor_directory_labels(f::AbstractString)
    label = factor_to_folder(f)
    parsed = parse_factor_folder_label(label)
    parsed === nothing && return [label]
    return unique([label, factor_to_scientific_folder(parsed), factor_to_folder(parsed)])
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
