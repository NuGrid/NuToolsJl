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
