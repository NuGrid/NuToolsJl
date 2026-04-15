using DataFrames

function read_initial_abundance(filepath)

    isotopes = String[]
    values   = Float64[]

    for line in eachline(filepath)

        line = strip(line)
        isempty(line) && continue

        parts = split(line)

        try
            Z = parse(Int, parts[1])

            if uppercase(parts[2]) == "PROT"
                isotope = "H-1"
                X = parse(Float64, parts[end])

            else
                raw = lowercase(parts[2])

                if length(parts) >= 3 && occursin(r"^\d+$", parts[3])
                    symbol = uppercase(raw)
                    A = parse(Int, parts[3])
                else
                    m = match(r"([a-z]+)(\d+)", raw)
                    m === nothing && continue

                    symbol = uppercase(m.captures[1])
                    A = parse(Int, m.captures[2])
                end

                X = parse(Float64, parts[end])
                isotope = symbol * "-" * string(A)
            end

            push!(isotopes, isotope)
            push!(values, X)

        catch
            continue
        end
    end

    return DataFrame(isotope=isotopes, X=values)
end


function read_iso_massf(filepath)

    isotopes = String[]
    values   = Float64[]

    for line in eachline(filepath)

        line = strip(line)

        isempty(line) && continue
        startswith(line, "#") && continue
        startswith(line, "H NUM") && continue

        parts = split(line)

        if length(parts) >= 6
            try
                X = parse(Float64, parts[5])
                raw = parts[end]

                if raw == "PROT"
                    isotope = "H-1"
                elseif raw == "NEUT"
                    isotope = "n"
                elseif length(parts) >= 7
                    isotope = uppercase(parts[6]) * "-" * parts[7]
                else
                    m = match(r"([A-Za-z]+)(\d+)", raw)
                    m === nothing && continue
                    isotope = uppercase(m.captures[1]) * "-" * m.captures[2]
                end

                push!(isotopes, isotope)
                push!(values, X)

            catch
                continue
            end
        end
    end

    return DataFrame(isotope=isotopes, X=values)
end


function read_trajectory(filepath)

    time = Float64[]
    T    = Float64[]
    rho  = Float64[]

    for line in eachline(filepath)

        line = strip(line)
        isempty(line) && continue
        startswith(line, "#") && continue

        parts = split(line)
        length(parts) < 3 && continue

        try
            push!(time, parse(Float64, parts[1]))
            push!(T,    parse(Float64, parts[2]))
            push!(rho,  parse(Float64, parts[3]))
        catch
            continue
        end
    end

    return DataFrame(time=time, T=T, rho=rho)
end
