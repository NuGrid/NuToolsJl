function read_iso_massf(filepath)
    isotopes = String[]
    values = Float64[]

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
    temperature = Float64[]
    density = Float64[]

    ageunit = "SEC"
    tunit = "T9K"
    rhounit = "CGS"

    for line in eachline(filepath)
        line = strip(line)

        isempty(line) && continue
        startswith(line, "#") && continue

        if occursin("=", line)
            key, value = strip.(split(line, "=", limit=2))
            key = uppercase(key)
            value = uppercase(value)

            key == "AGEUNIT" && (ageunit = value; continue)
            key == "TUNIT" && (tunit = value; continue)
            key == "RHOUNIT" && (rhounit = value; continue)
            continue
        end

        parts = split(line)
        length(parts) < 3 && continue

        try
            tval = parse(Float64, parts[1])
            Tval = parse(Float64, parts[2])
            rhoval = parse(Float64, parts[3])

            if ageunit in ("YRS", "YR", "YEAR", "YEARS")
                tval *= 365.25 * 24.0 * 3600.0
            elseif !(ageunit in ("SEC", "S"))
                throw(DomainError(ageunit, "unsupported trajectory AGEUNIT"))
            end

            if tunit == "T8K"
                Tval /= 10.0
            elseif tunit != "T9K"
                throw(DomainError(tunit, "unsupported trajectory TUNIT; expected T9K or T8K"))
            end

            if rhounit == "LOG"
                rhoval = 10.0^rhoval
            elseif rhounit != "CGS"
                throw(DomainError(rhounit, "unsupported trajectory RHOUNIT; expected CGS or LOG"))
            end

            push!(time, tval)
            push!(temperature, Tval)
            push!(density, rhoval)
        catch
            continue
        end
    end

    return DataFrame(time_s=time, temperature_T9=temperature, density_cgs=density)
end
