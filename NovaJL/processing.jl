struct Isotope
    Z::Int
    N::Int
    A::Int
    name::String
    X::Float64
end

struct AbundanceSet
    isotopes::Vector{Isotope}
end


function build_abundance_set(df)

    isotopes = Isotope[]

    for row in eachrow(df)

        parsed = parse_isotope(row.isotope)
        parsed === nothing && continue

        Z, N, A = parsed

        push!(isotopes,
            Isotope(Z, N, A, row.isotope, row.X)
        )
    end

    return AbundanceSet(isotopes)
end
