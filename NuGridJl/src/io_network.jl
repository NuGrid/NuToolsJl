# io_network.jl — the reaction network and the ppn input namelists.
#
# `networksetup.txt` carries two tables: the list of tracked isotopes and, below
# it, ~14k reaction rows.  Each reaction row records an active flag, the reacting
# species, the printed rate, a source label (VITAL, NACRR, STL01, …), a reaction
# type (p,g)/(a,p)/…, the REACLIB chapter, the applied rate multiplier and the
# Q-value — everything a sensitivity study needs to know which rate is in play.

"""
    Reaction

One network reaction: `index`, `active` flag, expanded `reactants`/`products`
(stoichiometric multiplicity preserved, filler species dropped), the `rtype`
code such as `\"(p,g)\"`, the rate `source` label, the printed `rate`, the applied
`multiplier`, the `qvalue`, and the REACLIB `chapter`.
"""
struct Reaction
    index::Int
    active::Bool
    reactants::Vector{Isotope}
    products::Vector{Isotope}
    rtype::String
    source::String
    rate::Float64
    multiplier::Float64
    qvalue::Float64
    chapter::Int
end

"""
    Network

Everything in `networksetup.txt`: the tracked `isotopes` and the `reactions`.
"""
struct Network
    isotopes::Vector{Isotope}
    reactions::Vector{Reaction}
end

# Rate/multiplier/Q columns are E-notation here, but tolerate D and the bare
# "1.234+18" (missing E) Fortran variants.
function _parse_network_num(value)
    text = strip(String(value))
    text = replace(text, r"[Dd]" => "E")
    if !occursin(r"[Ee]", text)
        text = replace(text, r"^([+-]?(?:\d+\.?\d*|\.\d+))([+-]\d+)$" => s"\1E\2")
    end
    return parse(Float64, text)
end

_expand(count, iso) = iso === nothing ? Isotope[] : fill(iso, max(count, 1))

const _NETWORK_RE = r"^\s*(\d+)\s+([TF])\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+->\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)"

function _parse_reaction(line)
    m = match(_NETWORK_RE, line)
    m === nothing && return nothing
    reactants = vcat(
        _expand(parse(Int, m.captures[3]), tryparse(Isotope, m.captures[4])),
        _expand(parse(Int, m.captures[5]), tryparse(Isotope, m.captures[6])),
    )
    products = vcat(
        _expand(parse(Int, m.captures[7]), tryparse(Isotope, m.captures[8])),
        _expand(parse(Int, m.captures[9]), tryparse(Isotope, m.captures[10])),
    )
    return Reaction(
        parse(Int, m.captures[1]),
        m.captures[2] == "T",
        reactants, products,
        String(m.captures[13]),
        String(m.captures[12]),
        _parse_network_num(m.captures[11]),
        _parse_network_num(m.captures[15]),
        _parse_network_num(m.captures[16]),
        parse(Int, m.captures[14]),
    )
end

# Isotope-table rows look like "   5   HE  4     T    4.     2.   1"; the symbol
# spans a variable number of tokens (and can itself be "F", fluorine), so key off
# the trailing numeric columns MASSNR / PROTON / ISOMER instead of the T/F flag.
function _parse_isotope_row(line)
    p = split(line)
    length(p) < 4 && return nothing
    tryparse(Int, p[1]) === nothing && return nothing
    a = tryparse(Float64, p[end - 2]); a === nothing && return nothing
    z = tryparse(Float64, p[end - 1]); z === nothing && return nothing
    isom = tryparse(Float64, p[end]); isom === nothing && return nothing
    return Isotope(Int(round(z)), Int(round(a)), max(Int(round(isom)) - 1, 0))
end

"""
    read_network(path) -> Network

Parse `networksetup.txt` into its isotope list and `Reaction` vector.
"""
function read_network(path::AbstractString)
    isfile(path) || throw(ArgumentError("networksetup file not found: $path"))
    isotopes = Isotope[]
    reactions = Reaction[]
    in_isotopes = false
    for raw in eachline(path)
        if occursin("TABLE OF ISOTOPES", raw)
            in_isotopes = true; continue
        elseif occursin("REACTION NETWORK", raw)
            in_isotopes = false; continue
        end
        if in_isotopes
            iso = _parse_isotope_row(raw)
            iso === nothing || push!(isotopes, iso)
        else
            r = _parse_reaction(raw)
            r === nothing || push!(reactions, r)
        end
    end
    return Network(isotopes, reactions)
end

# --- reaction labels -----------------------------------------------------------

"Nuclear-physics species name, mass first: `\"26Al\"`, `\"p\"`, `\"n\"`, `\"g\"`."
function reaction_species_name(iso::Isotope)
    iso == Isotope(1, 1, 0) && return "p"
    iso == Isotope(0, 1, 0) && return "n"
    tag = iso.isomer == 0 ? "" : (iso.isomer == 1 ? "m" : string("m", iso.isomer))
    string(iso.A, element_symbol(iso.Z), tag)
end

_heaviest(isos) = isempty(isos) ? nothing : isos[argmax(mass_number.(isos))]
_rtype_inner(rtype) = strip(rtype, ['(', ')', ' '])

"""
    label(r::Reaction) -> String

Compact reaction label, e.g. `\"25Mg(p,g)26Al\"`: heaviest reactant, the reaction
code, heaviest product.
"""
function label(r::Reaction)
    heavy_in = _heaviest(r.reactants)
    heavy_out = _heaviest(r.products)
    (heavy_in === nothing || heavy_out === nothing) && return string(r.rtype)
    string(reaction_species_name(heavy_in), "(", _rtype_inner(r.rtype), ")",
           reaction_species_name(heavy_out))
end

const _GREEK = Dict('g' => "\\gamma", 'a' => "\\alpha", 'v' => "\\nu", 'b' => "\\beta")

_greekify(code) = join(get(_GREEK, lowercase(c), string(c)) for c in code)

"""
    latex_label(r::Reaction) -> LaTeXString

`label` typeset for plots, with superscript masses and γ/α/ν substitutions:
`25Mg(p,g)26Al` → `{}^{25}\\mathrm{Mg}(p,\\gamma){}^{26}\\mathrm{Mg}`-style.
"""
function latex_label(r::Reaction)
    heavy_in = _heaviest(r.reactants)
    heavy_out = _heaviest(r.products)
    (heavy_in === nothing || heavy_out === nothing) && return LaTeXString(r.rtype)
    inner = join((_greekify(strip(s)) for s in split(_rtype_inner(r.rtype), ',')), ",")
    LaTeXString(string(
        _latex_species(heavy_in), "(", inner, ")", _latex_species(heavy_out),
    ))
end

function _latex_species(iso::Isotope)
    iso == Isotope(1, 1, 0) && return "p"
    iso == Isotope(0, 1, 0) && return "n"
    tag = iso.isomer == 0 ? "" : "m"
    string("{}^{", iso.A, tag, "}\\mathrm{", element_symbol(iso.Z), "}")
end

# --- ppn input namelists -------------------------------------------------------

"Coerce a Fortran namelist value string to Bool / Int / Float64 / stripped String."
function _coerce_namelist(value)
    v = strip(value)
    up = uppercase(v)
    (up == "T" || up == ".TRUE.") && return true
    (up == "F" || up == ".FALSE.") && return false
    if (startswith(v, '\'') && endswith(v, '\'')) || (startswith(v, '"') && endswith(v, '"'))
        return String(v[2:end-1])
    end
    n = tryparse(Int, v)
    n !== nothing && return n
    f = tryparse(Float64, replace(v, r"[Dd]" => "E"))
    f !== nothing && return f
    return String(v)
end

"""
    read_namelist(path) -> Dict{Symbol,Any}

Parse a ppn `&group … /` Fortran namelist (`ppn_frame/physics/solver.input`) into
a dictionary, stripping inline `!` comments and coercing values to native types.
"""
function read_namelist(path::AbstractString)
    out = Dict{Symbol,Any}()
    for raw in eachline(path)
        line = strip(raw)
        (isempty(line) || startswith(line, "!") || startswith(line, "&") ||
            startswith(line, "/") || startswith(line, "#")) && continue
        occursin('=', line) || continue
        code = split(line, '!'; limit = 2)[1]
        occursin('=', code) || continue
        key, value = strip.(split(code, '='; limit = 2))
        isempty(key) && continue
        out[Symbol(key)] = _coerce_namelist(value)
    end
    return out
end
