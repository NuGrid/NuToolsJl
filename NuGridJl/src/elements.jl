# elements.jl — element and isotope identity.
#
# nuppn text output is "stringly typed": a species is a right-justified token
# such as "HE  4", "PROT", "NEUT", or "AL 26", and every downstream consumer
# re-parses that string.  We lift those tokens once into a small value type so
# the rest of NuGridJl can dispatch on (Z, A, isomer) rather than on formatting.

"Element symbols indexed by proton number; `ELEMENT_SYMBOLS[Z + 1]`, with `Z = 0` the neutron."
const ELEMENT_SYMBOLS = (
    "n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",
)

"Neutron/proton shell closures, drawn as guide lines on the nuclear chart."
const MAGIC_NUMBERS = (2, 8, 20, 28, 50, 82, 126)

"Proton number for an element symbol (case-insensitive), or `nothing` if unknown."
function proton_number(symbol::AbstractString)
    isempty(symbol) && return nothing
    canonical = uppercasefirst(lowercase(symbol))
    idx = findfirst(==(canonical), ELEMENT_SYMBOLS)
    idx === nothing ? nothing : idx - 1
end

"Element symbol for proton number `Z`; falls back to `\"Z=…\"` outside the table."
function element_symbol(Z::Integer)
    0 <= Z < length(ELEMENT_SYMBOLS) ? ELEMENT_SYMBOLS[Z + 1] : string("Z=", Z)
end

"""
    Isotope(Z, A, isomer = 0)
    Isotope(token::AbstractString)

A nuclear species identified by proton number `Z`, mass number `A`, and isomeric
level `isomer` (0 = ground state, 1 = first metastable, …).

The string constructor accepts every spelling NuGrid tooling emits — `"HE  4"`,
`"he4"`, `"He-4"`, `"4He"`, mass-first REACLIB labels with an isomer flag such as
`"26Alg"`/`"26Alm"`, and the specials `"PROT"`/`"p"` and `"NEUT"`/`"n"`. Filler
tokens (`"OOOOO"`, blank) are *not* isotopes; use [`tryparse`](@ref) if you need
to detect them without an error.
"""
struct Isotope
    Z::Int
    A::Int
    isomer::Int
end

Isotope(Z::Integer, A::Integer) = Isotope(Int(Z), Int(A), 0)

_strip_separators(token::AbstractString) = replace(strip(token), r"[\s_\-]" => "")

# A bare "symbol + mass" core in either order, e.g. "He4"/"Mg24" or "4He"/"24Mg".
function _parse_core(core::AbstractString)
    m = match(r"^([A-Za-z]+)(\d+)$", core)
    if m !== nothing
        z = proton_number(m.captures[1])
        z === nothing || return (z, parse(Int, m.captures[2]))
    end
    m = match(r"^(\d+)([A-Za-z]+)$", core)
    if m !== nothing
        z = proton_number(m.captures[2])
        z === nothing || return (z, parse(Int, m.captures[1]))
    end
    return nothing
end

function Base.tryparse(::Type{Isotope}, token::AbstractString)
    t = strip(token)
    isempty(t) && return nothing
    up = uppercase(t)

    (t == "n" || up == "NEUT") && return Isotope(0, 1, 0)
    (t == "p" || up == "PROT") && return Isotope(1, 1, 0)
    all(==('O'), up) && return nothing              # "OOOOO" placeholder

    core = _strip_separators(t)
    isempty(core) && return nothing

    za = _parse_core(core)
    za !== nothing && return Isotope(za[1], za[2], 0)

    # Retry treating a trailing g/m as an isomer flag ("26Alm", "Al26m").
    last = lowercase(core[end])
    if last in ('g', 'm')
        za = _parse_core(core[1:prevind(core, lastindex(core))])
        za !== nothing && return Isotope(za[1], za[2], last == 'm' ? 1 : 0)
    end
    return nothing
end

function Isotope(token::AbstractString)
    iso = tryparse(Isotope, token)
    iso === nothing && throw(ArgumentError("could not parse isotope from \"$token\""))
    return iso
end

proton_number(iso::Isotope) = iso.Z
mass_number(iso::Isotope) = iso.A
neutron_number(iso::Isotope) = iso.A - iso.Z
element_symbol(iso::Isotope) = element_symbol(iso.Z)

"Human-readable name, e.g. `\"Al-26\"`, `\"Al-26m\"`, `\"n\"`, `\"p\"`."
function isotope_name(iso::Isotope)
    iso == Isotope(0, 1, 0) && return "n"
    iso == Isotope(1, 1, 0) && return "p"
    suffix = iso.isomer == 0 ? "" : (iso.isomer == 1 ? "m" : string("m", iso.isomer))
    string(element_symbol(iso.Z), "-", iso.A, suffix)
end

"LaTeX label for plots, e.g. `L\"{}^{26}\\mathrm{Al}\"`."
function latex_label(iso::Isotope)
    iso == Isotope(0, 1, 0) && return L"n"
    iso == Isotope(1, 1, 0) && return L"p"
    tag = iso.isomer == 0 ? "" : (iso.isomer == 1 ? "m" : string("m", iso.isomer))
    LaTeXString(string("{}^{", iso.A, tag, "}\\mathrm{", element_symbol(iso.Z), "}"))
end

Base.show(io::IO, iso::Isotope) = print(io, isotope_name(iso))

# --- stability -----------------------------------------------------------------
# Stable mass numbers per element (Z => sorted A's); ported from the NuGrid
# reference table.  Elements with no stable isotope are simply absent.
const _STABLE_TABLE = (
    (1, (1, 2)), (2, (3, 4)), (3, (6, 7)), (4, (9,)), (5, (10, 11)),
    (6, (12, 13)), (7, (14, 15)), (8, (16, 17, 18)), (9, (19,)),
    (10, (20, 21, 22)), (11, (23,)), (12, (24, 25, 26)), (13, (27,)),
    (14, (28, 29, 30)), (15, (31,)), (16, (32, 33, 34, 36)), (17, (35, 37)),
    (18, (36, 38, 40)), (19, (39, 40, 41)), (20, (40, 42, 43, 44, 46, 48)),
    (21, (45,)), (22, (46, 47, 48, 49, 50)), (23, (50, 51)),
    (24, (50, 52, 53, 54)), (25, (55,)), (26, (54, 56, 57, 58)), (27, (59,)),
    (28, (58, 60, 61, 62, 64)), (29, (63, 65)), (30, (64, 66, 67, 68, 70)),
    (31, (69, 71)), (32, (70, 72, 73, 74, 76)), (33, (75,)),
    (34, (74, 76, 77, 78, 80, 82)), (35, (79, 81)),
    (36, (78, 80, 82, 83, 84, 86)), (37, (85, 87)), (38, (84, 86, 87, 88)),
    (39, (89,)), (40, (90, 91, 92, 94, 96)), (41, (93,)),
    (42, (92, 94, 95, 96, 97, 98, 100)), (44, (96, 98, 99, 100, 101, 102, 104)),
    (45, (103,)), (46, (102, 104, 105, 106, 108, 110)), (47, (107, 109)),
    (48, (106, 108, 110, 111, 112, 113, 114, 116)), (49, (113, 115)),
    (50, (112, 114, 115, 116, 117, 118, 119, 120, 122, 124)), (51, (121, 123)),
    (52, (120, 122, 123, 124, 125, 126, 128, 130)), (53, (127,)),
    (54, (124, 126, 128, 129, 130, 131, 132, 134, 136)), (55, (133,)),
    (56, (130, 132, 134, 135, 136, 137, 138)), (57, (138, 139)),
    (58, (136, 138, 140, 142)), (59, (141,)),
    (60, (142, 143, 144, 145, 146, 148, 150)),
    (62, (144, 147, 148, 149, 150, 152, 154)), (63, (151, 153)),
    (64, (152, 154, 155, 156, 157, 158, 160)), (65, (159,)),
    (66, (156, 158, 160, 161, 162, 163, 164)), (67, (165,)),
    (68, (162, 164, 166, 167, 168, 170)), (69, (169,)),
    (70, (168, 170, 171, 172, 173, 174, 176)), (71, (175, 176)),
    (72, (174, 176, 177, 178, 179, 180)), (73, (180, 181)),
    (74, (180, 182, 183, 184, 186)), (75, (185, 187)),
    (76, (184, 186, 187, 188, 189, 190, 192)), (77, (191, 193)),
    (78, (190, 192, 194, 195, 196, 198)), (79, (197,)),
    (80, (196, 198, 199, 200, 201, 202, 204)), (81, (203, 205)),
    (82, (204, 206, 207, 208)), (83, (209,)), (90, (232,)), (92, (235, 238)),
)

const STABLE_ISOTOPES = Dict{Int,Vector{Int}}(z => collect(as) for (z, as) in _STABLE_TABLE)

"Whether `iso` is a ground-state stable isotope (isomers are never stable)."
is_stable(iso::Isotope) =
    iso.isomer == 0 && haskey(STABLE_ISOTOPES, iso.Z) && iso.A in STABLE_ISOTOPES[iso.Z]

"Is proton or neutron number `n` a nuclear magic number?"
is_magic(n::Integer) = n in MAGIC_NUMBERS
