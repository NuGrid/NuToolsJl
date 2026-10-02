#!/usr/bin/env python3
"""Convert a raw STARLIB release (as downloaded from the STARLIB github/site) into
the NPDATA/starlib reaction-data + sunet species files that starlib.F90 reads.

Raw STARLIB files have no leading "n reactions = X, n temp coordinates = Y" header
(starlib.F90's first read assumes one), so this script counts reactions and prepends
a header in the exact fixed-column format starlib.F90's '(13x,I5,23x,I4)' read expects.
It also builds a matching sunet (species Z/A) file, reusing Z/A from any existing
sunet files for species already known, and resolving new species itself.

Usage:
  python3 tools/build_starlib_files.py <raw_starlib_file> <tag>

Writes:
  NPDATA/starlib/starlib_<tag>.txt
  NPDATA/starlib/sunet_<tag>.dat
"""
import re
import sys
from pathlib import Path

NPDATA_STARLIB = Path(__file__).resolve().parents[3] / "physics" / "NPDATA" / "starlib"
EXISTING_SUNETS = ["sunet_taly_012025.dat", "sunet_mc10_mc13_082022.dat"]

# special-case species that don't follow the plain "<element><mass>" pattern
SPECIAL_SPECIES = {
    "n": (0, 1),
    "p": (1, 1),
    "d": (1, 2),
    "t": (1, 3),
    # additional excited/isomeric-state labels for Al-26 beyond the ground ("al-6")
    # and first isomer ("al*6") already present in the existing sunet files.
    "al01": (13, 26),
    "al02": (13, 26),
    "al03": (13, 26),
}

ELEMENT_Z = {
    sym.lower(): z
    for z, sym in enumerate(
        [
            "n", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si",
            "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
            "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc",
            "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
            "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
            "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
            "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
            "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg",
        ]
    )
}

SPECIES_RE = re.compile(r"^([a-z]+)(\d+)$")


def build_header(nrates: int, nt9: int) -> str:
    field1 = str(nrates)
    if len(field1) > 5:
        raise ValueError(f"nrates={nrates} needs more than 5 digits; format can't hold it")
    field2 = str(nt9)
    if len(field2) > 4:
        raise ValueError(f"nt9={nt9} needs more than 4 digits; format can't hold it")
    return "n reactions =" + field1.rjust(5) + ", n temp coordinates = " + field2.ljust(4) + "\n"


def extract_species(raw_path: Path):
    """Return (n_reactions, n_t9_coords, ordered unique species list)."""
    species_seen = []
    species_set = set()
    n_reactions = 0
    n_t9_coords = None
    with raw_path.open() as f:
        for line in f:
            try:
                float(line[:12])
                continue  # a T9/rate/uncertainty data row
            except ValueError:
                pass
            n_reactions += 1
            pos = 2 + 3
            row_coords = 0
            for _ in range(6):
                tok = line[pos:pos + 5].strip()
                pos += 5
                if tok and tok not in species_set:
                    species_set.add(tok)
                    species_seen.append(tok)
    # count t9 coords by re-scanning until the second reaction header
    with raw_path.open() as f:
        next(f)
        count = 0
        for line in f:
            try:
                float(line[:12])
                count += 1
            except ValueError:
                break
        n_t9_coords = count
    return n_reactions, n_t9_coords, species_seen


def load_known_species():
    known = {}
    for fn in EXISTING_SUNETS:
        path = NPDATA_STARLIB / fn
        if not path.is_file():
            continue
        with path.open() as f:
            next(f)
            for line in f:
                parts = line.split()
                if len(parts) >= 3:
                    known[parts[0]] = (int(parts[1]), int(parts[2]))
    return known


def resolve_species(species, known):
    resolved = {}
    unresolved = []
    for name in species:
        if name in known:
            resolved[name] = known[name]
            continue
        if name in SPECIAL_SPECIES:
            resolved[name] = SPECIAL_SPECIES[name]
            continue
        m = SPECIES_RE.match(name)
        if m:
            sym, mass = m.group(1), int(m.group(2))
            z = ELEMENT_Z.get(sym)
            if z is not None:
                resolved[name] = (z, mass)
                continue
        unresolved.append(name)
    return resolved, unresolved


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)
    raw_path = Path(sys.argv[1]).resolve()
    tag = sys.argv[2]

    n_reactions, n_t9_coords, species = extract_species(raw_path)
    print(f"reactions: {n_reactions}, t9 coords: {n_t9_coords}, unique species: {len(species)}")

    known = load_known_species()
    resolved, unresolved = resolve_species(species, known)
    if unresolved:
        print("ERROR: could not resolve Z/A for the following species:")
        for name in unresolved:
            print(" ", repr(name))
        print("Add them to SPECIAL_SPECIES or extend ELEMENT_Z, then re-run.")
        sys.exit(1)

    out_reaction = NPDATA_STARLIB / f"starlib_{tag}.txt"
    out_sunet = NPDATA_STARLIB / f"sunet_{tag}.dat"

    header = build_header(n_reactions, n_t9_coords)
    with raw_path.open() as src, out_reaction.open("w") as dst:
        dst.write(header)
        for line in src:
            dst.write(line)
    print(f"wrote {out_reaction}")

    with out_sunet.open("w") as f:
        f.write(f"{len(species)}\n")
        for name in species:
            z, a = resolved[name]
            f.write(f"{name:>5}  {z:>3}  {a:>3} \n")
    print(f"wrote {out_sunet}")


if __name__ == "__main__":
    main()
