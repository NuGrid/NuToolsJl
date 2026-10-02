#!/usr/bin/env python3
"""
list_network_rates.py

Output the list of active nuclear reactions and their rate sources for a given
nova case. Intended for sharing with collaborators (e.g. hydrodynamics groups)
so they know exactly which rate database was used for each reaction.

Usage:
  python3 tools/list_network_rates.py --nova ne_nova_1.15_12_X_weiss_mixed
  python3 tools/list_network_rates.py --nova ne_nova_1.15_12_X_weiss_mixed --run runs_star_opt2_novital
  python3 tools/list_network_rates.py --nova ne_nova_1.15_12_X_weiss_mixed --run runs --out my_rates.csv

Output:
  <nova_case>/results/network_rate_list_<run>.csv   (CSV)
  <nova_case>/results/network_rate_list_<run>.txt   (human-readable table)
"""

import re
import sys
import os
import csv
import argparse
from pathlib import Path


# ── Source label descriptions ────────────────────────────────────────────────

SOURCE_NAMES = {
    "STL02": "STARLIB (Sallaska et al. 2013, ETR25 Monte Carlo rates)",
    "NACRR": "NACRE II recommended (Xu et al. 2013)",
    "NACRL": "NACRE II lower limit",
    "NACRU": "NACRE II upper limit",
    "ILI01": "Iliadis et al. 2001 (shell-model, charged particle)",
    "JINAC": "JINA REACLIB v2 (Cyburt et al. 2010)",
    "JINAR": "JINA REACLIB v2 reverse",
    "JINAV": "JINA REACLIB v2 (VITAL supplement)",
    "RVRSE": "Reverse rate via detailed balance",
    "VITAL": "VITAL internal rate (vital.F90; pp-chain / compound)",
    "NETB1": "Nuclear structure (beta decays, Nubase/ENSDF)",
    "ODA94": "Oda et al. 1994 (weak rates)",
    "LMP00": "Langanke & Martínez-Pinedo 2000 (stellar weak rates)",
    "KADON": "KADoNiS neutron capture database",
    "BASEL": "Basel (Hauser-Feshbach, Rauscher & Thielemann 2000)",
    "KADR ": "KADoNiS recommended",
    "CF88 ": "Caughlan & Fowler 1988",
    "FKTH ": "Fowler, Caughlan & Zimmermann 1975",
    "FFW85": "Fuller, Fowler & Newman 1985 (stellar weak rates)",
}

# ── Reaction-type channel notation ───────────────────────────────────────────

RTYPE_CHANNEL = {
    "(p,g)": ("p",  "γ"),
    "(p,a)": ("p",  "α"),
    "(p,n)": ("p",  "n"),
    "(a,g)": ("α",  "γ"),
    "(a,p)": ("α",  "p"),
    "(a,n)": ("α",  "n"),
    "(n,g)": ("n",  "γ"),
    "(n,p)": ("n",  "p"),
    "(n,a)": ("n",  "α"),
    "(g,p)": ("γ",  "p"),
    "(g,a)": ("γ",  "α"),
    "(g,n)": ("γ",  "n"),
    "(+,g)": ("β+", "ν"),
    "(-,g)": ("β-", "ν̄"),
    "(e,v)": ("e-", "ν"),
    "(v,v)": None,          # compound / multi-body: handled separately
}

# ── Species formatting ────────────────────────────────────────────────────────

_SP_RE = re.compile(r'^([A-Z]+)\s+(\d+)$')

def parse_species(sp: str):
    """Return (symbol, mass_number) or None for OOOOO."""
    sp = sp.strip()
    if sp in ("OOOOO", ""):
        return None
    if sp == "NEUT":
        return ("n", 0)
    if sp == "PROT":
        return ("H", 1)
    m = _SP_RE.match(sp)
    if m:
        raw, mass = m.group(1), int(m.group(2))
        sym = raw[0] + raw[1:].lower() if len(raw) > 1 else raw
        return (sym, mass)
    return (sp, 0)


def fmt_sp(sp: str) -> str:
    """Format species as '12C', 'p', 'n', '4He'. Returns '' for OOOOO."""
    parsed = parse_species(sp)
    if parsed is None:
        return ""
    sym, mass = parsed
    if sym == "n" and mass == 0:
        return "n"
    if mass == 1 and sym == "H":
        return "p"
    return f"{mass}{sym}"


def fmt_sp_plain(sp: str) -> str:
    """Same but replaces Greek with ASCII: α→a, γ→g, β→b, ν→nu."""
    s = fmt_sp(sp)
    return s  # No Greek in species names; only in channel labels


# ── Reaction notation ─────────────────────────────────────────────────────────

def build_notation(n1, sp1, n2, sp2, n3, sp3, n4, sp4, rtype):
    """Build standard A(b,c)D notation or a fallback arrow notation."""
    channel = RTYPE_CHANNEL.get(rtype)

    if channel is None:
        # compound / (v,v) or unknown: arrow notation
        lhs_parts = []
        if int(n1) > 0 and sp1.strip() != "OOOOO":
            prefix = f"{n1}×" if int(n1) > 1 else ""
            lhs_parts.append(prefix + fmt_sp(sp1))
        if int(n2) > 0 and sp2.strip() != "OOOOO":
            prefix = f"{n2}×" if int(n2) > 1 else ""
            lhs_parts.append(prefix + fmt_sp(sp2))
        rhs_parts = []
        if int(n3) > 0 and sp3.strip() != "OOOOO":
            prefix = f"{n3}×" if int(n3) > 1 else ""
            rhs_parts.append(prefix + fmt_sp(sp3))
        if int(n4) > 0 and sp4.strip() != "OOOOO":
            prefix = f"{n4}×" if int(n4) > 1 else ""
            rhs_parts.append(prefix + fmt_sp(sp4))
        lhs = " + ".join(lhs_parts) if lhs_parts else "?"
        rhs = " + ".join(rhs_parts) if rhs_parts else "?"
        return f"{lhs} → {rhs}"

    proj, eject = channel

    # For photodisintegration (γ,x): target is sp1, products are sp3+sp4
    if proj == "γ":
        target = fmt_sp(sp1)
        products = [s for s in (fmt_sp(sp3), fmt_sp(sp4)) if s]
        residual = products[0] if products else "?"
        return f"{target}(γ,{eject}){residual}"

    # For decays (+,g), (-,g), (e,v): sp1 is parent, sp3 is daughter
    if proj in ("β+", "β-", "e-"):
        parent = fmt_sp(sp1)
        daughter = fmt_sp(sp3)
        return f"{parent}({proj},{eject}){daughter}"

    # Standard 2-body: target is sp1, projectile matches proj label, product is sp3
    target = fmt_sp(sp1)
    residual = fmt_sp(sp3)
    if not residual:
        residual = fmt_sp(sp4)
    return f"{target}({proj},{eject}){residual}"


# ── Number parsing ───────────────────────────────────────────────────────────

_FORT_FLOAT_RE = re.compile(r'^([+-]?\d+\.?\d*)([+-]\d+)$')

def parse_float(s: str) -> float:
    """Parse a Fortran float that may have no 'E' before a 3-digit exponent."""
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        m = _FORT_FLOAT_RE.match(s)
        if m:
            return float(f"{m.group(1)}E{m.group(2)}")
        raise ValueError(f"Cannot parse float: {s!r}")


# ── networksetup.txt parsing ──────────────────────────────────────────────────

RXN_RE = re.compile(
    r'^\s*(\d+)\s+([TF])\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+->\s+(\d+)\s+(.{5})\s+\+\s+(\d+)\s+(.{5})\s+'
    r'(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\S+)\s+(\S+)'
)


def parse_networksetup(path: Path) -> list[dict]:
    rows = []
    for line in path.read_text().splitlines():
        m = RXN_RE.match(line)
        if not m:
            continue
        (idx, act, n1, sp1, n2, sp2, n3, sp3, n4, sp4,
         q_val, src, rtype, ilabb, rfac, bind_diff) = m.groups()
        rows.append({
            "index":    int(idx),
            "active":   act == "T",
            "n1": n1, "sp1": sp1.strip(),
            "n2": n2, "sp2": sp2.strip(),
            "n3": n3, "sp3": sp3.strip(),
            "n4": n4, "sp4": sp4.strip(),
            "Q_MeV":    parse_float(q_val),
            "source":   src,
            "rtype":    rtype,
        })
    return rows


# ── Output formatting ─────────────────────────────────────────────────────────

def write_csv(rows: list[dict], path: Path):
    fieldnames = ["index", "reaction", "rtype", "Q_MeV", "source", "source_description"]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({
                "index":              r["index"],
                "reaction":           r["notation"],
                "rtype":              r["rtype"],
                "Q_MeV":             f"{r['Q_MeV']:.4f}",
                "source":             r["source"],
                "source_description": SOURCE_NAMES.get(r["source"], r["source"]),
            })


def write_text(rows: list[dict], path: Path, nova: str, run_dir: str):
    src_counts = {}
    for r in rows:
        src_counts[r["source"]] = src_counts.get(r["source"], 0) + 1

    with path.open("w") as f:
        f.write(f"Network rate source list\n")
        f.write(f"Nova case : {nova}\n")
        f.write(f"Run dir   : {run_dir}\n")
        f.write(f"Reactions : {len(rows)} active\n")
        f.write("=" * 90 + "\n\n")

        f.write(f"{'#':>5}  {'Reaction':<38}  {'Type':<7}  {'Q (MeV)':>9}  {'Source'}\n")
        f.write("-" * 90 + "\n")
        for r in rows:
            f.write(
                f"{r['index']:>5}  {r['notation']:<38}  {r['rtype']:<7}  "
                f"{r['Q_MeV']:>9.4f}  {r['source']}\n"
            )

        f.write("\n" + "=" * 90 + "\n")
        f.write("Rate source legend\n")
        f.write("-" * 90 + "\n")
        for src in sorted(src_counts, key=lambda s: -src_counts[s]):
            desc = SOURCE_NAMES.get(src, "(unknown)")
            f.write(f"  {src:<6}  ({src_counts[src]:4d} reactions)  {desc}\n")


# ── Entry point ───────────────────────────────────────────────────────────────

def nova_cases_dir() -> Path:
    here = Path(__file__).resolve().parent
    return here.parent / "nova_cases"


def main():
    parser = argparse.ArgumentParser(
        description="List active reactions and rate sources for a nova case."
    )
    parser.add_argument("--nova",    required=True,
                        help="Nova case directory name under nova_cases/")
    parser.add_argument("--run",     default="runs",
                        help="Run directory to read (default: runs)")
    parser.add_argument("--out",     default=None,
                        help="Output file stem (no extension). "
                             "Default: <nova_case>/results/network_rate_list_<run>")
    parser.add_argument("--active-only", action="store_true", default=True,
                        help="Include only T-flagged reactions (default: True)")
    parser.add_argument("--include-vital", action="store_true", default=False,
                        help="Include VITAL block reactions (index <= 117, default: omit)")
    args = parser.parse_args()

    nova_dir = nova_cases_dir() / args.nova
    if not nova_dir.is_dir():
        sys.exit(f"ERROR: nova case not found: {nova_dir}")

    nwsetup = nova_dir / args.run / "baseline" / "networksetup.txt"
    if not nwsetup.is_file():
        sys.exit(f"ERROR: networksetup.txt not found: {nwsetup}")

    print(f"Parsing {nwsetup} ...")
    all_rows = parse_networksetup(nwsetup)

    rows = [r for r in all_rows if r["active"]]
    if not args.include_vital:
        rows = [r for r in rows if r["index"] > 117]

    for r in rows:
        r["notation"] = build_notation(
            r["n1"], r["sp1"], r["n2"], r["sp2"],
            r["n3"], r["sp3"], r["n4"], r["sp4"],
            r["rtype"],
        )

    # Output paths
    results_dir = nova_dir / "results"
    results_dir.mkdir(exist_ok=True)
    stem = args.out or str(results_dir / f"network_rate_list_{args.run}")
    csv_path  = Path(stem).with_suffix(".csv")
    txt_path  = Path(stem).with_suffix(".txt")

    write_csv(rows, csv_path)
    write_text(rows, txt_path, args.nova, args.run)

    print(f"  {len(rows)} active reactions written")
    print(f"  CSV : {csv_path}")
    print(f"  Text: {txt_path}")

    # Brief source summary
    src_counts = {}
    for r in rows:
        src_counts[r["source"]] = src_counts.get(r["source"], 0) + 1
    print("\nSource breakdown:")
    for src, n in sorted(src_counts.items(), key=lambda x: -x[1]):
        print(f"  {src:<6}  {n:4d}  {SOURCE_NAMES.get(src, '')}")


if __name__ == "__main__":
    main()
