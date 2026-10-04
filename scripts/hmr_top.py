#!/usr/bin/env python3
"""
Hydrogen mass repartitioning (HMR) written into the topology.

  hmr_top.py <topol.top> [--factor 3]

Every hydrogen in a solute molecule gets factor x its mass (1.008 -> 3.024 Da)
and the same amount is taken from the heavy atom it is bonded to, so each
molecule keeps its total mass. Together with constraints = h-bonds this allows
dt = 4 fs.

Why the topology and not the mdp: `mass-repartition-factor` only exists from
GROMACS 2024 on, and applying both would scale the hydrogens twice. Doing it
here works with any GROMACS version and puts the masses in a file you can read
and commit.

Rules
  - water (any moleculetype with [ settles ]) is left alone, as is usual for HMR
  - only files next to topol.top are edited (pdb2gmx chain .itp files,
    CHARMM-GUI toppar/); force-field files under GMXLIB are never touched
  - a file that already carries the marker, or any hydrogen heavier than
    2 Da (e.g. CHARMM-GUI's own HMR option), is not repartitioned again

No third-party dependencies.
"""
from __future__ import annotations

import argparse
import pathlib
import re
import sys

MARKER = "; HMR applied by scripts/hmr_top.py"
H_MAX = 1.1        # anything lighter is a hydrogen
HMR_MIN = 2.0      # anything this heavy is an already-repartitioned hydrogen


def includes(top: pathlib.Path, seen: set[pathlib.Path]) -> list[pathlib.Path]:
    """topol.top plus every #include that resolves to a local file."""
    out = []
    if top in seen or not top.exists():
        return out
    seen.add(top)
    out.append(top)
    for line in top.read_text().splitlines():
        m = re.match(r'\s*#include\s+"([^"]+)"', line)
        if m:
            inc = (top.parent / m.group(1)).resolve()
            if inc.exists():
                out.extend(includes(inc, seen))
    return out


def split_sections(lines: list[str]):
    """Yield (section_name, start_index) for every [ section ] header."""
    for i, line in enumerate(lines):
        m = re.match(r"\s*\[\s*(\w+)\s*\]", line)
        if m:
            yield m.group(1).lower(), i


def data_rows(lines, start, end):
    for i in range(start + 1, end):
        s = lines[i].split(";")[0].strip()
        if s and not s.startswith("#"):
            yield i, s.split()


def process(path: pathlib.Path, factor: float) -> tuple[int, list[str]]:
    lines = path.read_text().splitlines()
    if any(l.startswith(MARKER) for l in lines):
        return 0, [f"{path.name}: already repartitioned (marker found), skipped"]
    secs = list(split_sections(lines)) + [("__end__", len(lines))]
    # group sections by moleculetype
    mols, cur = [], None
    for (name, i), (_, nxt) in zip(secs, secs[1:]):
        if name == "moleculetype":
            cur = {"name": None, "atoms": None, "bonds": [], "settles": False, "hdr": (i, nxt)}
            mols.append(cur)
            for _, f in data_rows(lines, i, nxt):
                cur["name"] = f[0]
                break
        elif cur is not None and name == "atoms":
            cur["atoms"] = (i, nxt)
        elif cur is not None and name == "bonds":
            cur["bonds"].append((i, nxt))
        elif cur is not None and name == "settles":
            cur["settles"] = True
    changed, notes = 0, []
    for mol in mols:
        if mol["settles"] or mol["atoms"] is None:
            continue
        mass, row_of = {}, {}
        for i, f in data_rows(lines, *mol["atoms"]):
            if len(f) < 8:
                sys.exit(f"{path.name}: [ atoms ] row without a mass column in "
                         f"{mol['name']} -- cannot repartition:\n  {lines[i]}")
            nr = int(f[0])
            mass[nr], row_of[nr] = float(f[7]), i
        hyd = {n for n, m in mass.items() if m < H_MAX}
        if not hyd:
            continue
        if any(HMR_MIN <= m < 4.5 and n not in hyd for n, m in mass.items()):
            notes.append(f"{path.name}:{mol['name']}: hydrogens already heavy, skipped")
            continue
        partner = {}
        for sec in mol["bonds"]:
            for _, f in data_rows(lines, *sec):
                a, b = int(f[0]), int(f[1])
                if a in hyd and b not in hyd:
                    partner.setdefault(a, []).append(b)
                elif b in hyd and a not in hyd:
                    partner.setdefault(b, []).append(a)
        total0 = sum(mass.values())
        for h in sorted(hyd):
            heavy = partner.get(h, [])
            if len(heavy) != 1:
                sys.exit(f"{path.name}:{mol['name']}: hydrogen {h} is bonded to "
                         f"{len(heavy)} heavy atoms -- refusing to guess")
            dm = mass[h] * (factor - 1)
            mass[h] += dm
            mass[heavy[0]] -= dm
            if mass[heavy[0]] < mass[h]:
                sys.exit(f"{path.name}:{mol['name']}: atom {heavy[0]} would end up "
                         f"lighter than its hydrogens")
        for n, m in mass.items():
            i = row_of[n]
            code, _, comment = lines[i].partition(";")
            f = code.split()
            f[7] = f"{m:.5f}"
            row = (f"{f[0]:>6} {f[1]:>10} {f[2]:>6} {f[3]:>6} {f[4]:>6} {f[5]:>6} "
                   f"{f[6]:>10} {f[7]:>10}")
            lines[i] = " ".join([row] + f[8:]) + (f" ;{comment}" if comment else "")
        assert abs(sum(mass.values()) - total0) < 1e-3
        changed += len(hyd)
        notes.append(f"{path.name}:{mol['name']}: {len(hyd)} hydrogens repartitioned "
                     f"(total mass {total0:.2f} Da unchanged)")
    if changed:
        lines.insert(0, f"{MARKER} (factor {factor:g}); do not apply again")
        path.write_text("\n".join(lines) + "\n")
    return changed, notes


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("top")
    ap.add_argument("--factor", type=float, default=3.0)
    a = ap.parse_args()
    top = pathlib.Path(a.top).resolve()
    root = top.parent
    total = 0
    for f in includes(top, set()):
        if root not in f.parents and f != top:
            continue                      # never edit files outside the build dir
        n, notes = process(f, a.factor)
        total += n
        for x in notes:
            print(f"  {x}")
    if total == 0:
        print("  no hydrogens repartitioned (already done, or none found)")
    print(f"HMR: {total} hydrogens, factor {a.factor:g}")


if __name__ == "__main__":
    main()
