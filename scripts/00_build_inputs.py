#!/usr/bin/env python3
"""
Build the starting structure of one system: input/<system>/<system>.pdb

This is the step the rest of the pipeline cannot check for you. It does the
mechanical part reproducibly and prints everything you need to judge the
scientific part (which dimer, which RNA, which gaps).

Subcommands
-----------
  fetch <PDBID>          download input/_pdb/<PDBID>.pdb from RCSB
  show  <PDBID|file>     chains, residue ranges, sequences, gaps, BIOMT operators
  build <system>         assemble the chains listed in config/systems.tsv,
                         apply RNA mutations from config/rna_designs.tsv, clean,
                         and write input/<system>/<system>.pdb

What `build` does, in order
  1. reads the source named in the registry's source_pdb column:
       <PDBID>          input/_pdb/<PDBID>.pdb (fetched if missing)
       model            input/<system>/model.pdb (AlphaFold / grafted model)
       from:<system>    input/<other>/<other>.pdb, already built
       graft:<system>   single-chain model input/<system>/model.pdb superposed on
                        the protein of an already built complex, whose RNA and
                        Mg2+ are taken over (scCP + pac, S9)
  2. first MODEL only; one alternate location per residue; MSE -> MET;
     hydrogens, waters and every other HETATM dropped (Mg2+ kept for RNA systems)
  3. assembles the `chains` column. The first chain is taken as deposited; each
     further chain is placed with the BIOMT operator (REMARK 350) that gives it
     the most contacts with what is already assembled. That is how a C/C dimer
     is recovered from a capsid asymmetric unit, where the second C is a
     symmetry mate that is not in the file.
  4. mutates the RNA in place: shared sugar-phosphate and ring atoms are kept,
     atoms that differ are removed so CHARMM-GUI (or pdbfixer) rebuilds them
  5. reports chain breaks and missing heavy atoms, then writes the PDB

No third-party dependencies.
"""
from __future__ import annotations

import argparse
import math
import pathlib
import sys
import urllib.request

REPO = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))
import _registry  # noqa: E402

INPUT = REPO / "input"
PDB_CACHE = INPUT / "_pdb"
DESIGNS = REPO / "config" / "rna_designs.tsv"

# --------------------------------------------------------------------------
# residue vocabulary
# --------------------------------------------------------------------------

AA_HEAVY = {
    "ALA": "N CA C O CB",
    "ARG": "N CA C O CB CG CD NE CZ NH1 NH2",
    "ASN": "N CA C O CB CG OD1 ND2",
    "ASP": "N CA C O CB CG OD1 OD2",
    "CYS": "N CA C O CB SG",
    "GLN": "N CA C O CB CG CD OE1 NE2",
    "GLU": "N CA C O CB CG CD OE1 OE2",
    "GLY": "N CA C O",
    "HIS": "N CA C O CB CG ND1 CD2 CE1 NE2",
    "ILE": "N CA C O CB CG1 CG2 CD1",
    "LEU": "N CA C O CB CG CD1 CD2",
    "LYS": "N CA C O CB CG CD CE NZ",
    "MET": "N CA C O CB CG SD CE",
    "PHE": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ",
    "PRO": "N CA C O CB CG CD",
    "SER": "N CA C O CB OG",
    "THR": "N CA C O CB OG1 CG2",
    "TRP": "N CA C O CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2",
    "TYR": "N CA C O CB CG CD1 CD2 CE1 CE2 CZ OH",
    "VAL": "N CA C O CB CG1 CG2",
}
AA_HEAVY = {k: set(v.split()) for k, v in AA_HEAVY.items()}
for _h in ("HSD", "HSE", "HSP", "HID", "HIE", "HIP"):
    AA_HEAVY[_h] = AA_HEAVY["HIS"]
AA_ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
          "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
          "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
          "TYR": "Y", "VAL": "V"}

# RNA residue names seen in the wild -> one letter
RNA_NAMES = {}
for _b, _alts in {"A": "A RA RA5 RA3 ADE", "C": "C RC RC5 RC3 CYT",
                  "G": "G RG RG5 RG3 GUA", "U": "U RU RU5 RU3 URA"}.items():
    for _n in _alts.split():
        RNA_NAMES[_n] = _b

SUGAR_PHOSPHATE = {"P", "OP1", "OP2", "OP3", "O1P", "O2P", "O3P", "O5'", "C5'",
                   "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}
PURINE_RING = ["N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"]
PYRIMIDINE_RING = ["N1", "C2", "O2", "N3", "C4", "C5", "C6"]
BASE_HEAVY = {
    "A": set(PURINE_RING) | {"N6"},
    "G": set(PURINE_RING) | {"O6", "N2"},
    "C": set(PYRIMIDINE_RING) | {"N4"},
    "U": set(PYRIMIDINE_RING) | {"O4"},
}
PURINES = {"A", "G"}
# glycosidic frame kept across a purine <-> pyrimidine swap, so the new base
# is rebuilt in the same anti orientation as the old one
PUR_TO_PYR = {"N9": "N1", "C4": "C2", "C8": "C6"}
PYR_TO_PUR = {v: k for k, v in PUR_TO_PYR.items()}

# Modified nucleotides seen in MS2 soaking structures -> parent base, and the
# extra atoms to strip (2BU1 carries 5-bromouracil at -5).
MODIFIED_NT = {"5BU": ("U", {"BR"})}

WATERS = {"HOH", "WAT", "H2O", "DOD", "TIP", "TIP3", "SOL"}
MG_NAMES = {"MG", "MG2"}


def kind(resname: str) -> str:
    if resname in AA_HEAVY or resname == "MSE":
        return "protein"
    if resname in RNA_NAMES or resname in MODIFIED_NT:
        return "rna"
    if resname in WATERS:
        return "water"
    if resname in MG_NAMES:
        return "mg"
    return "het"


# --------------------------------------------------------------------------
# PDB I/O
# --------------------------------------------------------------------------

class Atom:
    __slots__ = ("name", "resname", "chain", "resseq", "icode", "x", "y", "z",
                 "occ", "b", "element")

    def __init__(self, line: str):
        self.name = line[12:16].strip()
        self.resname = line[17:20].strip()
        self.chain = line[21]
        self.resseq = int(line[22:26])
        self.icode = line[26]
        self.x, self.y, self.z = (float(line[30:38]), float(line[38:46]),
                                  float(line[46:54]))
        self.occ = float(line[54:60] or 1.0) if line[54:60].strip() else 1.0
        self.b = float(line[60:66]) if line[60:66].strip() else 0.0
        el = line[76:78].strip() if len(line) >= 78 else ""
        self.element = el or "".join(c for c in self.name if c.isalpha())[:1]
        self.name = self.name.replace("*", "'")   # pre-remediation sugar names

    @property
    def res(self) -> tuple:
        return (self.chain, self.resseq, self.icode)

    def xyz(self) -> tuple:
        return (self.x, self.y, self.z)

    def copy(self) -> "Atom":
        a = Atom.__new__(Atom)
        for s in Atom.__slots__:
            setattr(a, s, getattr(self, s))
        return a


def read_pdb(path: pathlib.Path) -> tuple[list[Atom], list[tuple], list[str]]:
    """Atoms (first model, one altloc), BIOMT operators, and notes."""
    atoms: list[Atom] = []
    biomt: dict[int, list] = {}
    in_bm1 = False
    chosen_alt: dict[tuple, str] = {}
    notes: list[str] = []
    for line in path.read_text(errors="ignore").splitlines():
        rec = line[:6]
        if rec == "ENDMDL":
            notes.append("multiple MODELs: only the first was read")
            break
        if line.startswith("REMARK 350"):
            if "BIOMOLECULE:" in line:
                in_bm1 = line.split("BIOMOLECULE:")[1].strip() == "1"
            elif "BIOMT" in line and in_bm1:
                f = line.split()
                row, op = int(f[2][-1]), int(f[3])
                biomt.setdefault(op, [None, None, None])[row - 1] = \
                    [float(v) for v in f[4:8]]
            continue
        if rec not in ("ATOM  ", "HETATM"):
            continue
        alt = line[16]
        a = Atom(line)
        if alt != " ":
            key = a.res
            chosen_alt.setdefault(key, alt)
            if alt != chosen_alt[key]:
                continue
        atoms.append(a)
    if chosen_alt:
        notes.append(f"{len(chosen_alt)} residue(s) with alternate locations: "
                     f"kept the first altloc of each")
    ops = [tuple(tuple(r) for r in m) for _, m in sorted(biomt.items())
           if all(r is not None for r in m)]
    return atoms, ops, notes


def resolve_source(src: str, system: str | None = None) -> pathlib.Path:
    if src == "model":
        p = INPUT / system / "model.pdb"
        if not p.exists():
            sys.exit(f"source is 'model': put your model at {p.relative_to(REPO)}")
        return p
    if src.startswith("from:"):
        other = src[5:]
        p = INPUT / other / f"{other}.pdb"
        if not p.exists():
            sys.exit(f"build {other} first ({p.relative_to(REPO)} missing)")
        return p
    if src.lower() in ("cocrystal", "todo", "-", ""):
        sys.exit(f"source_pdb for {system} is '{src}': replace it in "
                 f"config/systems.tsv with the actual PDB ID")
    p = pathlib.Path(src)
    if p.suffix in (".pdb", ".ent") and p.exists():
        return p
    return fetch(src)


def fetch(pdbid: str) -> pathlib.Path:
    pdbid = pdbid.upper()
    PDB_CACHE.mkdir(parents=True, exist_ok=True)
    p = PDB_CACHE / f"{pdbid}.pdb"
    if p.exists():
        return p
    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    print(f"fetching {url}", file=sys.stderr)
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            p.write_bytes(r.read())
    except Exception as e:  # noqa: BLE001
        sys.exit(f"could not download {pdbid} ({e}). Download it by hand into "
                 f"{p.relative_to(REPO)}")
    return p


def fmt_atom(serial: int, a: Atom, record: str = "ATOM  ") -> str:
    name = a.name if len(a.name) == 4 or len(a.element) == 2 else f" {a.name}"
    return (f"{record}{serial % 100000:5d} {name:<4} {a.resname:>3} {a.chain}"
            f"{a.resseq:4d}{a.icode}   {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
            f"{a.occ:6.2f}{a.b:6.2f}          {a.element:>2}")


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def apply(op: tuple, a: Atom) -> Atom:
    n = a.copy()
    x, y, z = a.xyz()
    n.x = op[0][0] * x + op[0][1] * y + op[0][2] * z + op[0][3]
    n.y = op[1][0] * x + op[1][1] * y + op[1][2] * z + op[1][3]
    n.z = op[2][0] * x + op[2][1] * y + op[2][2] * z + op[2][3]
    return n


def rot_angle(op: tuple) -> float:
    tr = op[0][0] + op[1][1] + op[2][2]
    return math.degrees(math.acos(max(-1.0, min(1.0, (tr - 1) / 2))))


def dist(a: Atom, b: Atom) -> float:
    return math.dist(a.xyz(), b.xyz())


def reps(atoms: list[Atom]) -> list[Atom]:
    """One or two representative atoms per residue, for fast contact counts."""
    return [a for a in atoms if a.name in ("CA", "P", "C4'")]


def contacts(xs: list[Atom], ys: list[Atom], cut: float = 8.0) -> tuple[int, float]:
    n, dmin = 0, 1e9
    for a in xs:
        for b in ys:
            d = dist(a, b)
            if d < cut:
                n += 1
            if d < dmin:
                dmin = d
    return n, dmin


IDENTITY = ((1.0, 0.0, 0.0, 0.0), (0.0, 1.0, 0.0, 0.0), (0.0, 0.0, 1.0, 0.0))


# --------------------------------------------------------------------------
# cleaning and checks
# --------------------------------------------------------------------------

def residues(atoms: list[Atom]) -> list[list[Atom]]:
    out: list[list[Atom]] = []
    for a in atoms:
        if out and out[-1][0].res == a.res:
            out[-1].append(a)
        else:
            out.append([a])
    return out


def clean(atoms: list[Atom], keep_mg: bool) -> tuple[list[Atom], list[str]]:
    kept, dropped, notes, modified = [], {}, [], {}
    for a in atoms:
        k = kind(a.resname)
        if a.element.upper() in ("H", "D"):
            continue
        if k == "protein":
            if a.resname == "MSE":
                a.resname = "MET"
                if a.name == "SE":
                    a.name, a.element = "SD", "S"
            kept.append(a)
        elif k == "rna":
            if a.resname in MODIFIED_NT:
                parent, extra = MODIFIED_NT[a.resname]
                modified[a.resname] = parent
                if a.name in extra:
                    continue
                a.resname = parent
            a.resname = RNA_NAMES[a.resname]
            if a.name == "OP3":       # 5'-terminal phosphate oxygen, not in CHARMM
                continue
            kept.append(a)
        elif k == "mg" and keep_mg:
            a.resname, a.name, a.element = "MG", "MG", "MG"
            kept.append(a)
        elif k != "water":
            dropped[a.resname] = dropped.get(a.resname, 0) + 1
    for m, p in modified.items():
        notes.append(f"modified nucleotide {m} reverted to {p}")
    if dropped:
        notes.append("dropped HETATM residues: " +
                     ", ".join(f"{k}({v} atoms)" for k, v in sorted(dropped.items())))
    return kept, notes


def chain_report(atoms: list[Atom]) -> list[str]:
    """Sequence, numbering range, breaks and missing heavy atoms, per chain."""
    lines = []
    chains: dict[str, list[list[Atom]]] = {}
    for r in residues(atoms):
        chains.setdefault(r[0].chain, []).append(r)
    for ch, rs in chains.items():
        k = kind(rs[0][0].resname)
        if k == "mg":
            lines.append(f"  chain {ch}: {len(rs)} Mg2+")
            continue
        seq = "".join(AA_ONE.get(r[0].resname, "X") if k == "protein"
                      else RNA_NAMES.get(r[0].resname, "X") for r in rs)
        lines.append(f"  chain {ch}: {k:<7} {len(rs):4d} res  "
                     f"{rs[0][0].resseq}..{rs[-1][0].resseq}  {seq}")
        occ = sum(a.occ for r in rs for a in r) / sum(len(r) for r in rs)
        if occ < 0.99:
            lines.append(f"    mean occupancy {occ:.2f} -- partially occupied (disorder "
                         f"or incomplete soaking); use a fully occupied copy if one exists, "
                         f"otherwise record it in docs/DECISIONS.md")
        link = ("C", "N", 2.0) if k == "protein" else ("O3'", "P", 2.1)
        for prev, cur in zip(rs, rs[1:]):
            a = next((x for x in prev if x.name == link[0]), None)
            b = next((x for x in cur if x.name == link[1]), None)
            if a is None or b is None or dist(a, b) > link[2]:
                gap = cur[0].resseq - prev[0].resseq - 1
                d = f"{dist(a, b):.1f} A" if a and b else "atom missing"
                lines.append(f"    BREAK after {prev[0].resname}{prev[0].resseq}"
                             f" ({gap} missing residue(s), {d})"
                             " -- model it before simulating")
        for r in rs:
            names = {a.name for a in r}
            if k == "protein":
                want = AA_HEAVY[r[0].resname] - {"OXT"}
            else:
                want = (SUGAR_PHOSPHATE - {"OP3", "O3P", "OP1", "OP2", "O1P", "O2P",
                                           "P"}) | BASE_HEAVY[r[0].resname]
            miss = sorted(want - names)
            if miss:
                lines.append(f"    missing heavy atoms in {r[0].resname}{r[0].resseq}:"
                             f" {' '.join(miss)}")
    return lines


# --------------------------------------------------------------------------
# RNA mutation
# --------------------------------------------------------------------------

def mutate_residue(res: list[Atom], new: str) -> list[Atom]:
    old = res[0].resname
    if old == new:
        return res
    out = []
    for a in res:
        if a.name in SUGAR_PHOSPHATE:
            pass
        elif (old in PURINES) == (new in PURINES):
            ring = PURINE_RING if new in PURINES else PYRIMIDINE_RING
            if a.name not in ring:
                continue
        elif old in PURINES:
            if a.name not in PUR_TO_PYR:
                continue
            a.name = PUR_TO_PYR[a.name]
        else:
            if a.name not in PYR_TO_PUR:
                continue
            a.name = PYR_TO_PUR[a.name]
        a.resname = new
        out.append(a)
    return out


def load_designs() -> dict[str, tuple[str, str]]:
    if not DESIGNS.exists():
        return {}
    out = {}
    for line in DESIGNS.read_text().splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("name\t"):
            continue
        f = [x.strip() for x in line.split("\t")]
        out[f[0]] = (f[1], f[2])
    return out


def mutate_chain(atoms: list[Atom], chain: str, target: str) -> tuple[list[Atom], list[str]]:
    rs = [r for r in residues(atoms) if r[0].chain == chain]
    if not rs:
        sys.exit(f"rna_designs.tsv names chain {chain}, which is not in the build")
    native = "".join(r[0].resname for r in rs)
    target = target.upper().replace("T", "U")
    if len(target) != len(native):
        sys.exit(f"RNA chain {chain}: target has {len(target)} nt but {len(native)} are "
                 f"resolved\n  native {native}\n  target {target}\n"
                 f"Write the target over the resolved residues only ('.' keeps one).")
    mut = {}
    for r, t in zip(rs, target):
        if t != "." and t != r[0].resname:
            if t not in BASE_HEAVY:
                sys.exit(f"invalid nucleotide '{t}' in target for chain {chain}")
            mut[r[0].res] = t
    out = []
    for r in residues(atoms):
        out.extend(mutate_residue(r, mut[r[0].res]) if r[0].res in mut else r)
    marks = "".join("|" if t in (".", n) else "*" for n, t in zip(native, target))
    shown = "".join(n if t == "." else t for n, t in zip(native, target))
    notes = [f"RNA chain {chain}: {len(mut)} mutation(s) "
             f"(residues {rs[0][0].resseq}..{rs[-1][0].resseq})",
             f"  native {native}", f"         {marks}", f"  built  {shown}"]
    return out, notes


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------

def cmd_show(src: str) -> None:
    path = resolve_source(src)
    atoms, ops, notes = read_pdb(path)
    atoms, n2 = clean(atoms, keep_mg=True)
    print(f"{path.relative_to(REPO) if path.is_relative_to(REPO) else path}")
    for n in notes + n2:
        print(f"  note: {n}")
    print("\n".join(chain_report(atoms)))
    print(f"  BIOMT operators (biomolecule 1): {len(ops)}")
    by_chain: dict[str, list[Atom]] = {}
    for a in reps(atoms):
        by_chain.setdefault(a.chain, []).append(a)
    if len(ops) > 1 and by_chain:
        print("  closest symmetry partners (contacts = representative atoms < 8 A):")
        for ch, xs in by_chain.items():
            best = []
            for ch2, ys in by_chain.items():
                for i, op in enumerate(ops):
                    moved = [apply(op, y) for y in ys]
                    n, dmin = contacts(xs, moved)
                    if dmin > 3.0 and n:
                        best.append((n, ch2, i + 1, rot_angle(op)))
            for n, ch2, i, ang in sorted(best, reverse=True)[:3]:
                print(f"    {ch} <- {ch2} op{i:<3} {n:4d} contacts  rotation {ang:5.1f} deg")


def cmd_build(system: str, keep_mg: bool | None) -> None:
    row = _registry.get(system, _registry.load())
    if row["source_pdb"].startswith("graft:"):
        return cmd_graft(system, row)
    spec = [c.strip() for c in row["chains"].split(",") if c.strip()]
    if not spec or not all(len(c) == 1 for c in spec):
        sys.exit(f"chains '{row['chains']}' for {system} cannot be built "
                 f"automatically; prepare input/{system}/{system}.pdb by hand")
    has_rna = row["has_rna"] == "yes"
    if keep_mg is None:
        keep_mg = has_rna
    path = resolve_source(row["source_pdb"], system)
    atoms, ops, notes = read_pdb(path)
    atoms, n2 = clean(atoms, keep_mg)
    notes += n2
    if not ops:
        ops = [IDENTITY]

    by_chain: dict[str, list[Atom]] = {}
    for a in atoms:
        by_chain.setdefault(a.chain, []).append(a)
    for c in spec:
        if c not in by_chain or all(kind(a.resname) == "mg" for a in by_chain[c]):
            sys.exit(f"chain {c} not found in {path.name} "
                     f"(has {', '.join(sorted(by_chain))}) -- run `show` first")

    # output chain IDs: protein A, B, C...; RNA R, S, T...; crystal Mg2+ M
    next_id = {"protein": iter("ABCDEFGHIJKL"), "rna": iter("RSTUVW")}
    assembled: list[Atom] = []
    used: set[tuple] = set()
    log = []
    for k, c in enumerate(spec):
        poly = [a for a in by_chain[c] if kind(a.resname) != "mg"]
        if k == 0:
            choice = (0, IDENTITY, 0, 0.0)
        else:
            cand = []
            for i, op in enumerate(ops):
                if (c, i) in used:
                    continue
                n, dmin = contacts(reps(assembled), [apply(op, a) for a in reps(poly)])
                if dmin > 3.0:              # < 3 A between CA/P means overlap, not contact
                    cand.append((n, -i, op))
            if not cand or max(cand)[0] == 0:
                sys.exit(f"no symmetry copy of chain {c} touches the assembly; check "
                         f"chains in config/systems.tsv against `show` output")
            n, mi, op = max(cand, key=lambda t: (t[0], t[1]))
            choice = (-mi, op, n, rot_angle(op))
        i, op, n, ang = choice
        used.add((c, i))
        new_id = next(next_id[kind(poly[0].resname)])
        moved = [apply(op, a) for a in by_chain[c]]
        for a in moved:
            a.chain = "M" if kind(a.resname) == "mg" else new_id
        polym = [a for a in moved if a.chain != "M"]
        mgs = [a for a in moved if a.chain == "M"]
        assembled.extend(polym)
        assembled.extend(mgs)
        log.append(f"chain {c} -> {new_id}: "
                   + ("as deposited" if k == 0 else
                      f"BIOMT op{i + 1}, {n} contacts, rotation {ang:.1f} deg"
                      + ("  (two-fold)" if abs(ang - 180) < 5 else "")))

    # crystallographic Mg2+ that ended up nowhere near the assembly are not ours
    polymer = [a for a in assembled if a.chain != "M"]
    mg = [a for a in assembled if a.chain == "M"
          and min(dist(a, b) for b in polymer) < 6.0]
    for j, a in enumerate(mg, 1):
        a.resseq, a.icode = j, " "
    if mg:
        log.append(f"{len(mg)} crystallographic Mg2+ kept within 6 A (chain M)")

    designs = load_designs()
    if system in designs:
        ch, target = designs[system]
        if target.upper().startswith("TODO"):
            sys.exit(f"config/rna_designs.tsv: sequence for {system} is still TODO")
        if target.lower() != "native":
            polymer, mnotes = mutate_chain(polymer, ch, target)
            log += mnotes
    elif has_rna and row["source_pdb"] not in ("model",) and \
            not row["source_pdb"].startswith("from:"):
        log.append("WARNING: RNA system with no row in config/rna_designs.tsv -- "
                   "native RNA kept")

    write_output(system, row, path, polymer, mg, notes, log)


def write_output(system: str, row: dict, path: pathlib.Path, polymer: list[Atom],
                 mg: list[Atom], notes: list[str], log: list[str]) -> None:
    out_dir = INPUT / system
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / f"{system}.pdb"
    lines = [f"REMARK   1 BruSafe input for {system}",
             f"REMARK   1 source {path.name}  chains {row['chains']}"]
    lines += [f"REMARK   1 {x}" for x in log]
    serial, last = 0, None
    for a in polymer + mg:
        if last is not None and a.chain != last:
            lines.append("TER")
        serial += 1
        lines.append(fmt_atom(serial, a, "HETATM" if a.chain == "M" else "ATOM  "))
        last = a.chain
    lines += ["TER", "END"]
    out.write_text("\n".join(lines) + "\n")

    print(f"{system}  <-  {path.relative_to(REPO) if path.is_relative_to(REPO) else path}")
    for x in notes + log:
        print(f"  {x}")
    print("\n".join(chain_report(polymer + mg)))
    print(f"wrote {out.relative_to(REPO)}  ({serial} heavy atoms)")
    if row["has_rna"] == "yes":
        print("next: CHARMM-GUI Solution Builder with this PDB (docs/TOPOLOGIES.md),\n"
              f"      unzip its gromacs/ folder into input/{system}/gromacs/, then\n"
              f"      ./scripts/01_prepare.sh {system} --source charmm-gui")
    else:
        print(f"next: ./scripts/01_prepare.sh {system}")


# --------------------------------------------------------------------------
# graft: single-chain model + RNA from a built reference complex (S9)
# --------------------------------------------------------------------------

def align(a: str, b: str) -> list[tuple[int, int]]:
    """Global alignment (Needleman-Wunsch, linear gaps), returns aligned index
    pairs. Enough to place the two CP copies of a single-chain construct on
    the two chains of a dimer; the linker falls out as an insertion."""
    n, m, gap = len(a), len(b), -2
    sc = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        sc[i][0] = i * gap
    for j in range(1, m + 1):
        sc[0][j] = j * gap
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sc[i][j] = max(sc[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1),
                           sc[i - 1][j] + gap, sc[i][j - 1] + gap)
    pairs, i, j = [], n, m
    while i > 0 and j > 0:
        if sc[i][j] == sc[i - 1][j - 1] + (2 if a[i - 1] == b[j - 1] else -1):
            pairs.append((i - 1, j - 1)); i -= 1; j -= 1
        elif sc[i][j] == sc[i - 1][j] + gap:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def jacobi4(a: list[list[float]]) -> tuple[list[float], list[list[float]]]:
    """Eigen-decomposition of a symmetric 4x4 matrix (cyclic Jacobi)."""
    a = [r[:] for r in a]
    v = [[1.0 if i == j else 0.0 for j in range(4)] for i in range(4)]
    for _ in range(100):
        off = sum(a[i][j] ** 2 for i in range(4) for j in range(4) if i != j)
        if off < 1e-18:
            break
        for p in range(3):
            for q in range(p + 1, 4):
                if abs(a[p][q]) < 1e-30:
                    continue
                th = (a[q][q] - a[p][p]) / (2 * a[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1); s_ = t * c
                for k in range(4):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p], a[k][q] = c * akp - s_ * akq, s_ * akp + c * akq
                for k in range(4):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k], a[q][k] = c * apk - s_ * aqk, s_ * apk + c * aqk
                for k in range(4):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p], v[k][q] = c * vkp - s_ * vkq, s_ * vkp + c * vkq
    return [a[i][i] for i in range(4)], v


def superpose(mobile: list[tuple], target: list[tuple]) -> tuple:
    """Least-squares rigid fit of mobile onto target (Horn quaternions).
    Returns a 3x4 operator usable with apply()."""
    n = len(mobile)
    cm = [sum(p[k] for p in mobile) / n for k in range(3)]
    ct = [sum(p[k] for p in target) / n for k in range(3)]
    S = [[0.0] * 3 for _ in range(3)]
    for p, q in zip(mobile, target):
        for i in range(3):
            for j in range(3):
                S[i][j] += (p[i] - cm[i]) * (q[j] - ct[j])
    (xx, xy, xz), (yx, yy, yz), (zx, zy, zz) = S
    N = [[xx + yy + zz, yz - zy, zx - xz, xy - yx],
         [yz - zy, xx - yy - zz, xy + yx, zx + xz],
         [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
         [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]
    w, v = jacobi4(N)
    k = max(range(4), key=lambda i: w[i])
    q0, q1, q2, q3 = (v[i][k] for i in range(4))
    R = [[q0*q0 + q1*q1 - q2*q2 - q3*q3, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
         [2*(q1*q2 + q0*q3), q0*q0 - q1*q1 + q2*q2 - q3*q3, 2*(q2*q3 - q0*q1)],
         [2*(q1*q3 - q0*q2), 2*(q2*q3 + q0*q1), q0*q0 - q1*q1 - q2*q2 + q3*q3]]
    t = [ct[i] - sum(R[i][j] * cm[j] for j in range(3)) for i in range(3)]
    return tuple(tuple(R[i]) + (t[i],) for i in range(3))


def rmsd_after(op: tuple, mobile: list[tuple], target: list[tuple]) -> list[float]:
    out = []
    for p, q in zip(mobile, target):
        x = tuple(sum(op[i][j] * p[j] for j in range(3)) + op[i][3] for i in range(3))
        out.append(math.dist(x, q))
    return out


def cmd_graft(system: str, row: dict) -> None:
    """Place a single-chain protein model on the protein chains of an already
    built reference complex and take that complex's RNA (and Mg2+)."""
    ref_sys = row["source_pdb"].split(":", 1)[1]
    ref_path = INPUT / ref_sys / f"{ref_sys}.pdb"
    model_path = INPUT / system / "model.pdb"
    if not ref_path.exists():
        sys.exit(f"build {ref_sys} first ({ref_path.relative_to(REPO)} missing)")
    if not model_path.exists():
        sys.exit(f"put the single-chain model at {model_path.relative_to(REPO)}")

    ref, _, _ = read_pdb(ref_path)
    model, _, notes = read_pdb(model_path)
    model, n2 = clean(model, keep_mg=False)
    notes += n2
    prot_model = [a for a in model if kind(a.resname) == "protein"]
    if len({a.chain for a in prot_model}) != 1:
        sys.exit(f"{model_path.name}: expected ONE protein chain (the single-chain "
                 f"construct), found {sorted({a.chain for a in prot_model})}")
    ref_prot = [a for a in ref if kind(a.resname) == "protein"]
    ref_rest = [a for a in ref if kind(a.resname) in ("rna", "mg")]
    if not ref_rest:
        sys.exit(f"{ref_sys} has no RNA to graft")

    def ca_list(atoms):
        rs = [r for r in residues(atoms) if any(a.name == "CA" for a in r)]
        return ("".join(AA_ONE.get(r[0].resname, "X") for r in rs),
                [next(a for a in r if a.name == "CA") for r in rs])

    mseq, mca = ca_list(prot_model)
    rseq, rca = ca_list(ref_prot)        # chain A then chain B, in file order
    pairs = align(mseq, rseq)
    if len(pairs) < 0.8 * len(rseq):
        sys.exit(f"only {len(pairs)} of {len(rseq)} reference residues align to the "
                 f"model -- is this the right construct?")
    mob = [mca[i].xyz() for i, _ in pairs]
    tgt = [rca[j].xyz() for _, j in pairs]
    op = superpose(mob, tgt)
    d = rmsd_after(op, mob, tgt)
    rms_all = math.sqrt(sum(x * x for x in d) / len(d))
    core = [k for k, x in enumerate(d) if x < 2.0]     # refit on the rigid core
    if len(core) >= 0.5 * len(d):
        op = superpose([mob[k] for k in core], [tgt[k] for k in core])
        d = rmsd_after(op, mob, tgt)
    rms_core = math.sqrt(sum(d[k] ** 2 for k in core) / max(1, len(core)))

    placed = [apply(op, a) for a in prot_model]
    for a in placed:
        a.chain = "A"
    rna = [a.copy() for a in ref_rest if kind(a.resname) == "rna"]
    mg = [a.copy() for a in ref_rest if kind(a.resname) == "mg"]
    aligned = {i for i, _ in pairs}
    insert = [k for k in range(len(mseq)) if k not in aligned]
    clashes = sum(1 for a in placed for b in rna if dist(a, b) < 2.0)

    # per-residue deviation, worst first, so a mis-folded copy is obvious
    worst = sorted(((x, mca[pairs[k][0]]) for k, x in enumerate(d)),
                   key=lambda t: -t[0])[:5]
    log = [f"single-chain model {model_path.name} superposed on {ref_sys} protein "
           f"({len(pairs)} CA pairs)",
           f"CA RMSD all aligned {rms_all:.2f} A, rigid core ({len(core)} CA) "
           f"{rms_core:.2f} A",
           "largest deviations: " + ", ".join(f"{a.resname}{a.resseq} {x:.1f} A"
                                              for x, a in worst),
           f"model residues with no counterpart in the dimer (linker etc.): "
           f"{len(insert)}" + (f", {mseq[insert[0]:insert[-1] + 1]!s:.40}"
                               if insert else ""),
           f"RNA and Mg2+ taken from {ref_sys} unchanged",
           f"protein-RNA heavy-atom pairs closer than 2.0 A: {clashes}"]
    if rms_core > 2.0:
        log.append("WARNING: core RMSD > 2 A -- the two CP copies in the model do "
                   "not sit like the native dimer; fix the model before simulating")
    if clashes:
        log.append("WARNING: the linker or a loop overlaps the RNA -- inspect and "
                   "rebuild that region before CHARMM-GUI")
    write_output(system, row, model_path, placed + rna, mg, notes, log)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fetch"); p.add_argument("pdbid")
    p = sub.add_parser("show");  p.add_argument("source")
    p = sub.add_parser("build"); p.add_argument("system")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--keep-mg", dest="keep_mg", action="store_true", default=None,
                   help="keep crystallographic Mg2+ (default for RNA systems)")
    g.add_argument("--no-mg", dest="keep_mg", action="store_false")
    a = ap.parse_args()
    if a.cmd == "fetch":
        print(fetch(a.pdbid))
    elif a.cmd == "show":
        cmd_show(a.source)
    else:
        cmd_build(a.system, a.keep_mg)


if __name__ == "__main__":
    main()
