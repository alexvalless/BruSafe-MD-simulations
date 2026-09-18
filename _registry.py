#!/usr/bin/env python3
"""
Registry reader and campaign planner for the BruSafe MD campaign.

Subcommands
-----------
  field <system> <column>     print one registry field (used by lib.sh)
  list [--tier N]             print system names
  plan  --machines M [--tier N]
                              expand systems to replicas, balance across
                              machines by estimated GPU cost, emit job lists
  status --runs DIR           scan run directories and report progress

The planner exists because the realistic failure mode with ten unmanaged
desktops is not "we ran out of compute" -- it is "two boxes ran the same
replica and one system never launched at all".

No third-party dependencies.
"""

from __future__ import annotations

import argparse
import pathlib
import re
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
REGISTRY = REPO / "config" / "systems.tsv"


# --------------------------------------------------------------------------
# registry
# --------------------------------------------------------------------------

def load(path: pathlib.Path = REGISTRY) -> list[dict]:
    if not path.exists():
        sys.exit(f"registry not found: {path}")
    header: list[str] | None = None
    rows: list[dict] = []
    for raw in path.read_text().splitlines():
        line = raw.rstrip("\n")
        if not line.strip():
            continue
        if line.startswith("#"):
            continue
        if header is None:
            # first non-comment line is the header
            header = [h.strip() for h in line.split("\t")]
            if header[0] != "name":
                sys.exit("registry header must start with a 'name' column")
            continue
        fields = line.split("\t")
        if len(fields) != len(header):
            sys.exit(f"malformed row (expected {len(header)} columns): {line[:60]}")
        rows.append(dict(zip(header, (f.strip() for f in fields))))
    if not rows:
        sys.exit("registry is empty")
    return rows


def get(name: str, rows: list[dict]) -> dict:
    for r in rows:
        if r["name"] == name:
            return r
    sys.exit(f"system not in registry: {name}")


# --------------------------------------------------------------------------
# planning
# --------------------------------------------------------------------------

def expand(rows: list[dict], max_tier: int) -> list[dict]:
    """One job per replica, with an estimated cost in ns."""
    jobs = []
    for r in rows:
        if int(r["tier"]) > max_tier:
            continue
        ns = int(r["ns"])
        for rep in range(1, int(r["replicas"]) + 1):
            jobs.append({
                "system": r["name"],
                "tier": int(r["tier"]),
                "replica": rep,
                "ns": ns,
                # RNA systems carry more atoms, so cost more per ns.
                "cost": ns * (1.6 if r["has_rna"] == "yes" else 1.0),
                "consumer": r["consumer"],
            })
    return jobs


def assign(jobs: list[dict], machines: int) -> dict[int, list[dict]]:
    """Longest-processing-time-first: sort by cost desc, always feed the
    least-loaded machine. Simple, and within 4/3 of optimal makespan."""
    loads = [0.0] * machines
    buckets: dict[int, list[dict]] = {i: [] for i in range(machines)}
    for job in sorted(jobs, key=lambda j: (-j["cost"], j["system"], j["replica"])):
        m = min(range(machines), key=lambda i: (loads[i], i))
        buckets[m].append(job)
        loads[m] += job["cost"]
    return buckets


def emit_plan(buckets: dict[int, list[dict]], nsday: float) -> None:
    total = 0.0
    print(f"{'machine':<9}{'jobs':<6}{'cost (ns-eq)':<14}{'est. days':<11}")
    print("-" * 72)
    for m, jobs in buckets.items():
        cost = sum(j["cost"] for j in jobs)
        total += cost
        print(f"gpu{m:<6}{len(jobs):<6}{cost:<14.0f}{cost / nsday:<11.1f}")
    print("-" * 72)
    makespan = max(sum(j["cost"] for j in jobs) for jobs in buckets.values()) / nsday
    print(f"aggregate {total:.0f} ns-eq | makespan {makespan:.1f} days "
          f"at {nsday:.0f} ns/day/GPU")
    print()
    for m, jobs in buckets.items():
        print(f"# ---- gpu{m} ----")
        for j in jobs:
            print(f"gpu{m}\t{j['system']}\t{j['replica']}\t{j['ns']}\t{j['consumer']}")


# --------------------------------------------------------------------------
# status
# --------------------------------------------------------------------------

PERF = re.compile(r"^\s*Performance:\s+([\d.]+)")
STEP = re.compile(r"^\s+Step\s+Time\s*$")


def scan(runs: pathlib.Path, rows: list[dict]) -> None:
    print(f"{'system':<20}{'rep':<5}{'stage':<14}{'progress':<12}{'ns/day':<9}")
    print("-" * 72)
    missing = []
    for r in rows:
        for rep in range(1, int(r["replicas"]) + 1):
            d = runs / r["name"] / f"rep{rep}"
            if not d.exists():
                missing.append(f"{r['name']}/rep{rep}")
                print(f"{r['name']:<20}{rep:<5}{'NOT STARTED':<14}{'-':<12}{'-':<9}")
                continue
            stage, prog, perf = inspect(d, int(r["ns"]))
            print(f"{r['name']:<20}{rep:<5}{stage:<14}{prog:<12}{perf:<9}")
    print("-" * 72)
    if missing:
        print(f"{len(missing)} replica(s) never launched:")
        for m in missing:
            print(f"  {m}")
    else:
        print("every registered replica has a run directory")


def inspect(d: pathlib.Path, target_ns: int) -> tuple[str, str, str]:
    if (d / "prod.gro").exists():
        log = d / "prod.log"
        perf = "-"
        if log.exists():
            for line in log.read_text(errors="ignore").splitlines():
                m = PERF.match(line)
                if m:
                    perf = m.group(1)
        return "DONE", f"{target_ns} ns", perf
    for stage, gro in (("prod", "prod.log"), ("npt_free", "npt_free.gro"),
                       ("npt", "npt.gro"), ("nvt", "nvt.gro"), ("em", "em.gro")):
        if (d / gro).exists():
            if stage == "prod":
                ns = last_ns(d / "prod.log")
                pct = f"{ns:.0f}/{target_ns} ns" if ns is not None else "running"
                return "prod", pct, "-"
            return stage, "equil", "-"
    return "empty", "-", "-"


def last_ns(log: pathlib.Path) -> float | None:
    """Last simulation time reported in an mdrun log, in ns."""
    if not log.exists():
        return None
    lines = log.read_text(errors="ignore").splitlines()
    for i in range(len(lines) - 1, 0, -1):
        if STEP.match(lines[i]):
            parts = lines[i + 1].split() if i + 1 < len(lines) else []
            if len(parts) >= 2:
                try:
                    return float(parts[1]) / 1000.0
                except ValueError:
                    return None
    return None


# --------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("field"); p.add_argument("system"); p.add_argument("column")
    p = sub.add_parser("list");  p.add_argument("--tier", type=int, default=3)
    p = sub.add_parser("plan")
    p.add_argument("--machines", type=int, required=True)
    p.add_argument("--tier", type=int, default=1)
    p.add_argument("--nsday", type=float, default=250.0,
                   help="measured ns/day per GPU from bench_gpu.sh")
    p = sub.add_parser("status"); p.add_argument("--runs", required=True)

    a = ap.parse_args()
    rows = load()

    if a.cmd == "field":
        row = get(a.system, rows)
        if a.column not in row:
            sys.exit(f"no such column: {a.column}")
        print(row[a.column])
    elif a.cmd == "list":
        for r in rows:
            if int(r["tier"]) <= a.tier:
                print(r["name"])
    elif a.cmd == "plan":
        jobs = expand(rows, a.tier)
        if not jobs:
            sys.exit("no jobs at that tier")
        emit_plan(assign(jobs, a.machines), a.nsday)
    elif a.cmd == "status":
        scan(pathlib.Path(a.runs), rows)


if __name__ == "__main__":
    try:
        main()
    except BrokenPipeError:
        # piping into head/less is normal usage
        try:
            sys.stdout.close()
        finally:
            sys.exit(0)
