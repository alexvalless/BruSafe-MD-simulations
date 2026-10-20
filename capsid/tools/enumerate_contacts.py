#!/usr/bin/env python3
"""
Contact-schedule enumeration for the MS2 T=3 capsid (PDB 2MS2).

Pipeline
  1. Parse chains A, B, C of the asymmetric unit (heavy atoms only).
  2. Apply the 60 BIOMT operators (REMARK 350) -> 180 chains.
  3. Dimers: A/B within each ASU (60) and C/C across ASUs (30)  -> 90 assembly units.
  4. Dimer-dimer interface weight = number of distinct residue pairs with a
     heavy-atom distance < CUTOFF.  An interface "edge" exists if weight >= EDGE_MIN.
  5. Growth paths on the 90-node dimer graph; contact schedule c_k = number of NEW
     interfaces formed when the k-th dimer is added (k = 2..90, 89 events);
     degeneracy proxy g_k = number of boundary dimers that could be added at step k.

Usage: python3 -I enumerate_contacts.py 2MS2.pdb outdir
"""
import sys, json, collections
import numpy as np
from scipy.spatial import cKDTree

CUTOFF = 4.5
EDGE_MIN_DEFAULT = 5
N_PATHS = 4000
rng = np.random.default_rng(20261019)


def parse(path):
    atoms, biomt = [], collections.defaultdict(dict)
    for line in open(path):
        if line.startswith("ATOM") and line[76:78].strip() != "H":
            atoms.append((line[21], int(line[22:26]), float(line[30:38]),
                          float(line[38:46]), float(line[46:54])))
        elif line.startswith("REMARK 350   BIOMT"):
            row = int(line[18]) - 1
            idx = int(line[19:23])
            vals = [float(x) for x in line[23:].split()]
            biomt[idx][row] = vals
    chains = collections.defaultdict(list)
    for c, r, x, y, z in atoms:
        chains[c].append((r, x, y, z))
    ops = []
    for idx in sorted(biomt):
        M = np.array([biomt[idx][i][:3] for i in range(3)])
        t = np.array([biomt[idx][i][3] for i in range(3)])
        ops.append((M, t))
    return chains, ops


def build_capsid(chains, ops):
    res, xyz, cid = [], [], []
    labels = []
    for oi, (M, t) in enumerate(ops):
        for c in "ABC":
            a = np.array(chains[c])
            p = a[:, 1:4] @ M.T + t
            res.append(a[:, 0].astype(int))
            xyz.append(p)
            lab = len(labels)
            labels.append((oi, c))
            cid.append(np.full(len(a), lab))
    return np.concatenate(res), np.vstack(xyz), np.concatenate(cid), labels


def chain_contacts(res, xyz, cid):
    tree = cKDTree(xyz)
    pairs = tree.query_pairs(CUTOFF, output_type="ndarray")
    a, b = cid[pairs[:, 0]], cid[pairs[:, 1]]
    keep = a != b
    pairs, a, b = pairs[keep], a[keep], b[keep]
    lo, hi = np.minimum(a, b), np.maximum(a, b)
    ra = np.where(a < b, res[pairs[:, 0]], res[pairs[:, 1]])
    rb = np.where(a < b, res[pairs[:, 1]], res[pairs[:, 0]])
    key = np.stack([lo, hi, ra, rb], axis=1)
    key = np.unique(key, axis=0)                      # unique residue pairs
    out = collections.Counter()
    for l, h in zip(key[:, 0], key[:, 1]):
        out[(int(l), int(h))] += 1
    atoms = collections.Counter()
    for l, h in zip(lo, hi):
        atoms[(int(l), int(h))] += 1
    return out, atoms


def main(pdb, outdir):
    chains, ops = parse(pdb)
    assert len(ops) == 60, len(ops)
    res, xyz, cid, labels = build_capsid(chains, ops)
    nchain = len(labels)
    print(f"chains in capsid: {nchain}; heavy atoms: {len(xyz)}")
    rescont, atomcont = chain_contacts(res, xyz, cid)

    idx = {lab: i for i, lab in enumerate(labels)}
    # --- dimers
    dimers = []                                        # (type, chain_i, chain_j)
    # A/B partners: for each A chain the B chain with the largest interface
    # (in 2MS2 the quasi-2-fold partner lies in a different ASU than the A chain)
    achains = [idx[(oi, "A")] for oi in range(60)]
    bchains = [idx[(oi, "B")] for oi in range(60)]
    usedB = set()
    for a in achains:
        cands = [(rescont.get((min(a, b), max(a, b)), 0), b) for b in bchains if b not in usedB]
        w, b = max(cands)
        assert w > 100, w
        dimers.append(("AB", a, b))
        usedB.add(b)
    cchains = [idx[(oi, "C")] for oi in range(60)]
    used = set()
    for c in cchains:
        if c in used:
            continue
        best, bw = None, 0
        for d in cchains:
            if d == c or d in used:
                continue
            w = rescont.get((min(c, d), max(c, d)), 0)
            if w > bw:
                best, bw = d, w
        assert best is not None
        dimers.append(("CC", c, best))
        used.update([c, best])
    assert len(dimers) == 90, len(dimers)
    chain2dimer = {}
    for di, (_, i, j) in enumerate(dimers):
        chain2dimer[i] = di
        chain2dimer[j] = di
    assert len(chain2dimer) == 180
    cc_w = [rescont[(min(i, j), max(i, j))] for t, i, j in dimers if t == "CC"]
    ab_w = [rescont[(min(i, j), max(i, j))] for t, i, j in dimers if t == "AB"]
    print(f"intra-dimer residue contacts: AB {np.mean(ab_w):.1f} (n=60), CC {np.mean(cc_w):.1f} (n=30)")

    # --- dimer-dimer weights
    W = collections.Counter()
    WA = collections.Counter()
    for (i, j), w in rescont.items():
        di, dj = chain2dimer[i], chain2dimer[j]
        if di != dj:
            W[(min(di, dj), max(di, dj))] += w
            WA[(min(di, dj), max(di, dj))] += atomcont[(i, j)]
    ws = sorted(W.values())
    print("number of dimer pairs with any contact:", len(ws))
    print("weight histogram (residue pairs):", np.histogram(ws, bins=[0, 3, 5, 10, 20, 40, 80, 1000])[0].tolist())

    results = {}
    for emin in (3, 5, 10, 20):
        edges = [e for e, w in W.items() if w >= emin]
        deg = collections.Counter()
        for a, b in edges:
            deg[a] += 1
            deg[b] += 1
        results[emin] = (len(edges), collections.Counter(deg.values()))
        print(f"EDGE_MIN={emin}: edges={len(edges)} degree distribution={dict(sorted(results[emin][1].items()))}")

    emin = EDGE_MIN_DEFAULT
    edges = {e: w for e, w in W.items() if w >= emin}
    nbr = collections.defaultdict(dict)
    for (a, b), w in edges.items():
        nbr[a][b] = w
        nbr[b][a] = w
    types = [d[0] for d in dimers]
    etype = collections.Counter()
    eweights = collections.defaultdict(list)
    for (a, b), w in edges.items():
        k = "-".join(sorted([types[a], types[b]]))
        etype[k] += 1
        eweights[k].append(w)
    print("edge classes:", dict(etype))
    for k, v in eweights.items():
        print(f"  {k}: mean residue-pair weight {np.mean(v):.1f}  range {min(v)}-{max(v)}")
    total_edges = len(edges)
    print("total interfaces (edges):", total_edges)

    # --- growth paths
    def grow(strategy):
        start = int(rng.integers(90))
        cluster = {start}
        c_sched, g_sched = [], []
        while len(cluster) < 90:
            cand = collections.Counter()
            for n in cluster:
                for m in nbr[n]:
                    if m not in cluster:
                        cand[m] += 1               # contacts it would form
            g_sched.append(len(cand))
            items = list(cand.items())
            if strategy == "random":
                m, c = items[int(rng.integers(len(items)))]
            elif strategy == "max":
                mx = max(c for _, c in items)
                pool = [x for x in items if x[1] == mx]
                m, c = pool[int(rng.integers(len(pool)))]
            else:
                mn = min(c for _, c in items)
                pool = [x for x in items if x[1] == mn]
                m, c = pool[int(rng.integers(len(pool)))]
            c_sched.append(c)
            cluster.add(m)
        return np.array(c_sched), np.array(g_sched)

    sched = {}
    for strat in ("random", "max", "min"):
        C, G = zip(*[grow(strat) for _ in range(N_PATHS)])
        C, G = np.array(C), np.array(G)
        sched[strat] = dict(c_mean=C.mean(0), c_sd=C.std(0), g_mean=G.mean(0))
        assert np.all(C.sum(1) == total_edges), (strat, C.sum(1)[:5])
        print(f"{strat:6s}: first={C.mean(0)[0]:.2f} last={C.mean(0)[-1]:.2f} sum={C.mean(0).sum():.1f}")

    import os
    os.makedirs(outdir, exist_ok=True)
    k = np.arange(2, 91)
    with open(f"{outdir}/contact_schedule.csv", "w") as f:
        f.write("k_dimers_after_addition,c_random_mean,c_random_sd,c_max,c_min,g_boundary_random_mean,g_boundary_max,g_boundary_min\n")
        for i, kk in enumerate(k):
            f.write(f"{kk},{sched['random']['c_mean'][i]:.3f},{sched['random']['c_sd'][i]:.3f},"
                    f"{sched['max']['c_mean'][i]:.3f},{sched['min']['c_mean'][i]:.3f},"
                    f"{sched['random']['g_mean'][i]:.2f},{sched['max']['g_mean'][i]:.2f},{sched['min']['g_mean'][i]:.2f}\n")
    json.dump(dict(structure="2MS2", cutoff_A=CUTOFF, edge_min_residue_pairs=emin,
                   total_interfaces=total_edges, edge_classes=dict(etype),
                   mean_weight_by_class={k: float(np.mean(v)) for k, v in eweights.items()},
                   contact_schedule_random=sched["random"]["c_mean"].tolist(),
                   contact_schedule_max=sched["max"]["c_mean"].tolist(),
                   contact_schedule_min=sched["min"]["c_mean"].tolist(),
                   boundary_sites_random=sched["random"]["g_mean"].tolist()),
              open(f"{outdir}/contact_schedule.json", "w"), indent=1)
    print("wrote", outdir)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
