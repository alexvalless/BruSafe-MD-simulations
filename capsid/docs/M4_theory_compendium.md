
# Module 4 — Assembly: Theoretical Compendium

**BruSafe / CAPSID · iGEM 2026 · the letter A in CAPSID**

*A single coherent derivation of every equation Module 4 rests on, in the order the physics
requires, with each parameter traced to a Crossref-verified source and each numerical claim
reproduced from running code.*

---

## 0. What this document is, and what it is not

Module 4 answers one question: **given a pool of MS2 coat-protein dimers, what fraction
closes into well-formed T=3 capsids, and at what dimer concentration is that fraction
maximal?**

That question has been attacked at four different levels of description, each of which
exists in this project as running code rather than as prose. This compendium is the
theory behind all four, assembled as one argument:

| Layer | Level of description | State variable | Answers |
|---|---|---|---|
| **L1** | Deterministic rate equations | Concentration of each partial shell size | Yield vs concentration and time |
| **L2** | Classical nucleation theory | Cluster size as a continuous coordinate | *Where* in parameter space assembly is possible at all |
| **L3** | Coarse-grained Monte Carlo | Positions and orientations of 90 subunits | The geometry L1 must assume |
| **L4** | Stochastic co-assembly with cargo | Occupancy of each face plus a threading path | How packaging signals change the answer |

This document is not a summary of the four existing Module 4 documents. It is the
derivation chain, written once, with the layers' points of contact and their one published
disagreement made explicit. **§9 is the wet-lab contract**: which measurement enters,
which parameter it moves, and what the model is obliged to do when the number disagrees.

Each layer, the bridge between them and the fitting routine ship as **runnable notebooks**
(`M4_L1_rate_equations.ipynb`, `M4_L2_nucleation_phase_diagram.ipynb`,
`M4_L3_coarsegrained_MC.ipynb`, `M4_L4_packaging_signal_coassembly.ipynb`,
`M4_bridge_CG_to_ODE_to_bench.ipynb`, `M4_FIT_to_bench_data.ipynb`). Each writes out the
module it needs and runs top to bottom on a fresh machine; see their `README.md`. The experimental designs, literature landscape and Spanish-language
team material live in their own artifacts and are referenced, not repeated.

Every equation below is either derived here from a stated physical assumption or lifted
verbatim from a cited paper. Where a number appears, it was read from the code that
produces it, not from memory — the project rule is that no value reaches a deliverable
without resolving against its source.

---

## 1. Notation and the constants shared by all four layers

Throughout, a **subunit** means one coat-protein dimer, the obligate structural unit of
MS2 established crystallographically by Golmohammadi and colleagues and refined for
recombinant assembly by Ni and colleagues. A T=3 shell contains

$$N = 90 \text{ dimers} = 180 \text{ coat-protein monomers}$$

Free energies are in kcal/mol and the thermal scale is

$$RT = 0.593\ \mathrm{kcal/mol} \qquad (T = 310\ \mathrm{K})$$

so 1 kcal/mol is 1.69 $k_BT$. Concentrations are in µM unless stated. The sign convention
is thermodynamic: $\Delta G < 0$ favours the bound state, and a dissociation constant maps
to a binding free energy by

$$\Delta G_{bind} = RT \ln K_d$$

with $K_d$ expressed relative to the standard state. This single relation is what lets a
measured affinity enter a kinetic model, and it is used in L1, L4 and Module 5 alike.

Two energies recur and must be kept distinct:

- $g_{assoc} > 0$ — the **association penalty** per subunit added: the translational and
  rotational entropy a free dimer surrenders on joining a shell. Zlotnick's foundational
  analysis of icosahedral assembly identifies this term as the origin of the nucleation
  barrier.
- $g_{contact} < 0$ — the free energy of **one subunit-subunit contact**. Ceres and Zlotnick
  measured it directly for hepatitis B capsids and found it remarkably weak, a few kcal/mol,
  which is precisely why assembly is reversible enough to be productive.

Assembly is possible because a subunit joining an existing shell forms *more than one*
contact, so a single penalty is paid against several favourable contacts. Everything in
Module 4 is an elaboration of that one sentence.

---

## 2. Layer L1 — the nucleation-elongation cascade

### 2.1 Upstream: the monomer-dimer equilibrium

Module 4 takes a total coat-protein concentration from Module 3 and must convert it to
*free dimer*, since the dimer is the assembling species. With dimerization

$$2\,M \rightleftharpoons D, \qquad K_d^{dim} = \frac{[M]^2}{[D]}$$

and total protein $P_{tot} = [M] + 2[D]$, eliminating $[D]$ gives a quadratic in $[M]$
whose physical root is

$$[M] = \frac{-K_d^{dim} + \sqrt{(K_d^{dim})^2 + 8 K_d^{dim} P_{tot}}}{4}$$

MS2 dimerization is tight — the folded dimer is the thermodynamic ground state of the
isolated chain — so with $K_d^{dim} = 0.01$ µM and $P_{tot} = 8$ µM the partition is
$[D] = 3.90$ µM, **97.5% of protein mass already dimeric**. Loosening $K_d^{dim}$ by two
orders of magnitude to 1.0 µM only drops that to 77.9%. The conclusion transfers to the
model as a licence: within the plausible range for MS2, Module 4 may take
$[D]_0 \approx P_{tot}/2$ and ignore the monomer pool. Module 3's steady state
($m_{ss} = 38.5$ nM mRNA, 7.96 µM monomer-equivalent) therefore enters Module 4 as
$[D]_0 = 3.98$ µM.

### 2.2 The state space and the master equation

Let $x_j(t)$ be the molar concentration of shells containing exactly $j$ subunits, with
$x_1$ the free dimer. Growth proceeds one subunit at a time — the assumption that no two
partial shells fuse, justified for capsid assembly by the strong asymmetry between the
abundant free subunit and the rare intermediate. Then the flux through step $j \to j{+}1$ is

$$J_j = f_j\, x_j x_1 - b_{j+1}\, x_{j+1}$$

and the cascade is the $N$-dimensional system

$$\frac{dx_j}{dt} = J_{j-1} - J_j$$

closed by $J_0 \equiv 0$, $J_N \equiv 0$ (a complete shell is inert), and by the free-subunit
balance which must count every subunit consumed anywhere in the ladder:

$$\frac{dx_1}{dt} = -\sum_{j=1}^{N-1} J_j - J_1$$

Total mass $\sum_j j\,x_j$ is conserved by construction; in the production runs it holds to
$1.5 \times 10^{-15}$ relative error, which is the numerical check that the flux bookkeeping
is right.

### 2.3 The contact ladder — where the geometry enters

The forward rate is a diffusion-limited on-rate, taken size-independent:

$$f_j = k_a \quad \text{for all } j$$

All the structure lives in the reverse rates, and they are **not free parameters**. A
subunit added to make the $j$-mer forms $c(j)$ new contacts, and the free-energy change of
that addition is

$$\Delta G_{add}(j) = g_{assoc} + c(j)\, g_{contact}$$

The production model encodes $c(j)$ as a three-regime step function reflecting the
lattice: an isolated subunit landing on a small cluster has one neighbour; a subunit
filling the edge of a growing cap has two; a subunit closing one of the final vertices is
surrounded on three sides.

$$c(j) = 1 \ \ (j \leq n_c), \qquad c(j) = 2 \ \ (n_c < j < N-8), \qquad c(j) = 3 \ \ (j \geq N-8)$$

With the literature energies $g_{assoc} = +6.0$ and $g_{contact} = -4.0$ kcal/mol this gives
a ladder that is **uphill before the nucleus and downhill after it**:

| Regime | Contacts c(j) | ΔG(add) | Meaning |
|---|---|---|---|
| $j \leq 4$ | 1 | **+2.0** kcal/mol | Nucleation is rare — this is the barrier |
| $4 < j < 82$ | 2 | **−2.0** kcal/mol | Elongation is favourable and fast |
| $j \geq 82$ | 3 | **−6.0** kcal/mol | Closure is nearly irreversible |

### 2.4 Detailed balance fixes the reverse rates

At equilibrium each step must satisfy $f_j x_j x_1 = b_{j+1} x_{j+1}$, so the ratio of
rates is the equilibrium constant of that step:

$$\frac{b_{j+1}}{f_j} = K_d(j) = \exp\!\left(\frac{\Delta G_{add}(j)}{RT}\right)$$

hence

$$b_{j+1} = k_a \exp\!\left(\frac{\Delta G_{add}(j)}{RT}\right)$$

**This is the defensible core of Module 4.** Ninety reverse rates are determined by two
energies and a contact count. Nothing in the ladder can be tuned to produce a desired
yield curve without moving $g_{assoc}$ or $g_{contact}$ out of their independently measured
ranges. A model that fitted 90 off-rates could reproduce any observation and would predict
nothing.

The cumulative profile follows by summation. With $B(j) = \sum_{i<j} c(i)$ the total bonds
in a $j$-mer,

$$\Delta G(j) = j\, g_{assoc} + B(j)\, g_{contact}$$

whose maximum defines the critical nucleus. For the step ladder that maximum sits at
$j = n_c = 4$ by construction — an *input*. Section 5 replaces it with an emergent value.

### 2.5 What L1 predicts: the concentration optimum and the trap

Integrating the cascade to $t = 240$ min at 37 °C over initial dimer concentration:

| [D]₀ (µM) | Closed capsid | Free dimer | Partial shells | Mean partial size |
|---|---|---|---|---|
| 0.50 | 7.9% | 88.1% | 4.0% | 20 |
| 1.00 | 41.9% | 54.4% | 3.7% | 27 |
| 2.00 | 77.4% | 21.6% | 1.0% | 20 |
| **2.48** | **94.8%** | 2.4% | 2.8% | 76 |
| 3.00 | 52.5% | 1.1% | 46.4% | 65 |
| 5.00 | 0.1% | 0.7% | 99.2% | 48 |
| 10.00 | 0.0% | 0.3% | 99.7% | 29 |

**More dimer is not better.** The optimum is interior, at $[D]_0 = 2.48$ µM. The mechanism
of the collapse is visible in the columns: between 2.5 and 5 µM the yield falls from 95% to
nil *while free dimer stays below 1%*. The protein has not failed to react — it has all
reacted, into partial shells that can no longer find the subunits they need. That is
**kinetic trapping**, and it is the physically interesting content of the layer. Zlotnick
and Endres established the mechanism; Hagan and Chandler's simulations connect it to the
same rate-equation structure used here.

The existence of an interior optimum is also the model's principal falsifiable claim.
A monotonically saturating yield curve would break the detailed-balance construction and
both scales built on it.

---

## 3. Layer L2 — classical nucleation theory as the widening layer

L1 answers "how much, when" for one parameter set. L2 answers a different question — **for
which contact energies and concentrations is productive assembly possible at all** — by
coarse-graining the ladder into a continuous droplet free energy. The gain is a phase
diagram; the cost is that discreteness and geometry are lost.

### 3.1 The droplet free energy

Treat cluster size $n$ as continuous. Growth is driven by a per-subunit chemical potential
difference and opposed by the unsatisfied bonds at the cluster rim. Bulk scales as $n$,
rim as the perimeter, which for a two-dimensional cap on a shell scales as $n^{2/3}$ in the
CNT convention:

$$G(n) = -n\,\Delta\mu + \sigma a_1 n^{2/3}$$

with $\sigma$ an edge tension and $a_1$ a geometric prefactor. The drive combines the
contact energy with the entropic cost of removing a subunit from solution:

$$\Delta\mu = -g_b\,\bar{c} - RT \ln\!\left(\frac{C_0}{C}\right)$$

where $g_b < 0$ is the energy per contact, $\bar{c}$ the net contacts a subunit gains on
joining, and $C_0$ a reference concentration. The logarithm is the supersaturation: **below
a threshold concentration $\Delta\mu \leq 0$ and no cluster is stable at any size.**

### 3.2 Barrier and critical nucleus

$G(n)$ rises, peaks, and falls. Setting $dG/dn = 0$:

$$-\Delta\mu + \tfrac{2}{3}\sigma a_1 n^{-1/3} = 0 \;\Longrightarrow\; n^{*} = \left(\frac{2\sigma a_1}{3\Delta\mu}\right)^{3}$$

and substituting back,

$$G^{*} = \frac{4}{27}\,\frac{(\sigma a_1)^3}{\Delta\mu^{2}}$$

Both scalings are informative in themselves: the critical nucleus shrinks as the inverse
cube of the drive and the barrier as its inverse square. **Doubling supersaturation cuts the
barrier fourfold and the nucleus eightfold** — which is why capsid yield is so violently
concentration-sensitive, and why the optimum in L1 is sharp rather than broad.

### 3.3 Three regimes, and a phase diagram

Classifying each point of the $(g_b, C)$ plane by the barrier it produces:

| Regime | Condition |
|---|---|
| no assembly | $\Delta\mu \leq 0$, or $G^{*} > 35$ kcal/mol |
| kinetic trap | $G^{*} < 1$ kcal/mol |
| productive | otherwise |

The two thresholds are timescale statements, not new physics: a barrier above
$\sim 35$ kcal/mol is not crossed in an incubation, and a barrier below $\sim 1$ kcal/mol is
no barrier at all, so nucleation outruns elongation and every subunit starts a shell that
none can finish. Over $g_b \in [-1, -5]$ kcal/mol per contact and
$C \in [0.3, 300]$ µM, with $\bar{c} = 1.0$, $\sigma = 1.6$ kcal/mol, $a_1 = 1.5$,
$C_0 = 500$ µM:

- **no assembly 40%** of the plane,
- **productive 27%**,
- **kinetic trap 33%**.

A representative productive point, $g_b = -1.5$ kcal/mol and $C = 120$ µM, gives
$\Delta\mu = 0.65$ kcal/mol, $n^{*} = 14.7$ dimers and $G^{*} = 4.79$ kcal/mol $= 8.1\,k_BT$
— a barrier crossed often enough to nucleate, rarely enough not to exhaust the pool.

**The headline is that productive assembly occupies only about a quarter of the accessible
parameter plane, and is bracketed on both sides by failure.** Assembly conditions are not a
threshold to exceed but a window to hit. Note that L2's $n^{*} \approx 15$ is much larger
than L1's assumed $n_c = 4$; the two are not directly comparable, because L2's continuum
rim term overestimates the penalty for the very small clusters where discreteness dominates.
Reconciling them is exactly what L3 is for.

---

## 4. Layer L3 — coarse-grained Monte Carlo: measuring what L1 assumes

L1 and L2 both take the geometry as given: a contact count, an edge tension. L3 puts 90
subunits in a box with explicit positions, orientations and a directional interaction, and
lets the geometry emerge.

### 4.1 A patchy-particle potential with built-in curvature

Each subunit is a hard sphere of diameter $\sigma$ decorated with attractive patches in the
Kern-Frenkel form: the interaction is switched on only when both partners point at each
other within an angular tolerance $\delta$. For a pair at separation $r$ with unit vector
$\hat{d}$,

| Configuration of the pair | Pair energy u(i,j) |
|---|---|
| $r < \sigma$ (hard core) | $+\infty$ |
| $r > r_{well}$ (out of range) | $0$ |
| $\sigma \leq r \leq r_{well}$, patches aligned | $-\varepsilon\left(\frac{1}{2} + \frac{1}{2}\chi_{ij}\right)$ |
| $\sigma \leq r \leq r_{well}$, patches not aligned | $0$ |

"Patches aligned" means, with $\hat{p}$ the patch vectors in the laboratory frame,

$$\max_k\, \hat{p}_i^{(k)}\!\cdot\!\hat{d} > \cos\delta \quad \text{and} \quad \max_l\, \hat{p}_j^{(l)}\!\cdot\!(-\hat{d}) > \cos\delta$$

with $\delta = 28°$ and $r_{well} = 1.2\,\sigma$.

The factor $\chi_{ij}$ is what makes this a *capsid* model rather than a generic
self-assembling fluid. A flat sheet is the ground state of a naive patchy sphere; a closed
shell requires a preferred dihedral angle between neighbours. Rewarding anti-alignment of
the engaged patch vectors,

$$\chi_{ij} = \tfrac{1}{2}\left(1 + \hat{p}_i^{(k)}\!\cdot\!\left(-\hat{p}_j^{(l)}\right)\right)$$

builds the spontaneous curvature in energetically, so growing caps close instead of
extending. The patch geometry itself is laid out as a curved five-patch wedge, the local
environment of a T=3 lattice site.

### 4.2 Sampling

Configurations are sampled by Metropolis Monte Carlo in the canonical ensemble, with
translation and rotation moves accepted at

$$P_{acc} = \min\left[1, \exp\!\left(-\frac{\Delta U}{k_BT}\right)\right]$$

and a cell list restricting neighbour search to $r < r_{well}$, which is what makes 90
particles over 1200 sweeps tractable. Clusters are defined by a bond threshold at
$-0.4\,\varepsilon$, i.e. a pair counts as bonded only if it is well inside the attractive
well rather than grazing it.

### 4.3 What L3 measures

Scanning volume fraction $\phi$ at $\varepsilon = 11$ $k_BT$:

| Volume fraction φ | Largest cluster (subunits) | Mean coordination | Well-formed fraction |
|---|---|---|---|
| 0.03 | 12.7 | 1.75 | 0.226 |
| 0.06 | 23.0 | 1.89 | 0.244 |
| 0.10 | 35.7 | 2.13 | 0.393 |
| 0.15 | 58.3 | 2.07 | 0.393 |
| 0.22 | 47.0 | 1.75 | 0.304 |
| 0.30 | 38.3 | 1.38 | 0.211 |

**The optimum in $\phi$ is interior — near 0.10–0.15 — and it is not the same statement as
L1's optimum in concentration, but it is the same physics.** Crowding first accelerates
nucleation, then produces more nuclei than there are subunits to complete them, and mean
coordination *falls* as the largest cluster shrinks: the signature of malformed aggregates
rather than shells.

L3 is where the assumption of L1 becomes a measurement, which is the subject of the next
section.

---

## 5. The bridge L3 → L1, and the disagreement we chose to publish

### 5.1 A one-parameter cap law

The bonds in a compact cap of $j$ subunits on a closed shell of $N$ can be written as a bulk
term minus the unsatisfied open perimeter. In a T=3 shell each dimer has five neighbours,
so a closed shell holds $\tfrac{5}{2}(N-1)$ bonds per the handshake count; the rim
correction must vanish at both $j = 1$ and $j = N$, which the geometric mean of occupied and
vacant fractions does. Hence

$$b(j) = \tfrac{5}{2}(j-1) - a\left[\sqrt{j\left(1 - \tfrac{j}{N}\right)} - s_1\right]$$

with $s_1$ chosen so $b(1) = 0$ exactly and a **single** shape parameter $a$ carrying the
rim penalty. Fitting $a$ to 2,547 clusters censused from the L3 trajectories at
$\varepsilon = 6, 7, 8$:

$$a = 6.75 \pm 0.37 \quad (\text{RMS residual } 4.8 \text{ bonds})$$

and the law extrapolates to $b(90) = 229$ against the geometric 225 — a 1.8% closure error
from a fit that never saw a complete shell. The differentiated ladder

$$c(j) = b(j+1) - b(j)$$

is the measured analogue of L1's step function. **It is not a step function.** It rises
continuously, from $c \approx 0.7$ at $j = 4$, through $c \approx 1.5$ at $j = 10$, to the
bulk value at large $j$.

### 5.2 An emergent critical nucleus

With a measured ladder, $n_c$ stops being an input. It is the argmax of the cumulative
profile:

$$n_c = \arg\max_j \left[\, j\, g_{assoc} + b(j)\, g_{contact} \,\right]$$

Substituting the measured $a$ at *unchanged* literature energies abolishes assembly
outright: the barrier becomes 24 kcal/mol $\approx 41\,k_BT$ and nothing nucleates in 240
minutes. That is not a bug but the layer's most instructive output — **the ladder shape and
the energy offsets are not independently transferable**, because $g_{assoc}$ was originally
calibrated against the step ladder.

Enforcing equal barrier height instead requires $g_{assoc} = 4.01$ kcal/mol — and at that
value the CG geometry independently returns

$$n_c = 4$$

the very value L1 assumed. A genuine cross-scale consistency check, on a quantity neither
layer was tuned to reproduce, and it passes. It also agrees with the 3–5 dimer nucleus
inferred from single-particle waiting times by Garmann and colleagues.

### 5.3 The disagreement

Matching barrier *height* is not matching ladder *shape*, and the shape is what sets the
concentration optimum:

| Model | Predicted optimal [D]₀ | Yield at optimum |
|---|---|---|
| Step ladder (production L1) | 2.50 µM | 95% |
| CG-measured ladder, $g_{assoc} \approx 4.0$ | 9.25 µM | 97% |

**The two scales disagree by 3.7× on the single number the wet lab will set.** Across the
literature box for the energies the CG-bridged optimum spans 7.8–15 µM and never falls
below 5.8 µM. Inverting the other way — demanding the step-ladder optimum of 2.5 µM from
the measured geometry — requires $g_{contact} = -5.55$ kcal/mol, outside Ceres and
Zlotnick's independently measured $-4.4$ to $-2.9$ interval.

Three sources of evidence cannot all be right. We publish the disagreement rather than
picking the flattering branch, because one benchtop concentration sweep adjudicates it:
a peak near 8–12 µM vindicates the CG geometry, near 2.5 µM indicts it, and no peak at all
falsifies the detailed-balance construction underneath both.

---

## 6. Layer L4 — packaging-signal-mediated co-assembly

L1–L3 model protein alone. The BruSafe particle assembles *around an mRNA carrying
engineered pac-sites*, and Dykeman, Stockley and Twarock showed that such signals change
assembly qualitatively, not just quantitatively. L4 is an independent re-implementation of
their framework on a 12-face study system, validated on exact invariants before any kinetics
were compared.

### 6.1 Two reaction orders, both at detailed balance

Two classes of event occur. A subunit can bind a packaging signal on the cargo — second
order, since it consumes free protein — or two adjacent bound subunits can form a
protein-protein bond — first order, an intramolecular closure. Detailed balance is imposed
separately on each. For a signal of affinity $\Delta G_R$ at cargo concentration $C$,

$$k_{2R} = \frac{k_{1R}}{C\,\exp\!\left(-\Delta G_R / RT\right)}$$

and for a protein-protein bond of energy $\Delta G_P$,

$$k_{2P} = k_{1P}\exp\!\left(\frac{\Delta G_P}{RT}\right)$$

Because subunit count scales with ensemble size at fixed volume, the on-rate must be
rescaled against the reference ensemble to hold concentration constant:

$$k_{1R} \to k_{1R}\,\frac{N_{RNA}}{N_{ref}}, \qquad N_{ref} = 3000$$

Events are drawn by the Gillespie algorithm — with total propensity $A$, the waiting time is

$$\tau = -\frac{\ln u}{A}, \qquad u \sim \mathcal{U}(0,1)$$

and propensities are held in a binary-indexed tree so that selecting an event costs
$O(\log N)$ rather than $O(N)$.

### 6.2 The trap is geometric, not energetic

The mechanism that makes this layer non-trivial: the cargo **threads the shell as a
self-avoiding walk on the face-adjacency graph**. A growing end can reach a face all of
whose neighbours are already visited — a dead end. The shell is then incomplete and cannot
be finished without reversing, and reversal means paying $k_{2P}$ to break protein bonds.
The trap therefore arises from the topology of the threading path, not from an energy
barrier. Two earlier versions of this engine were discarded for getting exactly this wrong:
one allowed free rerouting and produced near-complete yield in every condition (nothing ever
trapped), the other gated escape on an invented reserved-face rule and over-trapped.

### 6.3 Validation on exact invariants

Before comparing any kinetics, two combinatorial quantities of the published model were
reproduced exactly:

- **73** unique protein-only intermediates,
- **1,264** distinct threading paths, from 75,840 directed Hamiltonian paths divided by the
  order-60 icosahedral rotation group.

Two exact agreements on quantities that depend only on the geometry and the assembly rules
are the licence to compare rates at all.

### 6.4 What L4 predicts

| Finding | Numbers |
|---|---|
| Interior optimum in signal affinity, 12-PS cargo | 96.5% at $\Delta G_{PS} = -2.0$, falling to 54.3% at $-12$ |
| **Optimum reverses for a sparse 3-pac-site cargo** | 44.5% at $-2.0$, rising to **60.1% at $-9.0$** |
| Signal count is the largest unused lever | 20.3 / 44.5 / 60.8 / 96.5% at 1 / 3 / 6 / 12 sites (≈14 points per site at $n=3$) |
| Layout matters at matched mean affinity | 300 random cargoes span 63.3–97.5% (mean 80.2, sd 6.1) |
| Protein-protein window is narrow | 95.1% at $\Delta G_P = -2.0$; 42.8% at $-3.0$; 18.3% at $-4.0$ |

The reversal is the result that matters for BruSafe and is **not present in the source
paper**, which only ever simulated fully decorated genomes. For a densely decorated cargo,
very tight signals are harmful: they nucleate everywhere and freeze partial shells. For a
sparse cargo — three pac-sites, as in our cassette — there is no competing nucleation to
suppress, and tighter is better. The MS2 TR hairpin itself sits at $\Delta G = RT \ln K_d
= -11.6$ to $-12.3$ kcal/mol for $K_d = 1$–3 nM, i.e. at the tight end. **The theory
therefore supports keeping the wild-type recognition core rather than weakening it**, which
is what the de-repeated triple-pac-site design does.

What L4 does *not* reproduce is stated for the record: the paper's absolute yields (their
dimensionless conversion factor is described only as adjusted, never given, so absolute
yields are underdetermined by the published text), their $-3.0$ kcal/mol optimum, and their
nucleator-position preference. These are flagged, not tuned.

---

## 7. The parameter registry

Every quantity Module 4 consumes, its value, its literature bracket, and the measurement
that would pin it down. No value in this table is a free fit.

| Parameter | Value | Literature range | Identified by | Source |
|---|---|---|---|---|
| $g_{contact}$ | −4.0 kcal/mol | −4.4 to −2.9 | Position of the yield optimum | Ceres & Zlotnick 2002 |
| $g_{assoc}$ | +6.0 kcal/mol | 4 to 8 | Held fixed; absorbed into $g_{contact}$ if floated | Zlotnick 1994 |
| $k_a$ | 3.0 µM⁻¹min⁻¹ | 0.4 to 4 | Peak height; nucleation lag | Hagan & Chandler 2010 |
| $n_c$ | 4 dimers | 3 to 6 | Single-particle waiting times | Garmann 2022 |
| $N$ | 90 dimers | structural, not fitted | Crystallography | Zlotnick 1994; Ni 1995 |
| $K_d^{dim}$ | 0.01 µM | 0.001 to 1 | Analytical ultracentrifugation | Ni 1995 |
| $\sigma$ (L2 edge tension) | 1.6 kcal/mol | calibrated | Barrier height | derived, this work |
| $\bar{c}$ (L2 net contacts) | 1.0 | calibrated | Drive slope | derived, this work |
| $a$ (cap-law shape) | 6.75 ± 0.37 | measured | CG bond census, 2547 clusters | derived, this work |
| $\Delta G_P$ (L4) | −2.5 kcal/mol | −2.0 to −4.0 | Yield vs protein-protein energy | Dykeman 2013 |
| $\Delta G_{PS}$ (L4) | −2.0 to −12 kcal/mol | measured for TR | Filter binding, ITC | Dykeman 2013; Stockley 2016 |
| $RT$ | 0.593 kcal/mol | 310 K | — | — |

---

## 8. Limits of validity

Stated plainly, because a model whose limits are hidden cannot be trusted where it does work.

1. **L1 is a mean-field rate equation.** It tracks concentrations of size classes, not
   individual shells, so it cannot distinguish a well-formed 90-mer from a malformed one of
   the same size. Its yield is properly read as an **upper bound on well-formed particles**.
   Hagan's review of assembly modelling makes this limitation explicit for exactly this
   class of model.
2. **One-subunit-at-a-time growth** excludes fusion of partial shells. Justified while
   intermediates are rare; it is least safe deep in the trapped regime, which is where the
   model's most dramatic prediction lives.
3. **L2's continuum rim term** is wrong at small $n$, where discreteness dominates. Its
   $n^{*} \approx 15$ should be read as a scaling statement, not a prediction of the nucleus.
4. **L3 is a 90-particle model in reduced units**, with a patch geometry that approximates a
   T=3 lattice site rather than reproducing the MS2 interface. Its $\phi$ optimum transfers
   as a direction, not a molar concentration.
5. **L4 runs on a 12-face study system**, not the T=3 lattice. Its ordering and reversals
   transfer; its magnitudes do not. Porting it to 90 faces is the natural extension.
6. **No electrostatics.** Cargo charge, ionic strength and pH enter only implicitly through
   effective energies, whereas Chevreuil and colleagues show the assembly pathway itself can
   depend on them.
7. **Temperature enters only through $RT$.** The hydrophobic, entropy-driven character of the
   interface means the energies themselves are temperature-dependent; the temperature lever
   is capped accordingly.
8. **The two false retraction flags** in the corpus (`10.1371/journal.pcbi.1010652`,
   `10.1371/journal.pcbi.1009913`) are `new_version` updates re-verified against Crossref,
   not retractions, and must be cleared before the corpus is cited as evidence of screening.

---

## 9. Wet-lab validation: where this model meets the bench

Everything above is Mode A — prediction with no data in it. This section is the contract
between the model and the wet lab: which measurements enter, which parameter each one moves,
and what the model is obliged to do when the number comes back wrong. It is written for the
equipment actually available: a plate reader, a gel rig, a spectrophotometer, and **at most
three TEM sessions**, which is the binding constraint on the whole design.

### 9.1 The one experiment that matters most

**A concentration sweep of coat-protein dimer, read by turbidity, with TEM at three points.**

Every layer of Module 4 makes the same qualitative claim — that yield is **non-monotonic**
in dimer concentration, rising to an optimum and then collapsing into kinetic traps. L1
derives it from rate equations, L2 explains it as a moving barrier, L3 reproduces it from
geometry alone. A single sweep tests all three at once. If yield instead saturates
monotonically, the central claim of the module is wrong, and no amount of re-fitting rescues
it.

Design: assemble at $[D]_0 \in \{0.5, 1, 2, 3, 5\}$ µM (add 8 and 12 µM if protein allows),
240 min, three replicates, one buffer. Read $A_{340}$ continuously in the plate reader.

### 9.2 Why turbidity alone is not enough

Turbidity is a mass-squared signal: for particles small against the wavelength, scattering
per species goes as (mass)$^2$, so the plate reader reports

$$\tau(t) \;\propto\; \sum_j j^{2}\,x_j(t),$$

which a **partial shell answers almost as loudly as a closed one**. A 60-dimer malformed
aggregate scatters $(60/90)^2 = 44\%$ as much as a complete capsid. This is precisely the
kinetic-trap signature: above the optimum the model predicts turbidity that *plateaus high*
while the closed-capsid count *collapses to zero*. A plate reader on its own cannot tell that
apart from success, which is the entire reason TEM sessions must be spent, and spent in the
right places.

Native PAGE is the cheap partial fix: closed capsids run as a single sharp species, trapped
intermediates as a smear. It resolves the ambiguity qualitatively for every concentration in
the sweep, at no TEM cost.

### 9.3 Spending three TEM sessions

The three sessions are the scarce resource, so they are allocated to the three points where
the model's predictions differ most, not spread evenly across the sweep:

| Session | $[D]_0$ | Model predicts | What it discriminates |
|---|---|---|---|
| 1 | at the **turbidity peak** | mostly closed 27 nm shells | confirms the optimum is productive, not just scattering |
| 2 | **2× above the peak** | high turbidity, few closed shells, mass in large partial shells | the trap prediction — the single most falsifiable claim in the module |
| 3 | **~0.4× the peak** | few particles, most mass free dimer | confirms the low side fails by *not nucleating*, not by trapping |

Session 2 is the one that earns the module its keep. Do not sacrifice it to a replicate.

Count ≥ 200 particles per grid and score each as correct-diameter (25–29 nm) or not; the
fraction is what the model consumes. Record the count itself — it sets the binomial error
bar, and a fraction with no $n$ behind it cannot be weighted in a fit.

### 9.4 What to graph, and what each graph is for

1. **Turbidity traces, $A_{340}$ vs time, one curve per concentration.** Read the *lag* and
   the *plateau*, not just the endpoint. Lag length is a nucleation observable; plateau height
   is a mass observable. Generated by the bridge notebook as `M4_bridge_observables.png`.
2. **Yield curve: capsid fraction vs $[D]_0$, log-x.** The headline figure. The model's
   optimum position and height are read straight off it. `M4_L1_yield_curve.png`.
3. **TEM size histogram at the three sessions, mass-weighted, log-y.** The trap is visible
   as a mid-size hump that has no business being there. A linear y-axis hides it.
4. **Prior versus posterior yield curve with your measurement on it.** This is the figure for
   the wiki: it shows the model *moving* because of your data.
   `M4_bridge_correction.png` and `M4_FIT_prior_vs_posterior.png`.
5. **Profile-likelihood curve for the contact energy.** Turns a fitted number into an interval,
   and shows a judge that the fit was interrogated rather than declared. `M4_FIT_profile.png`.

### 9.5 How the model is adjusted afterwards — which observable moves which parameter

The fit is identifiable only because different observables carry different parameters. This
was checked before any data existed, by scanning each parameter across its literature range
and recording what moved:

| Observable | Parameter it constrains | Evidence it is the right pairing |
|---|---|---|
| **Position** of the yield optimum | $g_{\text{contact}}$ | across its literature range the optimum moves 12.0 → 0.44 µM |
| **Height** of the optimum, and lag length | $k_a$ | across its range the position stays at 2.52 µM while height moves 0.39 → 0.95 |
| Overall vertical offset | $s_{\text{rec}}$, a recovery scale | absorbs purification and grid-sampling losses that are not physics |
| TEM size distribution above the optimum | validates $n_c$ and the ladder | not fitted — a genuine out-of-sample test |

Because position and height are carried by different parameters, two numbers — the optimum's
location `D_MEAS` and its height `Y_MEAS` — invert to two energies. That inversion is the
corrective step in `M4_bridge_CG_to_ODE_to_bench.ipynb`; the fuller weighted fit over the whole
curve, with the profile interval, is `M4_FIT_to_bench_data.ipynb`.

**$n_c$ is deliberately held fixed at 4 in the fit.** From concentration data alone it is
confounded with $g_{\text{contact}}$ — a larger nucleus and a weaker contact shift the yield
curve the same way — so floating both would produce a tighter-looking fit that constrains
neither. The value 4 is the literature prior, not a measurement.

The bridge measures $n_c$ independently, from where the coarse-grained contact ladder turns
the free-energy barrier over, and reads $n_c = 8$. **That factor of two is a real
disagreement between the two routes, and we report it as one.** It is not reconciled by
tuning either number: the prior is a mean-field estimate for a $T{=}3$ shell, the bridge
value comes from one patch geometry at one bond energy, and neither has yet been checked
against bench data. What the disagreement does mean is that any yield-curve fit quoting a
contact energy must also quote the $n_c$ it assumed — which is why the profile interval in
`M4_FIT_to_bench_data.ipynb` is reported at fixed $n_c$ rather than marginalised over it.

### 9.6 What we do if the data disagrees

Stated in advance, so the response is not chosen after seeing the result:

- **Optimum in the wrong place** → re-fit $g_{\text{contact}}$. Expected; this is what the
  fitting notebook is for.
- **No optimum at all, monotonic saturation** → L1's mechanism is wrong for MS2. Report it.
  The nucleation-elongation picture would not survive that, and dressing it up would be worse
  than losing it.
- **Optimum present but TEM above it shows free dimer rather than partial shells** → the
  yield collapse is real but the *mechanism* is not trapping. L1 keeps its curve and loses its
  explanation.
- **Fitted $g_{\text{contact}}$ outside the −2.9 to −4.4 kcal/mol literature interval** → this
  has already happened. The bridge inverts to **−4.78 kcal/mol**, and the notebook prints the
  comparison as a plain `False`. It is reported as a disagreement (§5.3), not tuned away.

### 9.7 Data format

`M4_bench_data_template.csv` carries the column names both data-facing notebooks read:
`assay, D0_uM, t_min, replicate, value, value_sd, units, n_particles, notes`, with `assay` one
of `turbidity`, `nativePAGE`, `TEM`. Fill it in, point the fitting notebook at it, and every
figure in §9.4 regenerates. `M4_bench_data_example.csv` is synthetic and exists only so the
notebooks run before data exists; it must be replaced, not fitted.

---

## 10. What would falsify each layer

| Layer | Prediction | Falsified by |
|---|---|---|
| L1 | Interior optimum in $[D]_0$, with collapse at high concentration and free dimer staying low | A monotonically saturating yield curve |
| L1 | Trapped mass sits in large partial shells (78–79 dimers at 2.5 µM, 69–71 at 3.0, 48–50 at 5.0) | TEM showing free dimer, not partial shells, above the optimum |
| L2 | Productive assembly occupies a bounded window, not a threshold | Productive yield at arbitrarily high concentration |
| L3 → L1 | Optimum at 8–12 µM | Optimum at 2.5 µM (indicts CG geometry) or no peak (breaks both) |
| L4 | Tighter pac-sites help a sparse cargo | Yield falling monotonically with pac-site affinity at $n_{PS} = 3$ |

The concentration sweep is the single highest-information experiment in Module 4: turbidity
across 1–20 µM dimer, with TEM at three points, because a turbidity plateau *without* a
matching capsid count is the trap signature and turbidity alone cannot see it. The full design is §9.

---

## 11. Bibliography

All entries below are Crossref-verified; the full 70-record corpus with themes and
verification columns is `M4_corpus_verified.csv`.

- Asor et al. (2023). *Effect of ionic strength on the assembly of simian vacuolating virus capsid protein around poly(styrene sulfonate)*. The European Physical Journal E. [10.1140/epje/s10189-023-00363-x](https://doi.org/10.1140/epje/s10189-023-00363-x)
- Baschek et al. (2012). *Stochastic dynamics of virus capsid formation: direct versus hierarchical self-assembly*. BMC Biophysics. [10.1186/2046-1682-5-22](https://doi.org/10.1186/2046-1682-5-22)
- Borodavka et al. (2012). *Evidence that viral RNAs have evolved for efficient, two-stage packaging*. Proceedings of the National Academy of Sciences. [10.1073/pnas.1204357109](https://doi.org/10.1073/pnas.1204357109)
- Caspar et al. (1962). *Physical Principles in the Construction of Regular Viruses*. Cold Spring Harbor Symposia on Quantitative Biology. [10.1101/sqb.1962.027.001.005](https://doi.org/10.1101/sqb.1962.027.001.005)
- Ceres et al. (2002). *Weak Protein−Protein Interactions Are Sufficient To Drive Assembly of Hepatitis B Virus Capsids*. Biochemistry. [10.1021/bi0261645](https://doi.org/10.1021/bi0261645)
- Chevreuil et al. (2018). *Nonequilibrium self-assembly dynamics of icosahedral viral capsids packaging genome or polyelectrolyte*. Nature Communications. [10.1038/s41467-018-05426-8](https://doi.org/10.1038/s41467-018-05426-8)
- Depta et al. (2022). *Hierarchical Coarse-Grained Strategy for Macromolecular Self-Assembly: Application to Hepatitis B Virus-Like Particles*. International Journal of Molecular Sciences. [10.3390/ijms232314699](https://doi.org/10.3390/ijms232314699)
- Dykeman et al. (2011). *Simple Rules for Efficient Assembly Predict the Layout of a Packaged Viral RNA*. Journal of Molecular Biology. [10.1016/j.jmb.2011.02.039](https://doi.org/10.1016/j.jmb.2011.02.039)
- Endres et al. (2002). *Model-Based Analysis of Assembly Kinetics for Virus Capsids or Other Spherical Polymers*. Biophysical Journal. [10.1016/s0006-3495(02)75245-4](https://doi.org/10.1016/s0006-3495(02)75245-4)
- Ford et al. (2013). *Sequence-Specific, RNA–Protein Interactions Overcome Electrostatic Barriers Preventing Assembly of Satellite Tobacco Necrosis Virus Coat Protein*. Journal of Molecular Biology. [10.1016/j.jmb.2013.01.004](https://doi.org/10.1016/j.jmb.2013.01.004)
- Frechette et al. (2025). *Computer Simulations Show That Liquid–Liquid Phase Separation Enhances Self-Assembly*. ACS Nano. [10.1021/acsnano.5c08120](https://doi.org/10.1021/acsnano.5c08120)
- Garmann et al. (2022). *Single-particle studies of the effects of RNA–protein interactions on the self-assembly of RNA virus particles*. Proceedings of the National Academy of Sciences. [10.1073/pnas.2206292119](https://doi.org/10.1073/pnas.2206292119)
- Gartner et al. (2024). *Design Principles for Fast and Efficient Self-Assembly Processes*. Physical Review X. [10.1103/physrevx.14.021004](https://doi.org/10.1103/physrevx.14.021004)
- Guo et al. (2024). *Modeling Viral Capsid Assembly: A Review of Computational Strategies and Applications*. Cells. [10.3390/cells13242088](https://doi.org/10.3390/cells13242088)
- Hagan et al. (2009). *A theory for viral capsid assembly around electrostatic cores*. The Journal of Chemical Physics. [10.1063/1.3086041](https://doi.org/10.1063/1.3086041)
- Hagan et al. (2010). *Understanding the Concentration Dependence of Viral Capsid Assembly Kinetics—the Origin of the Lag Time and Identifying the Critical Nucleus Size*. Biophysical Journal. [10.1016/j.bpj.2009.11.023](https://doi.org/10.1016/j.bpj.2009.11.023)
- Hagan et al. (2011). *Mechanisms of kinetic trapping in self-assembly and phase transformation*. The Journal of Chemical Physics. [10.1063/1.3635775](https://doi.org/10.1063/1.3635775)
- Hagan et al. (2023). *Self-assembly coupled to liquid-liquid phase separation*. PLOS Computational Biology. [10.1371/journal.pcbi.1010652](https://doi.org/10.1371/journal.pcbi.1010652)
- Hamilton et al. (2024). *RNA-induced Allosteric Coupling Drives Viral Capsid Assembly*. PRX Life. [10.1103/prxlife.2.013012](https://doi.org/10.1103/prxlife.2.013012)
- Hashemi et al. (2021). *Optimizing the synthesis and purification of MS2 virus like particles*. Scientific Reports. [10.1038/s41598-021-98706-1](https://doi.org/10.1038/s41598-021-98706-1)
- Hemberg et al. (2006). *Stochastic Kinetics of Viral Capsid Assembly Based on Detailed Protein Structures*. Biophysical Journal. [10.1529/biophysj.105.076737](https://doi.org/10.1529/biophysj.105.076737)
- Hesketh et al. (2015). *Mechanisms of assembly and genome packaging in an RNA virus revealed by high-resolution cryo-EM*. Nature Communications. [10.1038/ncomms10113](https://doi.org/10.1038/ncomms10113)
- Jhaveri et al. (2024). *Discovering optimal kinetic pathways for self-assembly using automatic differentiation*. Proceedings of the National Academy of Sciences. [10.1073/pnas.2403384121](https://doi.org/10.1073/pnas.2403384121)
- Kivenson et al. (2010). *Mechanisms of Capsid Assembly around a Polymer*. Biophysical Journal. [10.1016/j.bpj.2010.04.035](https://doi.org/10.1016/j.bpj.2010.04.035)
- Kra et al. (2023). *Energetics and Kinetic Assembly Pathways of Hepatitis B Virus Capsids in the Presence of Antivirals*. ACS Nano. [10.1021/acsnano.3c03595](https://doi.org/10.1021/acsnano.3c03595)
- Lago et al. (2001). *Probing the kinetics of formation of the bacteriophage MS2 translational operator complex: identification of a protein conformer unable to bind RNA*. Journal of Molecular Biology. [10.1006/jmbi.2000.4355](https://doi.org/10.1006/jmbi.2000.4355)
- Le et al. (2021). *In Vitro Assembly of Virus-Like Particles and Their Applications*. Life. [10.3390/life11040334](https://doi.org/10.3390/life11040334)
- Lino et al. (2017). *Display of single-chain variable fragments on bacteriophage MS2 virus-like particles*. Journal of Nanobiotechnology. [10.1186/s12951-016-0240-7](https://doi.org/10.1186/s12951-016-0240-7)
- Mendoza et al. (2020). *Shape selection and mis-assembly in viral capsid formation by elastic frustration*. eLife. [10.7554/elife.52525](https://doi.org/10.7554/elife.52525)
- Michaels et al. (2017). *Kinetic constraints on self-assembly into closed supramolecular structures*. Scientific Reports. [10.1038/s41598-017-12528-8](https://doi.org/10.1038/s41598-017-12528-8)
- Mizrahi et al. (2022). *Packaging contests between viral RNA molecules and kinetic selectivity*. PLOS Computational Biology. [10.1371/journal.pcbi.1009913](https://doi.org/10.1371/journal.pcbi.1009913)
- Mohajerani et al. (2022). *Multiscale Modeling of Hepatitis B Virus Capsid Assembly and Its Dimorphism*. ACS Nano. [10.1021/acsnano.2c02119](https://doi.org/10.1021/acsnano.2c02119)
- Naskalska et al. (2024). *Virus-like particles derived from bacteriophage MS2 as antigen scaffolds and RNA protective shells*. Nanomedicine. [10.2217/nnm-2023-0362](https://doi.org/10.2217/nnm-2023-0362)
- Nguyen et al. (2005). *Elasticity theory and shape transitions of viral shells*. Physical Review E. [10.1103/physreve.72.051923](https://doi.org/10.1103/physreve.72.051923)
- Ni et al. (1995). *Crystal structure of the MS2 coat protein dimer: implications for RNA binding and virus assembly*. Structure. [10.1016/s0969-2126(01)00156-3](https://doi.org/10.1016/s0969-2126(01)00156-3)
- Nooraei et al. (2021). *Virus-like particles: preparation, immunogenicity and their roles as nanovaccines and drug nanocarriers*. Journal of Nanobiotechnology. [10.1186/s12951-021-00806-7](https://doi.org/10.1186/s12951-021-00806-7)
- Oliver et al. (2020). *Assembly of Capsids from Hepatitis B Virus Core Protein Progresses through Highly Populated Intermediates in the Presence and Absence of RNA*. ACS Nano. [10.1021/acsnano.0c03569](https://doi.org/10.1021/acsnano.0c03569)
- Pan et al. (2012). *Development of a microRNA delivery system based on bacteriophage MS2 virus‐like particles*. The FEBS Journal. [10.1111/j.1742-4658.2012.08512.x](https://doi.org/10.1111/j.1742-4658.2012.08512.x)
- Peabody et al. (1999). *Asymmetric Contributions to RNA Binding by the Thr45Residues of the MS2 Coat Protein Dimer*. Journal of Biological Chemistry. [10.1074/jbc.274.36.25403](https://doi.org/10.1074/jbc.274.36.25403)
- Peabody et al. (2008). *Immunogenic Display of Diverse Peptides on Virus-like Particles of RNA Phage MS2*. Journal of Molecular Biology. [10.1016/j.jmb.2008.04.049](https://doi.org/10.1016/j.jmb.2008.04.049)
- Perkett et al. (2014). *Using Markov state models to study self-assembly*. The Journal of Chemical Physics. [10.1063/1.4878494](https://doi.org/10.1063/1.4878494)
- Perlmutter et al. (2015). *Mechanisms of Virus Assembly*. Annual Review of Physical Chemistry. [10.1146/annurev-physchem-040214-121637](https://doi.org/10.1146/annurev-physchem-040214-121637)
- Pickett et al. (1993). *Encapsidation of heterologous RNAs by bacteriophage MS2 coat protein*. Nucleic Acids Research. [10.1093/nar/21.19.4621](https://doi.org/10.1093/nar/21.19.4621)
- Plevka et al. (2009). *Structure and stability of icosahedral particles of a covalent coat protein dimer of bacteriophage MS2*. Protein Science. [10.1002/pro.184](https://doi.org/10.1002/pro.184)
- Polles et al. (2013). *Mechanical and Assembly Units of Viral Capsids Identified via Quasi-Rigid Domain Decomposition*. PLoS Computational Biology. [10.1371/journal.pcbi.1003331](https://doi.org/10.1371/journal.pcbi.1003331)
- Powell et al. (2001). *Asymmetric interactions in the adenosine-binding pockets of the MS2 coat protein dimer*. BMC Molecular Biology. [10.1186/1471-2199-2-6](https://doi.org/10.1186/1471-2199-2-6)
- Prel et al. (2015). *Highly efficient in vitro and in vivo delivery of functional RNAs using new versatile MS2-chimeric retrovirus-like particles*. Molecular Therapy - Methods &amp; Clinical Development. [10.1038/mtm.2015.39](https://doi.org/10.1038/mtm.2015.39)
- Pumpens et al. (2016). *The True Story and Advantages of RNA Phage Capsids as Nanotools*. Intervirology. [10.1159/000449503](https://doi.org/10.1159/000449503)
- Qian et al. (2023). *Temporal control by cofactors prevents kinetic trapping in retroviral Gag lattice assembly*. Biophysical Journal. [10.1016/j.bpj.2023.06.021](https://doi.org/10.1016/j.bpj.2023.06.021)
- Roldão et al. (2012). *On the Effect of Thermodynamic Equilibrium on the Assembly Efficiency of Complex Multi-Layered Virus-Like Particles (VLP): the Case of Rotavirus VLP*. PLoS Computational Biology. [10.1371/journal.pcbi.1002367](https://doi.org/10.1371/journal.pcbi.1002367)
- Rolfsson et al. (2008). *RNA Packing Specificity and Folding during Assembly of the Bacteriophage MS2*. Computational and Mathematical Methods in Medicine. [10.1080/17486700802168445](https://doi.org/10.1080/17486700802168445)
- Rolfsson et al. (2016). *Direct Evidence for Packaging Signal-Mediated Assembly of Bacteriophage MS2*. Journal of Molecular Biology. [10.1016/j.jmb.2015.11.014](https://doi.org/10.1016/j.jmb.2015.11.014)
- Schaaf et al. (1999). *A Small 2′-OH- and Base-dependent Recognition Element Downstream of the Initiation Site in the RNA Encapsidation Signal Is Essential for Hepatitis B Virus Replication Initiation*. Journal of Biological Chemistry. [10.1074/jbc.274.53.37787](https://doi.org/10.1074/jbc.274.53.37787)
- Singh et al. (2003). *Observed Hysteresis of Virus Capsid Disassembly Is Implicit in Kinetic Models of Assembly*. Journal of Biological Chemistry. [10.1074/jbc.m211408200](https://doi.org/10.1074/jbc.m211408200)
- Stockley et al. (2013). *Packaging signals in single-stranded RNA viruses: nature’s alternative to a purely electrostatic assembly mechanism*. Journal of Biological Physics. [10.1007/s10867-013-9313-0](https://doi.org/10.1007/s10867-013-9313-0)
- Tan et al. (2015). *The Interface between Hepatitis B Virus Capsid Proteins Affects Self-Assembly, Pregenomic RNA Packaging, and Reverse Transcription*. Journal of Virology. [10.1128/jvi.03545-14](https://doi.org/10.1128/jvi.03545-14)
- Twarock et al. (2019). *RNA-Mediated Virus Assembly: Mechanisms and Consequences for Viral Evolution and Therapy*. Annual Review of Biophysics. [10.1146/annurev-biophys-052118-115611](https://doi.org/10.1146/annurev-biophys-052118-115611)
- Vishwakarma et al. (2024). *Deciphering the Thermal Stability of Bacteriophage MS2-Derived Virus-like Particle and Its Engineered Variant*. ACS Biomaterials Science &amp; Engineering. [10.1021/acsbiomaterials.4c00770](https://doi.org/10.1021/acsbiomaterials.4c00770)
- Wang et al. (2016). *Preparation and evaluation of MS2 bacteriophage-like particles packaging hepatitis E virus RNA*. FEMS Microbiology Letters. [10.1093/femsle/fnw221](https://doi.org/10.1093/femsle/fnw221)
- Wang et al. (2019). *A Multiscale Model for the Self-Assembly of Coat Proteins in Bacteriophage MS2*. Journal of Chemical Information and Modeling. [10.1021/acs.jcim.9b00514](https://doi.org/10.1021/acs.jcim.9b00514)
- Wei et al. (2024). *Hierarchical assembly is more robust than egalitarian assembly in synthetic capsids*. Proceedings of the National Academy of Sciences. [10.1073/pnas.2312775121](https://doi.org/10.1073/pnas.2312775121)
- Whitelam et al. (2015). *The Statistical Mechanics of Dynamic Pathways to Self-Assembly*. Annual Review of Physical Chemistry. [10.1146/annurev-physchem-040214-121215](https://doi.org/10.1146/annurev-physchem-040214-121215)
- Williams et al. (2024). *Effect of coat-protein concentration on the self-assembly of bacteriophage MS2 capsids around RNA*. Nanoscale. [10.1039/d3nr03292b](https://doi.org/10.1039/d3nr03292b)
- Yvinec et al. (2012). *First passage times in homogeneous nucleation and self-assembly*. The Journal of Chemical Physics. [10.1063/1.4772598](https://doi.org/10.1063/1.4772598)
- Zandi et al. (2004). *Origin of icosahedral symmetry in viruses*. Proceedings of the National Academy of Sciences. [10.1073/pnas.0405844101](https://doi.org/10.1073/pnas.0405844101)
- Zhang et al. (2006). *Simulation Study of the Contribution of Oligomer/Oligomer Binding to Capsid Assembly Kinetics*. Biophysical Journal. [10.1529/biophysj.105.072207](https://doi.org/10.1529/biophysj.105.072207)
- Zhou et al. (2025). *MS2 virus-like particles as a versatile platform for multi-disease vaccines: a review*. Frontiers in Immunology. [10.3389/fimmu.2025.1651594](https://doi.org/10.3389/fimmu.2025.1651594)
- Zlotnick et al. (1994). *To Build a Virus Capsid*. Journal of Molecular Biology. [10.1006/jmbi.1994.1473](https://doi.org/10.1006/jmbi.1994.1473)
- Zlotnick et al. (1999). *A Theoretical Model Successfully Identifies Features of Hepatitis B Virus Capsid Assembly*. Biochemistry. [10.1021/bi991611a](https://doi.org/10.1021/bi991611a)
