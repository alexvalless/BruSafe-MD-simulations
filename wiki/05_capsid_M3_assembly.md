# CAPSID M3 · Capsid assembly

Notebooks: `BruSafe_04_M3_capsid_assembly.ipynb`, `literature/BruSafe_04L_M3_literature_dG.ipynb`

## Model

The *T* = 3 MS2 capsid consists of 180 coat-protein subunits, which assemble from **90 dimers**:
60 A/B dimers and 30 C/C dimers, with 180 dimer–dimer contacts in the closed shell. M3 follows every
intermediate size with a nucleation–elongation cascade over dimers,

C<sub>1</sub> + C<sub>j</sub> ⇌ C<sub>j+1</sub>, j = 1 … 89,

with a closed shell at N = 90. Association is degenerate (several equivalent addition sites),
dissociation is not, and forward and reverse rates are tied by detailed balance. The kinetics and the
equilibrium distribution are therefore the same model. The equilibrium distribution follows from the
free-dimer concentration alone and is computed in logarithms, since C<sub>1</sub><sup>90</sup>
underflows double precision.

The notebook takes the soluble dimer pool (or the dimer flux) from M2 and returns `assembly_sink`
to M2, expressed in monomer units.

## Contact energy from a measured melting point

The contact energy ΔG<sub>contact</sub> is not a free parameter. A thermal denaturation midpoint of the
VLP is a disassembly midpoint, so at the melting temperature the assay concentration equals the
critical concentration c\*, and Eq. (4.8) of the framework document can be inverted for
ΔG<sub>contact</sub>. Pistono et al. (2024) report T<sub>m</sub> ≈ 66 °C for wild-type MS2 in a thermal
shift assay run at 5–20 µM monomer (2.5–10 µM dimer). The inversion gives

| Assay concentration (dimer) | ΔG<sub>contact</sub> (kcal mol⁻¹) |
|---|---|
| 2.5 µM | −4.05 |
| 5 µM (used in notebook 04) | −3.82 |
| 10 µM | −3.59 |

An earlier placeholder of −4.2 kcal mol⁻¹ had no source and is retained only as a comparison line in
the figures. The values for the mutant ladder are in `capsid/data/litdata/derived_M3_dG_from_Tm.csv`.

## Results

**Equilibrium.** At every contact energy considered, intermediates stay below about 1% of the mass at
equilibrium, so an assembly assay measures essentially a two-species mixture of free dimer and
capsid. c\* and the equilibrium capsid fraction (0.98–0.998 across the sensitivity scan) depend only
on the total contact count and are robust.

**Assembly is limited by kinetics and not by thermodynamics.** Every expression arm of M2 lies two to
four orders of magnitude above c\*, so a titration designed to bracket c\* would spend its whole range
on one side of it. Being above c\* is necessary but not sufficient: at 20 µM dimer, which is above
c\*, the cascade produces essentially no capsid in 18 h because of the nucleation lag.

![Nucleation lag](../capsid/figures/fig19_M3_nucleation_lag.png)

*Figure 1. Left: capsid fraction against time at three dimer concentrations; only the 500 µM bolus
nucleates within 18 h. Right: size distribution of intermediates at 500 µM at four times; the
distribution is bimodal, with free dimer and closed shells and very little in between.*

**The kinetic trap saturates.** In this cascade, stalled intermediates plateau at a few percent
and do not grow with concentration. A yield that is non-monotonic in concentration would require an
off-pathway or malformed species, which the cascade does not contain.

**Bolus and feed.** An in vitro reassembly is a bolus, whereas in the cell dimer arrives continuously
at the rate supplied by M2. At the contact energy of −3.82 kcal mol⁻¹, a constant feed of the same total
dimer gives a higher stalled fraction than a bolus (0.040 against 0.008) and a lower capsid fraction at
18 h (0.741 against 0.884). The reaction has not finished at harvest: 131 µM of free dimer remains
unconverted. A slow feed therefore does not reduce the stalled fraction in this model.

## Coupling to packaging

In MS2, the TR/*pac* hairpin binds the A/B dimer and stabilises the bent FG-loop conformer, so it
could lower the nucleation barrier. This was modelled as extra contact energy on the first addition
only, which is a qualitative coupling and not a calibrated one. In the fed, in vivo regime, the
free-dimer pool is far above what nucleation requires, and the capsid fraction at 18 h changes by
only 1.6 percentage points between 0×, 1× and 3× *pac* constructs. Comparing harvested VLP yield
between *pac* constructs would therefore not test the allosteric coupling. The coupling would be
visible only where nucleation is limiting, that is, in vitro at 20–50 µM purified dimer with and
without the TR hairpin, followed by SEC or light scattering.

## Contact schedule from the structure

The contact schedule (the number of new dimer–dimer interfaces formed when the k-th dimer joins the
cluster) and the number of sites at which the next dimer can be added were enumerated from the
asymmetric unit of PDB 2MS2 by `capsid/tools/enumerate_contacts.py`, and the results are shipped in
`capsid/data/structure/` (`contact_schedule.csv`, `contact_schedule.json`). The procedure is:

1. The 60 icosahedral operators (REMARK 350) are applied to chains A, B and C, giving 180 chains.
2. The chains are grouped into the 90 dimers of the model: 60 A/B dimers (each A chain with the B chain
   it shares the largest interface with) and 30 C/C dimers.
3. Two dimers are in contact if at least 5 pairs of residues have heavy atoms closer than 4.5 Å. This
   gives **180 interfaces**, 60 between two A/B dimers and 120 between an A/B and a C/C dimer; no
   interface joins two C/C dimers. The mean interface has 53 residue pairs (A/B–A/B) and 49 (A/B–C/C).
4. Growth paths on this graph of 90 dimers are generated 4000 times for each of three rules (the next
   dimer chosen at random, or greedily with the most or the fewest contacts), and the contacts formed
   at each step are recorded. Every path closes with exactly 180 interfaces.

| Step k | Contacts formed (random growth, mean) | Sites where the next dimer can be added (mean) |
|---|---|---|
| 2 | 1.00 | 4.0 |
| 10 | 1.43 | 10.9 |
| 30 | 1.69 | 19.8 |
| 50 | 1.87 | 22.5 |
| 70 | 2.28 | 17.9 |
| 80 | 2.85 | 10.8 |
| 90 | 4.00 | 1.0 |

The first addition forms one interface and the last forms four. Contacts per addition rise slowly,
from 1.0 to about 2.0 over the first 60 steps, and then steeply as the shell closes (mean 1.66 for steps 2–59
and 2.70 for steps 60–90). The number of sites where the next dimer can be added is 4 for the first
addition, peaks at 22.5 near k = 49 and falls to 1 for the last dimer. The placeholder of notebook 04 has
contacts rising linearly from 1 to 3 and a constant number of sites equal to 2, so the structure
changes the degeneracy by up to an order of magnitude and the contacts of the final steps by one.

![Contact schedule from 2MS2](../capsid/figures/fig29_M3_contact_schedule.png)

*Figure 2. Contact schedule and boundary sites of the T = 3 capsid enumerated from PDB 2MS2 (4.5 Å,
at least 5 residue pairs per interface), for random growth (mean ± SD over 4000 paths) and for two greedy
rules, compared with the placeholder of notebook 04. The greedy curves are not bounds of all possible
paths.*

The schedule is geometric: it counts residue pairs in contact in the crystal lattice, treats every
interface as equivalent, and does not say how the dimers actually assemble. The edge threshold (5
residue pairs) is a choice, and `enumerate_contacts.py` prints the number of interfaces for thresholds of
3, 5, 10 and 20 residue pairs.

## What depends on the schedule

**The notebook has not yet been rerun with this schedule.** The outputs of notebook 04 shown on this page
were computed with the placeholder contact schedule and degeneracy vector, which have the correct
total and are not counts from the structure. Varying both over nine combinations gives:

| Output | Range over the nine combinations | Status |
|---|---|---|
| c\* | identical | robust |
| equilibrium capsid fraction | 0.978 – 0.998 | robust |
| kinetic capsid fraction at 18 h | 0.003 – 0.894 | **not robust** |
| stalled fraction | 0.014 – 0.129 | **not robust** |

The nucleation lag, the stalled fraction and the intermediate size distribution should not be quoted
until the structure-derived schedule above replaces the placeholder in the notebook and the model is
rerun. The enumeration requires no bench work and is done; wiring it in is the remaining step.

## Sensitivity to the contact energy (notebook 04L)

Across the range of ΔG<sub>contact</sub> allowed by the assay concentration (−4.05 to −3.59
kcal mol⁻¹), the equilibrium capsid fraction changes by under 1 percentage point, whereas the kinetic
capsid fraction at 18 h changes by about 28 percentage points. Equilibrium quantities are therefore
insensitive to the uncertainty of the assay, and kinetic ones are not.

![Contact energy sweep](../capsid/figures/fig25_M3_literature_dG_sweep.png)

*Figure 3. M3 outputs across the contact energies obtained from the published melting point.*

## What a TEM count can establish

Closed shells and fragments past half closure (45 dimers or more) are counted as particles by TEM,
with a 27 ± 2 nm acceptance window. At the measured contact energy the model predicts that the
stalled material sits in small oligomers, below the size that a 27 nm field registers as a particle,
so TEM would see essentially no malformed particles (counting 200 particles with none malformed
would bound the fraction below 1.5%, and 1000 particles below 0.3%, at 95%). The kinetic trap is
therefore an SEC observable, a shoulder between the dimer and capsid peaks, and not a TEM
observable. The particle number that enters M5 comes from SEC with a protein assay (7.7% CV, against
12.7% for 200 counted particles).

## Experiments supported by the module

- Disassembly midpoint by urea or thermal ramp with SEC, which constrains ΔG<sub>contact</sub>
  directly and is independent of the dPCR arm.
- SEC or DLS on the induction titration, looking for a shoulder between the dimer and capsid peaks.
- In vitro reassembly at 20–50 µM dimer ± TR hairpin.
- TEM reported as the fraction of correct particles (at least 200 particles across at least five
  fields), with its binomial 95% interval, and not as a representative image.
