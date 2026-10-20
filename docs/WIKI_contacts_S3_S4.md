<!--
DRAFT, NOT YET FOR THE WIKI. The numbers come from contacts.py runs on the 30 ns trajectories of S3
(5 replicas) and S4 (3 replicas), analysed over ALL frames (0-30 ns), not with the protocol's
5 ns discard, and the outputs are not in this repository. The wiki rule is that every number has its
file here. Before this goes into wiki/02_md_methods.md (or a new page):
  1. Rerun:   python3 analysis/contacts.py S3_cp_operator --equil-ns 5 --rna-first -12 \
                  --key-residues 29,43,45,47,49,52,55,57,59,61,63,85,87
              (same for S4, and S9 if wanted) and commit analysis_out/contacts/<system>/.
     The RNA is residues 4-16 of the 19-mer, so the first nucleotide is -12, NOT -15.
  2. Replace the numbers below with the new ones, and use all 5 replicas of S4 if they exist.
  3. Look at replica 3 of S3 (it deviated) and say what is seen.
-->

## Protein–RNA contacts in the MS2 coat protein–operator complex

**What was done.** S3 is the coat protein dimer bound to the wild-type operator fragment (U at −5) and S4 to
the crystal RNA (the C(−5) variant, the *pac* hairpin of the campaign). The two systems differ in a single
nucleotide. For each frame we counted the residue pairs (protein residue × nucleotide) with any heavy-atom
distance below 4.5 Å and took the occupancy of each pair, the fraction of frames in which it is in contact.
Occupancies are means over replicas and error bars are the standard deviation **across replicas**.

**The complex stays together.** In the five S3 replicas the RNA remained bound for the whole 30 ns (bound
fraction 1.00; closest protein–RNA approach 2.57 Å in every replica).

| S3, mean ± SD across 5 replicas | |
|---|---|
| residue-pair contacts per frame | 39.4 ± 1.8 |
| protein residues in contact | 26.9 ± 1.3 |
| nucleotides in contact | 8.7 ± 0.5 |
| protein–RNA hydrogen bonds per frame | 11.4 ± 0.8 |
| residue pairs at the dimer interface | 176 ± 5.5 |

**The known contacts are recovered.** The structures of the coat protein bound to the operator (van den Worm
*et al.*, 1998) show that the hairpin binds the A/B dimer asymmetrically: Val29, Thr45, Ser47, Thr59 and Lys61
form the binding pocket, and Lys43, Arg49, Ser52, Asn55, Lys57, Glu63, Tyr85 and Asn87 contact the RNA. In S3
all 13 of these residues are among the 30 residues that contact the RNA most often (4.5 Å); with a 3.5 Å
cutoff (direct contacts) 11 of the 13 are among the top 15, Val29 (a hydrophobic contact) and Asn55 dropping
out. The asymmetry is reproduced: subunit A holds the pocket (Thr45, Ser47, Thr59, Lys61, Glu63, Tyr85 and
Asn87, occupancy ≥ 0.94), and subunit B contacts the backbone (Arg49, Ser52, Asn55, Lys57 and Lys61,
occupancy ≥ 0.84).

**S4 against S3, preliminary.** The pocket residues are kept in S4 (Thr45, Ser47, Thr59, Lys61, Glu63 and Tyr85
of subunit A, occupancy ≥ 0.95). Two contacts decrease: Asn87 (0.99 → 0.61 ± 0.35 in A and 0.85 → 0.50 ± 0.36
in B) and Tyr85 of subunit B (0.50 → 0.17 ± 0.15). In the crystal structures Asn87 and Tyr85 contact the base at
−5, which is the nucleotide that differs between S3 and S4, so a change there is expected; the data here do not
show that it is the cause. The spread between replicas is large and S4 has 3 replicas against 5.

**Limits.**
- The runs are 30 ns long and sample only the bound basin.
- Contacts are geometric. They are not binding energies and say nothing about affinity.
- Water-mediated bridges and Mg²⁺ contacts are not counted (solute-only trajectories).
- Replica 3 of S3 deviates from the others in several subunit-B contacts.
- The scrambled control was dropped, so no specificity comparison exists.

*Reference:* van den Worm SHE, Stonehouse NJ, Valegård K, Murray JB, Walton C, Fridborg K, Stockley PG,
Liljas L. Crystal structures of MS2 coat protein mutants in complex with wild-type RNA operator fragments.
*Nucleic Acids Res.* 1998;26(5):1345–1351.
