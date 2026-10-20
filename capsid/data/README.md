# Data

`litdata/` holds the literature data pack, `templates/` the empty CSV files for bench data, and
`M4_bench_data_template.csv` the template of the exploratory fitting notebook. `FILES.json` lists what
ships with the repository, what has to be supplied, and what the notebooks write. The tables below are
in `litdata/`.

## Literature data pack

Published measurements used by the literature variants (`notebooks/literature/`) and by notebook 04.
Every table that was transcribed from a paper carries its source in the last columns (`source_*`),
and all sources are listed in `lit_sources.csv`.

| File | Rows | Content | Used by |
|---|---|---|---|
| `lit_sources.csv` | 4 | one row per source: module, quantity, authors, year, journal, DOI, PMCID, table or section, and the conversion applied | all |
| `lit_M1_induction_design.csv` | 45 | yield (mg per mL of culture) of three recombinant proteins in a 15-run Box–Behnken design over temperature, time and IPTG; Gutiérrez-González et al. 2019, Table 2 | 01L |
| `lit_M1_fedbatch_timepoints.csv` | 3 | fed-batch volumetric production at 22/25/30 °C and 6/8/10 h, with OD₆₀₀ and DCW; Carnicer et al. 2026, Tables 4–5. Reference only: the fed-batch regime lies outside the shake-flask assumptions of the model | none |
| `lit_M2_partition_gels.csv` | 2 | soluble and insoluble target protein (arenin) at 30 and 37 °C with replicate SD; Hidalgo-Vázquez et al. 2025, Table 1 | 03L |
| `lit_M3_ms2_thermal_stability.csv` | 5 | melting midpoints of wild-type and mutant MS2 VLPs by thermal shift assay; Pistono et al. 2024 | 04L (through the derived table) |
| `derived_M2_arrhenius_from_arenin.csv` | 9 | quantities derived from the previous table: soluble fractions, flux ratio, `s50` ratio, ΔE/R, and the forward predictions at 25 and 18 °C | 03L |
| `derived_M3_dG_from_Tm.csv` | 15 | contact energy obtained by inverting each melting point at the stated assay concentrations (2.5, 5 and 10 µM dimer) | 04, 04L |
| `derived_M3_assembly_at_literature_dG.csv` | 4 | M3 outputs at the placeholder and at the three literature contact energies | 04L |

## Structure-derived data

`structure/` holds the contact schedule of the *T* = 3 capsid, enumerated from PDB 2MS2 by
`../tools/enumerate_contacts.py` (heavy-atom cutoff 4.5 Å, an interface when at least 5 residue pairs are
in contact, 4000 growth paths per rule). `contact_schedule.csv` has, for each step k = 2…90, the contacts
formed and the number of boundary sites for random growth and for the two greedy rules;
`contact_schedule.json` has the same series with the edge classes and the mean interface weights. The
figure is `../figures/fig29_M3_contact_schedule.png`, made by `../tools/plot_contact_schedule.py`.
Notebook 04 does not read these files yet.

## Notes on the columns

- `lit_M2_partition_gels.csv`: the columns `band_soluble`, `band_pellet` and their SDs hold
  concentrations in mg of target protein per litre of culture (BCA combined with densitometric
  purity; see the `units` column), not raw band volumes. `lactose_uM` is empty because the
  inducer was 0.5 mM IPTG (`inducer` column).
- `lit_M1_induction_design.csv`: the centre-point runs are the only replicated conditions; the
  notebook identifies them by their factor levels.
- `derived_M3_dG_from_Tm.csv`: the melting midpoint is treated as a disassembly midpoint, so that at
  `Tm` the assay concentration equals the critical concentration and Eq. (4.8) of the framework
  document can be inverted for the contact energy. The 5 µM row for wild-type MS2
  (−3.82 kcal mol⁻¹) is the value used by notebook 04.
- Files with a `derived_` prefix are computed from the published values with the equations of the
  framework document. `derived_M3_dG_from_Tm.csv` is read directly by notebooks 04 and 04L; the
  other two are included so that the derived numbers can be compared without rerunning the
  notebooks.
