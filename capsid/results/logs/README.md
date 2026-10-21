# Run logs

Cell-by-cell output (timing and printed results) of the runs that produced the numbers reported in
the write-up. The logs were produced with `tools/exec_nb.py`.

| Log | Notebook |
|---|---|
| `log_nb01v.txt` | `BruSafe_01_M1_expression_IPTG.ipynb` |
| `log_nb01L.txt` | `literature/BruSafe_01L_M1_literature_design.ipynb` |
| `log_nb02v.txt` | `BruSafe_02_M1_expression_lactose.ipynb` |
| `log_nb03u.txt` | `BruSafe_03_M2_soluble_partition.ipynb` |
| `log_nb03L.txt` | `literature/BruSafe_03L_M2_literature_arenin.ipynb` |
| `log_nb04u.txt` | `BruSafe_04_M3_capsid_assembly.ipynb` (default, placeholder schedule) |
| `log_nb04s.txt` | `BruSafe_04_M3_capsid_assembly.ipynb` with `SCHEDULE_SOURCE = "structure"` |
| `log_nb04L.txt` | `literature/BruSafe_04L_M3_literature_dG.ipynb` |
| `log_nb05u.txt` | `BruSafe_05_M5_absolute_quantification.ipynb` |

The logs are records of the original runs and have not been edited. The notebooks' explanatory text
was revised afterwards (the code and the numerical results were not changed), so printed messages in
the logs may differ in wording from the current notebooks. The log of notebook 04 predates the TEM
section (§12) of the notebook.

One correction postdates the run of notebook 01: the `Xmax` prior inside `build_model` was written in
OD600 units (6.0 ± 1.5, bounds 2–12) although the model state is biomass in gDCW L⁻¹. It now matches the
prior table (2.40 ± 0.60, bounds 0.8–4.8). The logged posterior for `Xmax` (2.399 ± 0.040) is set by the
data and is not expected to change appreciably, but the notebook has not been re-run with the corrected prior.

The two logs of notebook 04 were regenerated after §9b was added. They were run on a machine without the SF Pro
font, so the repeated `findfont` warnings of matplotlib were removed from them; the timings are those of that machine.
`log_nb04s.txt` is the same notebook with `SCHEDULE_SOURCE = "structure"`; §9b prints the same comparison in both.
