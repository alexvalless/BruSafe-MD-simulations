# Wiki source pages

Source text for the computational pages of the iGEM Tec-Chihuahua 2026 wiki. The pages are written
in Markdown, are self-contained, and refer to the figures in `../capsid/figures/` and to the
notebooks in `../capsid/notebooks/`. Every number quoted in them comes from a notebook, a log file or
a data table in this repository, and the source is named in the text.

| Page | Content |
|---|---|
| [01 · Overview](01_overview.md) | the two computational arms and how they connect |
| [02 · Molecular dynamics](02_md_methods.md) | force field, systems, protocol, convergence checks, reporting rules |
| [03 · M1 Expression](03_capsid_M1_expression.md) | IPTG and lactose induction models, identifiability, validation |
| [04 · M2 Soluble partition](04_capsid_M2_partition.md) | flux-dependent partition between folding and aggregation |
| [05 · M3 Capsid assembly](05_capsid_M3_assembly.md) | nucleation–elongation cascade, contact energy, kinetic trap |
| [06 · M5 Quantification](06_capsid_M5_quantification.md) | cargo copies per particle, error budget, scale-up ranking |
| [07 · Data and sources](07_data_and_sources.md) | literature data pack and its provenance |
| [08 · Reproducibility](08_reproducibility.md) | environment, how to run the notebooks, how to build the MD report pages |
| [09 · Limitations](09_limitations.md) | what has and has not been validated |

## Publishing notes

- The iGEM wiki does not allow external CDNs. The interactive MD pages produced by
  `analysis/md_report.py wiki` (see `docs/TUTORIAL.md`, §10b) are self-contained when built with the
  default `--js local` option, and their whole `html/` folder is uploaded.
- Figures are PNG files at 150 dpi. They should be uploaded with the wiki's file manager and the
  image links in these pages changed to the uploaded URLs.
- The figure captions in the pages are written to be used as they are.
- Page 09 should be published together with the results pages: the limitations are part of the
  results.
