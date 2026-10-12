# Reservoir subsidence

One compact [Jupyter notebook](subsidence.ipynb) for Geertsma–Van Opstal and Tempone, Fjær & Landrø (2010) surface subsidence.

```bash
pip install -r requirements.txt
jupyter lab subsidence.ipynb
```

Run the notebook from the repository root. It contains all equations and calculation functions, the Figure 3-style comparison, and an optional example using `data/synthetic_reservoir.csv`. The comparison is self-contained. Inputs are assumed valid; there are no tests or input-validation routines.

The calculated Tempone and Van Opstal profiles nearly overlap (about 88 mm at the centre). The published Figure 3 appears closer to 95 mm; the notebook explicitly records this unresolved discrepancy and its source discretization. It does not claim exact reproduction of the published curve.

All reservoir data are synthetic. The notebook focuses on surface displacement; Monte Carlo screening and the former package/test structure are omitted. Earlier versions remain in Git history.

Reference: [Tempone et al. (2010), DOI: 10.1016/j.apm.2010.02.025](https://doi.org/10.1016/j.apm.2010.02.025).
