# Reservoir subsidence

One compact [Jupyter notebook](subsidence.ipynb) for Geertsma–Van Opstal and Tempone, Fjær & Landrø (2010) surface subsidence.

```bash
pip install -r requirements.txt
jupyter lab subsidence.ipynb
```

Run the notebook from the repository root. It contains all equations and calculation functions, the Figure 3-style comparison, and an optional example using `data/synthetic_reservoir.csv`. The comparison is self-contained. Inputs are assumed valid; there are no tests or input-validation routines.

Interpreting the stated 800 m depth as the reservoir bottom gives a 750 m source centre and a 950 m basement (150 m beneath the bottom). The calculated centre subsidence is 95.101 mm for Tempone and 94.830 mm for Van Opstal, consistent with Figure 3. Their maximum profile difference is 0.271 mm. The notebook documents its 715-nucleus discretization; it does not claim recovery of the authors’ original cell coordinates or a point-by-point match to published numerical data.

All reservoir data are synthetic. The notebook focuses on surface displacement; Monte Carlo screening and the former package/test structure are omitted. Earlier versions remain in Git history.

Reference: [Tempone et al. (2010), DOI: 10.1016/j.apm.2010.02.025](https://doi.org/10.1016/j.apm.2010.02.025).
