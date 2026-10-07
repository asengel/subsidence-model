# Subsidence Modelling

Python implementation of **Geertsma–Van Opstal subsidence modelling**.

All examples in this repository use **synthetic data**.

![Total subsidence map](figures/total_subsidence_map.png)

## Overview

Pressure depletion increases effective stress in a reservoir and cause reservoir-rock compaction. That compaction is transferred through the surrounding rock mass and result in surface subsidence.

The workflow implemented here is:

```text
Pressure depletion
        ↓
Reservoir compaction
        ↓
Geertsma influence function
        ↓
Van Opstal rigid-basement correction
        ↓
Spatial superposition
        ↓
Surface subsidence
        ↓
Visualization and uncertainty analysis
```

The spatial model operates on individual source cells, allowing pressure depletion, thickness, compaction coefficient, depth, and reservoir geometry to vary spatially.

## Data

All datasets included in this repository are synthetic and are provided only to demonstrate the workflow.

No proprietary reservoir, well, field, client, or employer information is included.

---

## References

- Geertsma, J. (1973). *A Basic Theory of Subsidence Due to Reservoir Compaction: The Homogeneous Case.*
- Van Opstal, G. H. C. (1974). *The Effect of Base-Rock Rigidity on Subsidence Due to Reservoir Compaction.*
- Fokker, P. A. & Orlic, B. (2006). Semi-analytic modelling of subsidence.

---

## Disclaimer

This repository is intended for engineering demonstration, research, and portfolio purposes. It should not be used as the sole basis for operational, regulatory, safety-critical, or infrastructure decisions.
