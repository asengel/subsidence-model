"""Surface-grid helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd


def make_surface_grid(xmin_m, xmax_m, ymin_m, ymax_m, spacing_m):
    """Create a regular surface evaluation grid."""
    if xmax_m <= xmin_m or ymax_m <= ymin_m:
        raise ValueError("Grid maxima must be greater than minima.")
    if spacing_m <= 0:
        raise ValueError("spacing_m must be > 0.")

    x = np.arange(xmin_m, xmax_m + spacing_m, spacing_m, dtype=float)
    y = np.arange(ymin_m, ymax_m + spacing_m, spacing_m, dtype=float)
    xx, yy = np.meshgrid(x, y)
    return x, y, xx, yy


def results_dataframe(xx, yy, reservoir_results):
    """Create a tabular surface-grid result."""
    data = {
        "X_m": np.asarray(xx).ravel(),
        "Y_m": np.asarray(yy).ravel(),
    }

    for name, values in reservoir_results.items():
        data[f"Subsidence_{name}_m"] = np.asarray(values).ravel()

    df = pd.DataFrame(data)
    subsidence_columns = [c for c in df if c.startswith("Subsidence_")]
    df["Total_Subsidence_m"] = df[subsidence_columns].sum(axis=1)
    return df
