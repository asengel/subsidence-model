"""Reservoir compaction calculations."""

from __future__ import annotations

import numpy as np


def calculate_compaction_volume(
    cell_area_m2,
    thickness_m,
    compaction_coefficient_per_bar,
    pressure_depletion_bar,
):
    """Calculate cell-wise reservoir compaction volume [m³].

    Uses the linear relationship

        ΔV = A h C_m Δp

    with C_m in 1/bar and pressure depletion in bar.
    """
    area = np.asarray(cell_area_m2, dtype=float)
    thickness = np.asarray(thickness_m, dtype=float)
    cm = np.asarray(compaction_coefficient_per_bar, dtype=float)
    depletion = np.asarray(pressure_depletion_bar, dtype=float)

    if np.any(area < 0):
        raise ValueError("cell_area_m2 must be non-negative.")
    if np.any(thickness < 0):
        raise ValueError("thickness_m must be non-negative.")
    if np.any(cm < 0):
        raise ValueError("compaction_coefficient_per_bar must be non-negative.")
    if np.any(depletion < 0):
        raise ValueError("pressure_depletion_bar must be non-negative.")

    return area * thickness * cm * depletion
