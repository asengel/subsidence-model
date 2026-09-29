"""Plotting utilities."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_subsidence_map(
    xx,
    yy,
    subsidence_m,
    *,
    title="Predicted Total Reservoir-Induced Subsidence",
    output_path=None,
):
    """Plot a filled-contour subsidence map in millimetres."""
    z_mm = np.asarray(subsidence_m) * 1000.0

    fig, ax = plt.subplots(figsize=(9, 7))
    contour = ax.contourf(xx, yy, z_mm, levels=40, cmap="viridis")
    cbar = fig.colorbar(contour, ax=ax)
    cbar.set_label("Subsidence [mm]")

    ax.set_title(title)
    ax.set_xlabel("Easting [m]")
    ax.set_ylabel("Northing [m]")
    ax.set_aspect("equal")
    fig.tight_layout()

    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, bbox_inches="tight")

    return fig, ax


def plot_center_profile(
    x_coordinates_m,
    subsidence_grid_m,
    *,
    title="Subsidence Profile Through Model Centre",
    output_path=None,
):
    """Plot the horizontal profile through the centre of the model."""
    z = np.asarray(subsidence_grid_m)
    mid_row = z.shape[0] // 2

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x_coordinates_m, z[mid_row, :] * 1000.0)
    ax.set_xlabel("Easting [m]")
    ax.set_ylabel("Subsidence [mm]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=200, bbox_inches="tight")

    return fig, ax
