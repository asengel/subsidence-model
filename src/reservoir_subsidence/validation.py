"""Validation and discretization studies for the subsidence model."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .compaction import calculate_compaction_volume
from .grid import make_surface_grid
from .van_opstal import (
    calculate_subsidence,
    geertsma_influence,
    van_opstal_influence,
)


def integrated_axisymmetric_influence(
    model: str,
    *,
    reservoir_depth_m: float,
    basement_depth_m: float | None = None,
    integration_radius_m: float = 500_000.0,
    n_radial_points: int = 5_000,
    poisson_ratio: float = 0.25,
) -> float:
    """Numerically integrate an axisymmetric point-source influence function.

    The dimensionless result is

        integral_0^R 2*pi*r*G(r) dr

    which is also the surface-subsidence-volume / source-compaction-volume
    ratio for a unit compaction source.

    Parameters
    ----------
    model : {"geertsma", "van_opstal"}
        Influence function to integrate.
    reservoir_depth_m : float
        Positive source depth below surface [m].
    basement_depth_m : float, optional
        Positive rigid-basement depth [m]. Required for ``van_opstal``.
    integration_radius_m : float
        Maximum radial extent of the surface integration [m].
    n_radial_points : int
        Number of logarithmically distributed non-zero radial samples.
    poisson_ratio : float
        Poisson's ratio.
    """
    if reservoir_depth_m <= 0:
        raise ValueError("reservoir_depth_m must be > 0.")
    if integration_radius_m <= 0:
        raise ValueError("integration_radius_m must be > 0.")
    if n_radial_points < 100:
        raise ValueError("n_radial_points must be >= 100.")

    # Include r = 0 explicitly and logarithmically sample the long tail.
    r = np.concatenate(
        [
            np.array([0.0]),
            np.geomspace(
                max(reservoir_depth_m * 1e-8, 1e-6),
                integration_radius_m,
                n_radial_points,
            ),
        ]
    )
    y = np.zeros_like(r)

    model_key = model.lower()

    if model_key == "geertsma":
        influence = geertsma_influence(
            r,
            y,
            source_x_m=0.0,
            source_y_m=0.0,
            reservoir_depth_m=reservoir_depth_m,
            poisson_ratio=poisson_ratio,
        )
    elif model_key == "van_opstal":
        if basement_depth_m is None:
            raise ValueError(
                "basement_depth_m is required for the Van Opstal model."
            )
        influence = van_opstal_influence(
            r,
            y,
            source_x_m=0.0,
            source_y_m=0.0,
            reservoir_depth_m=reservoir_depth_m,
            basement_depth_m=basement_depth_m,
            poisson_ratio=poisson_ratio,
        )
    else:
        raise ValueError("model must be 'geertsma' or 'van_opstal'.")

    integrand = 2.0 * np.pi * r * influence
    return float(np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)))


def bowl_volume_convergence(
    *,
    reservoir_depth_m: float = 2_000.0,
    basement_depth_m: float = 2_200.0,
    radius_multipliers=(5, 10, 20, 50, 100, 250),
    n_radial_points: int = 5_000,
) -> pd.DataFrame:
    """Evaluate numerical bowl-volume ratios over increasing integration radius."""
    rows = []

    for multiplier in radius_multipliers:
        radius = float(multiplier) * reservoir_depth_m

        geertsma_ratio = integrated_axisymmetric_influence(
            "geertsma",
            reservoir_depth_m=reservoir_depth_m,
            integration_radius_m=radius,
            n_radial_points=n_radial_points,
        )
        van_opstal_ratio = integrated_axisymmetric_influence(
            "van_opstal",
            reservoir_depth_m=reservoir_depth_m,
            basement_depth_m=basement_depth_m,
            integration_radius_m=radius,
            n_radial_points=n_radial_points,
        )

        rows.append(
            {
                "integration_radius_over_depth": float(multiplier),
                "integration_radius_m": radius,
                "geertsma_volume_ratio": geertsma_ratio,
                "van_opstal_volume_ratio": van_opstal_ratio,
            }
        )

    return pd.DataFrame(rows)


def synthetic_source_cells(source_spacing_m: float) -> pd.DataFrame:
    """Create the two-reservoir synthetic model at a chosen source-cell spacing."""
    if source_spacing_m <= 0:
        raise ValueError("source_spacing_m must be > 0.")

    records = []

    reservoirs = [
        # id, centre x/y, ellipse axes, depth, thickness, Cm, depletion, basement factor
        ("A", 0.0, 0.0, 1800.0, 1200.0, 2100.0, 28.0, 1.1e-5, 180.0, 1.10),
        ("B", 3200.0, 1000.0, 1200.0, 900.0, 2500.0, 22.0, 0.9e-5, 145.0, 1.15),
    ]

    for (
        reservoir_id,
        cx,
        cy,
        semi_major,
        semi_minor,
        depth,
        thickness,
        cm,
        centre_depletion,
        basement_factor,
    ) in reservoirs:
        x_values = np.arange(
            cx - semi_major,
            cx + semi_major + 0.5 * source_spacing_m,
            source_spacing_m,
        )
        y_values = np.arange(
            cy - semi_minor,
            cy + semi_minor + 0.5 * source_spacing_m,
            source_spacing_m,
        )

        for x in x_values:
            for y in y_values:
                q = (
                    ((x - cx) / semi_major) ** 2
                    + ((y - cy) / semi_minor) ** 2
                )

                if q <= 1.0:
                    records.append(
                        {
                            "reservoir_id": reservoir_id,
                            "x_m": x,
                            "y_m": y,
                            "reservoir_depth_m": depth,
                            "thickness_m": thickness
                            * (0.90 + 0.10 * (1.0 - q)),
                            "compaction_coefficient_per_bar": cm,
                            "pressure_depletion_bar": centre_depletion
                            * (0.75 + 0.25 * (1.0 - q)),
                            "basement_depth_factor": basement_factor,
                            "source_cell_area_m2": source_spacing_m**2,
                        }
                    )

    data = pd.DataFrame(records)

    data["compaction_volume_m3"] = calculate_compaction_volume(
        data["source_cell_area_m2"],
        data["thickness_m"],
        data["compaction_coefficient_per_bar"],
        data["pressure_depletion_bar"],
    )
    data["basement_depth_m"] = (
        data["reservoir_depth_m"] * data["basement_depth_factor"]
    )

    return data


def source_grid_convergence(
    *,
    source_spacings_m=(400.0, 200.0, 100.0),
    surface_spacing_m: float = 200.0,
    observation_x_m: float = 0.0,
    observation_y_m: float = 0.0,
) -> pd.DataFrame:
    """Run the synthetic model at progressively finer source-cell spacing.

    The surface evaluation grid is held fixed so that changes reflect source
    discretization rather than a simultaneous change in output-grid resolution.
    """
    _, _, xx, yy = make_surface_grid(
        xmin_m=-4500.0,
        xmax_m=6500.0,
        ymin_m=-4000.0,
        ymax_m=4500.0,
        spacing_m=surface_spacing_m,
    )

    rows = []

    for spacing in source_spacings_m:
        data = synthetic_source_cells(float(spacing))
        total_surface_m = np.zeros_like(xx, dtype=float)
        observation_subsidence_m = 0.0

        for _, group in data.groupby("reservoir_id"):
            total_surface_m += calculate_subsidence(
                xx,
                yy,
                group["x_m"].to_numpy(),
                group["y_m"].to_numpy(),
                group["reservoir_depth_m"].to_numpy(),
                group["basement_depth_m"].to_numpy(),
                group["compaction_volume_m3"].to_numpy(),
            )

            observation_subsidence_m += float(
                calculate_subsidence(
                    np.array([observation_x_m]),
                    np.array([observation_y_m]),
                    group["x_m"].to_numpy(),
                    group["y_m"].to_numpy(),
                    group["reservoir_depth_m"].to_numpy(),
                    group["basement_depth_m"].to_numpy(),
                    group["compaction_volume_m3"].to_numpy(),
                )[0]
            )

        surface_cell_area_m2 = surface_spacing_m**2
        integrated_bowl_volume_m3 = float(
            total_surface_m.sum() * surface_cell_area_m2
        )
        total_compaction_m3 = float(data["compaction_volume_m3"].sum())

        rows.append(
            {
                "source_spacing_m": float(spacing),
                "source_cells": int(len(data)),
                "total_compaction_m3": total_compaction_m3,
                "maximum_subsidence_mm": float(total_surface_m.max() * 1000.0),
                "observation_subsidence_mm": float(
                    observation_subsidence_m * 1000.0
                ),
                "integrated_bowl_volume_m3": integrated_bowl_volume_m3,
                "bowl_to_compaction_ratio": (
                    integrated_bowl_volume_m3 / total_compaction_m3
                ),
            }
        )

    result = pd.DataFrame(rows).sort_values(
        "source_spacing_m", ascending=False
    ).reset_index(drop=True)

    finest = result.iloc[-1]

    for column in [
        "maximum_subsidence_mm",
        "observation_subsidence_mm",
        "integrated_bowl_volume_m3",
    ]:
        result[f"{column}_error_vs_finest_pct"] = (
            (result[column] - finest[column]).abs()
            / abs(finest[column])
            * 100.0
        )

    return result
