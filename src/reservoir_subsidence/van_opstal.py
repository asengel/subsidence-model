"""Geertsma and Geertsma–Van Opstal subsidence influence functions."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

VAN_OPSTAL_POISSON_RATIO = 0.25


@dataclass(frozen=True)
class VanOpstalConstants:
    """Van Opstal coefficient set used with ν = 0.25."""

    a11: float = 1.0
    a12: float = -0.778
    a21: float = 2.0
    a22: float = 2.8
    b11: float = -0.2
    b21: float = 4.0


DEFAULT_CONSTANTS = VanOpstalConstants()


def geertsma_influence(
    surface_x_m,
    surface_y_m,
    source_x_m,
    source_y_m,
    reservoir_depth_m,
    poisson_ratio=VAN_OPSTAL_POISSON_RATIO,
):
    """Return the Geertsma vertical influence [1/m²]."""
    if not 0.0 < poisson_ratio < 0.5:
        raise ValueError("poisson_ratio must satisfy 0 < ν < 0.5.")
    if reservoir_depth_m <= 0:
        raise ValueError("reservoir_depth_m must be > 0.")

    x = np.asarray(surface_x_m, dtype=float)
    y = np.asarray(surface_y_m, dtype=float)

    r2 = (x - source_x_m) ** 2 + (y - source_y_m) ** 2
    base = r2 + reservoir_depth_m**2

    return (
        (1.0 - poisson_ratio)
        / (2.0 * np.pi)
        * (2.0 * reservoir_depth_m)
        / (base * np.sqrt(base))
    )


def van_opstal_influence(
    surface_x_m,
    surface_y_m,
    source_x_m,
    source_y_m,
    reservoir_depth_m,
    basement_depth_m,
    poisson_ratio=VAN_OPSTAL_POISSON_RATIO,
    constants=DEFAULT_CONSTANTS,
):
    """Return the Geertsma–Van Opstal vertical influence [1/m²]."""
    if not np.isclose(poisson_ratio, VAN_OPSTAL_POISSON_RATIO):
        raise ValueError(
            "This coefficient set is used with poisson_ratio=0.25."
        )
    if reservoir_depth_m <= 0:
        raise ValueError("reservoir_depth_m must be > 0.")
    if basement_depth_m < reservoir_depth_m:
        raise ValueError(
            "basement_depth_m must be >= reservoir_depth_m."
        )

    x = np.asarray(surface_x_m, dtype=float)
    y = np.asarray(surface_y_m, dtype=float)

    C = float(reservoir_depth_m)
    K = float(basement_depth_m)
    nu = float(poisson_ratio)
    q = 3.0 - 4.0 * nu
    r2 = (x - source_x_m) ** 2 + (y - source_y_m) ** 2

    a11, a12 = constants.a11, constants.a12
    a21, a22 = constants.a21, constants.a22
    b11, b21 = constants.b11, constants.b21

    base = geertsma_influence(
        x, y, source_x_m, source_y_m, C, poisson_ratio=nu
    )

    a1 = a11 * (a21 * K - C) - 2.0 * a11 * K
    a2 = 2.0 * a11 * K - q**2 * a11 * (a21 * K + C)
    a3 = q * a11 * (a21 * K + 2.0 * K + C)
    a4 = q * a11 * (a21 * K + 2.0 * K - C)
    a5 = 6.0 * a11 * K * (a21 * K - C) ** 2
    a6 = (
        6.0 * a11 * K * (a21 * K + C)
        * (6.0 * K - (a21 * K + C))
    )
    a7 = 60.0 * a11 * K**2 * (a21 * K + C) ** 3

    a8 = a12 * (a22 * K - C) - 2.0 * a12 * K
    a9 = 2.0 * a12 * K - q**2 * a12 * (a22 * K + C)
    a10 = q * a12 * (a22 * K + 2.0 * K + C)
    a11_img = q * a12 * (a22 * K + 2.0 * K - C)
    a12_img = 6.0 * a12 * K * (a22 * K - C) ** 2
    a13 = (
        6.0 * a12 * K * (a22 * K + C)
        * (6.0 * K - (a22 * K + C))
    )
    a14 = 60.0 * a12 * K**2 * (a22 * K + C) ** 3

    n1 = r2 + (a21 * K - C) ** 2
    n2 = r2 + (a21 * K + C) ** 2
    n3 = r2 + ((a21 + 2.0) * K + C) ** 2
    n4 = r2 + ((a21 + 2.0) * K - C) ** 2

    sum1a = (
        a1 / (n1 * np.sqrt(n1))
        + a2 / (n2 * np.sqrt(n2))
        - a3 / (n3 * np.sqrt(n3))
        + a4 / (n4 * np.sqrt(n4))
        + a5 / (n1**2 * np.sqrt(n1))
        + a6 / (n2**2 * np.sqrt(n2))
        - a7 / (n2**3 * np.sqrt(n2))
    )

    n1 = r2 + (a22 * K - C) ** 2
    n2 = r2 + (a22 * K + C) ** 2
    n3 = r2 + ((a22 + 2.0) * K + C) ** 2
    n4 = r2 + ((a22 + 2.0) * K - C) ** 2

    sum1b = (
        a8 / (n1 * np.sqrt(n1))
        + a9 / (n2 * np.sqrt(n2))
        - a10 / (n3 * np.sqrt(n3))
        + a11_img / (n4 * np.sqrt(n4))
        + a12_img / (n1**2 * np.sqrt(n1))
        + a13 / (n2**2 * np.sqrt(n2))
        - a14 / (n2**3 * np.sqrt(n2))
    )

    b1 = b11 * K
    b2 = q**2 * b11 * K
    b3 = q * b11 * K
    b4 = 3.0 * b11 * K * (b21 * K - C) * (6.0 * K - (b21 * K - C))
    b5 = (
        3.0 * b11 * K * (b21 * K + C)
        * (6.0 * K - q**2 * (b21 * K + C))
        - 36.0 * b11 * K**3
    )
    b6 = 3.0 * q * b11 * K * (b21 * K + 2.0 * K - C) ** 2
    b7 = 3.0 * q * b11 * K * (b21 * K + 2.0 * K + C) ** 2
    b8 = 30.0 * b11 * K**2 * (b21 * K - C) ** 3
    b9 = (
        30.0 * b11 * K**2 * (b21 * K + C) ** 2
        * (12.0 * K - (b21 * K + C))
    )
    b10 = 420.0 * b11 * K**3 * (b21 * K + C) ** 4

    n1 = r2 + (b21 * K - C) ** 2
    n2 = r2 + (b21 * K + C) ** 2
    n3 = r2 + ((b21 + 2.0) * K - C) ** 2
    n4 = r2 + ((b21 + 2.0) * K + C) ** 2

    sum2 = (
        -b1 / (n1 * np.sqrt(n1))
        + b2 / (n2 * np.sqrt(n2))
        - b3 / (n3 * np.sqrt(n3))
        + b3 / (n4 * np.sqrt(n4))
        - b4 / (n1**2 * np.sqrt(n1))
        + b5 / (n2**2 * np.sqrt(n2))
        + b6 / (n3**2 * np.sqrt(n3))
        - b7 / (n4**2 * np.sqrt(n4))
        + b8 / (n1**3 * np.sqrt(n1))
        + b9 / (n2**3 * np.sqrt(n2))
        - b10 / (n2**4 * np.sqrt(n2))
    )

    correction = (1.0 - nu) / (2.0 * np.pi) * (sum1a + sum1b + sum2)
    return base + correction


def calculate_subsidence(
    surface_x_m,
    surface_y_m,
    source_x_m,
    source_y_m,
    reservoir_depth_m,
    basement_depth_m,
    compaction_volume_m3,
    poisson_ratio=VAN_OPSTAL_POISSON_RATIO,
):
    """Superimpose subsidence from multiple compacting source cells [m]."""
    sx = np.asarray(surface_x_m, dtype=float)
    sy = np.asarray(surface_y_m, dtype=float)

    x0 = np.atleast_1d(np.asarray(source_x_m, dtype=float))
    y0 = np.atleast_1d(np.asarray(source_y_m, dtype=float))
    depth = np.atleast_1d(np.asarray(reservoir_depth_m, dtype=float))
    basement = np.atleast_1d(np.asarray(basement_depth_m, dtype=float))
    compaction = np.atleast_1d(np.asarray(compaction_volume_m3, dtype=float))

    lengths = {len(x0), len(y0), len(depth), len(basement), len(compaction)}
    if len(lengths) != 1:
        raise ValueError("All source arrays must have the same length.")
    if np.any(compaction < 0):
        raise ValueError("compaction_volume_m3 must be non-negative.")

    result = np.zeros(np.broadcast(sx, sy).shape, dtype=float)

    for xi, yi, zi, kb, dv in zip(x0, y0, depth, basement, compaction):
        result += dv * van_opstal_influence(
            sx, sy, xi, yi, zi, kb, poisson_ratio=poisson_ratio
        )

    return result
