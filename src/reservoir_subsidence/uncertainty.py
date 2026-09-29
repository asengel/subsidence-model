"""Simplified analytical uncertainty screening for reservoir subsidence.

This module is intentionally separate from the spatial Van Opstal solver.
It uses a homogeneous Geertsma circular-reservoir centre-point approximation
to demonstrate uncertainty propagation in elastic/poroelastic parameters.

Within each Monte Carlo realization, one sampled E-nu-alpha_B triplet is
applied to all listed reservoirs. This implies fully correlated mechanical
properties across reservoirs in the current demonstration.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def constrained_modulus(youngs_modulus_pa, poisson_ratio):
    """Calculate the drained constrained modulus M [Pa].

    M = E(1-nu) / ((1+nu)(1-2nu))
    """
    E = np.asarray(youngs_modulus_pa, dtype=float)
    nu = np.asarray(poisson_ratio, dtype=float)

    if np.any(~np.isfinite(E)) or np.any(E <= 0):
        raise ValueError("youngs_modulus_pa must contain finite values > 0.")
    if np.any(~np.isfinite(nu)) or np.any((nu <= 0.0) | (nu >= 0.5)):
        raise ValueError("poisson_ratio must satisfy 0 < nu < 0.5.")

    return E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))


def circular_center_subsidence(
    pressure_depletion_pa,
    thickness_m,
    reservoir_depth_m,
    radius_m,
    youngs_modulus_pa,
    poisson_ratio,
    biot_coefficient,
):
    """Calculate simplified Geertsma centre-point subsidence [m].

    Subsidence is returned as positive downward.
    """
    p = np.asarray(pressure_depletion_pa, dtype=float)
    h = np.asarray(thickness_m, dtype=float)
    z = np.asarray(reservoir_depth_m, dtype=float)
    radius = np.asarray(radius_m, dtype=float)
    alpha = np.asarray(biot_coefficient, dtype=float)

    if np.any(~np.isfinite(p)) or np.any(p < 0):
        raise ValueError("pressure_depletion_pa must be finite and non-negative.")
    if np.any(~np.isfinite(h)) or np.any(h < 0):
        raise ValueError("thickness_m must be finite and non-negative.")
    if np.any(~np.isfinite(z)) or np.any(z <= 0):
        raise ValueError("reservoir_depth_m must be finite and > 0.")
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0):
        raise ValueError("radius_m must be finite and > 0.")
    if np.any(~np.isfinite(alpha)) or np.any((alpha <= 0.0) | (alpha > 1.0)):
        raise ValueError("biot_coefficient must satisfy 0 < alpha_B <= 1.")

    M = constrained_modulus(youngs_modulus_pa, poisson_ratio)
    delta_h = alpha * p * h / M
    geometry = 1.0 - z / np.sqrt(z**2 + radius**2)

    return delta_h * 2.0 * (1.0 - poisson_ratio) * geometry


def _validate_range(name, bounds):
    """Validate a finite increasing two-value interval."""
    if len(bounds) != 2:
        raise ValueError(f"{name} must contain exactly two values.")
    lo, hi = map(float, bounds)
    if not np.isfinite(lo) or not np.isfinite(hi):
        raise ValueError(f"{name} values must be finite.")
    if hi <= lo:
        raise ValueError(f"{name} upper bound must be greater than lower bound.")
    return lo, hi


def monte_carlo_center_subsidence(
    reservoir_data: pd.DataFrame,
    *,
    youngs_modulus_gpa_range=(8.0, 15.0),
    poisson_ratio_range=(0.20, 0.30),
    biot_coefficient_range=(0.70, 1.00),
    n_simulations: int = 10_000,
    random_seed: int = 42,
):
    """Run reproducible Monte Carlo screening with uniform independent inputs.

    Required deterministic reservoir columns:
        pressure_depletion_bar
        thickness_m
        reservoir_depth_m
        radius_m
    """
    required = {
        "pressure_depletion_bar",
        "thickness_m",
        "reservoir_depth_m",
        "radius_m",
    }
    missing = required - set(reservoir_data.columns)
    if missing:
        raise ValueError(f"Missing reservoir columns: {sorted(missing)}")

    if not isinstance(n_simulations, (int, np.integer)) or n_simulations <= 0:
        raise ValueError("n_simulations must be a positive integer.")

    e_lo, e_hi = _validate_range(
        "youngs_modulus_gpa_range", youngs_modulus_gpa_range
    )
    nu_lo, nu_hi = _validate_range(
        "poisson_ratio_range", poisson_ratio_range
    )
    a_lo, a_hi = _validate_range(
        "biot_coefficient_range", biot_coefficient_range
    )

    if e_lo <= 0:
        raise ValueError("youngs_modulus_gpa_range must be > 0.")
    if nu_lo <= 0 or nu_hi >= 0.5:
        raise ValueError("poisson_ratio_range must satisfy 0 < nu < 0.5.")
    if a_lo <= 0 or a_hi > 1.0:
        raise ValueError(
            "biot_coefficient_range must satisfy 0 < alpha_B <= 1."
        )

    for column in required:
        values = pd.to_numeric(
            reservoir_data[column], errors="coerce"
        ).to_numpy()
        if np.any(~np.isfinite(values)):
            raise ValueError(f"{column} contains non-finite values.")

    if np.any(reservoir_data["pressure_depletion_bar"].to_numpy() < 0):
        raise ValueError("pressure_depletion_bar must be non-negative.")
    if np.any(reservoir_data["thickness_m"].to_numpy() < 0):
        raise ValueError("thickness_m must be non-negative.")
    if np.any(reservoir_data["reservoir_depth_m"].to_numpy() <= 0):
        raise ValueError("reservoir_depth_m must be > 0.")
    if np.any(reservoir_data["radius_m"].to_numpy() <= 0):
        raise ValueError("radius_m must be > 0.")

    rng = np.random.default_rng(random_seed)

    E_gpa = rng.uniform(e_lo, e_hi, n_simulations)
    nu = rng.uniform(nu_lo, nu_hi, n_simulations)
    alpha = rng.uniform(a_lo, a_hi, n_simulations)

    total_m = np.zeros(n_simulations, dtype=float)

    for row in reservoir_data.itertuples(index=False):
        total_m += circular_center_subsidence(
            pressure_depletion_pa=row.pressure_depletion_bar * 1e5,
            thickness_m=row.thickness_m,
            reservoir_depth_m=row.reservoir_depth_m,
            radius_m=row.radius_m,
            youngs_modulus_pa=E_gpa * 1e9,
            poisson_ratio=nu,
            biot_coefficient=alpha,
        )

    equivalent_cm_per_bar = (
        alpha
        / constrained_modulus(E_gpa * 1e9, nu)
        * 1e5
    )

    return pd.DataFrame(
        {
            "youngs_modulus_gpa": E_gpa,
            "poisson_ratio": nu,
            "biot_coefficient": alpha,
            "equivalent_compaction_coefficient_per_bar": equivalent_cm_per_bar,
            "total_subsidence_mm": total_m * 1000.0,
        }
    )
