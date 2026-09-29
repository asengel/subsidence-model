import numpy as np

from reservoir_subsidence.validation import (
    bowl_volume_convergence,
    integrated_axisymmetric_influence,
    source_grid_convergence,
)


def test_geertsma_integrated_bowl_volume_ratio():
    ratio = integrated_axisymmetric_influence(
        "geertsma",
        reservoir_depth_m=2000.0,
        integration_radius_m=500_000.0,
        n_radial_points=2_000,
        poisson_ratio=0.25,
    )

    # Infinite-domain analytical target is 2(1-nu) = 1.5.
    assert np.isclose(ratio, 1.5, atol=0.01)


def test_van_opstal_integrated_bowl_volume_conservation():
    ratio = integrated_axisymmetric_influence(
        "van_opstal",
        reservoir_depth_m=2000.0,
        basement_depth_m=2200.0,
        integration_radius_m=500_000.0,
        n_radial_points=2_000,
        poisson_ratio=0.25,
    )

    # Rigid-basement formulation approaches volume conservation.
    assert np.isclose(ratio, 1.0, atol=0.005)


def test_bowl_volume_converges_with_integration_radius():
    result = bowl_volume_convergence(
        radius_multipliers=(10, 50, 250),
        n_radial_points=1_500,
    )

    geertsma_error = np.abs(
        result["geertsma_volume_ratio"].to_numpy() - 1.5
    )
    van_opstal_error = np.abs(
        result["van_opstal_volume_ratio"].to_numpy() - 1.0
    )

    assert geertsma_error[-1] < geertsma_error[0]
    assert van_opstal_error[-1] < van_opstal_error[0]


def test_source_grid_refinement_reduces_peak_error():
    result = source_grid_convergence(
        source_spacings_m=(400.0, 200.0, 100.0),
        surface_spacing_m=200.0,
    )

    errors = result[
        "maximum_subsidence_mm_error_vs_finest_pct"
    ].to_numpy()

    assert errors[1] < errors[0]
    assert errors[-1] == 0.0
