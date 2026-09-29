import numpy as np
import pandas as pd
import pytest

from reservoir_subsidence.uncertainty import (
    constrained_modulus,
    circular_center_subsidence,
    monte_carlo_center_subsidence,
)


def example_data():
    return pd.DataFrame(
        [{
            "pressure_depletion_bar": 150.0,
            "thickness_m": 25.0,
            "reservoir_depth_m": 2000.0,
            "radius_m": 1500.0,
        }]
    )


def test_constrained_modulus_positive():
    assert constrained_modulus(10e9, 0.25) > 0


def test_zero_depletion_zero_subsidence():
    result = circular_center_subsidence(
        0.0, 30.0, 2000.0, 1500.0, 10e9, 0.25, 0.8
    )
    assert result == 0.0


def test_biot_coefficient_is_bounded():
    with pytest.raises(ValueError):
        circular_center_subsidence(
            15e6, 30.0, 2000.0, 1500.0, 10e9, 0.25, 1.2
        )


def test_poisson_ratio_below_half():
    with pytest.raises(ValueError):
        constrained_modulus(10e9, 0.5)


def test_monte_carlo_reproducible():
    data = example_data()
    a = monte_carlo_center_subsidence(
        data, n_simulations=100, random_seed=7
    )
    b = monte_carlo_center_subsidence(
        data, n_simulations=100, random_seed=7
    )
    assert np.allclose(
        a["total_subsidence_mm"], b["total_subsidence_mm"]
    )


def test_equivalent_compaction_coefficient_positive():
    result = monte_carlo_center_subsidence(
        example_data(), n_simulations=100, random_seed=7
    )
    assert np.all(
        result["equivalent_compaction_coefficient_per_bar"] > 0
    )


def test_invalid_sampling_ranges_rejected():
    data = example_data()

    with pytest.raises(ValueError):
        monte_carlo_center_subsidence(
            data,
            poisson_ratio_range=(0.20, 0.50),
            n_simulations=10,
        )

    with pytest.raises(ValueError):
        monte_carlo_center_subsidence(
            data,
            biot_coefficient_range=(0.70, 1.20),
            n_simulations=10,
        )
