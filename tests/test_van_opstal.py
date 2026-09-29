import numpy as np
from reservoir_subsidence.van_opstal import (
    calculate_subsidence,
    geertsma_influence,
    van_opstal_influence,
)


def test_geertsma_radial_symmetry():
    r = 1000.0
    values = geertsma_influence(
        np.array([r, -r, 0, 0]),
        np.array([0, 0, r, -r]),
        0,
        0,
        2000,
    )
    assert np.allclose(values, values[0])


def test_van_opstal_radial_symmetry():
    r = 1000.0
    values = van_opstal_influence(
        np.array([r, -r, 0, 0]),
        np.array([0, 0, r, -r]),
        0,
        0,
        2000,
        2200,
    )
    assert np.allclose(values, values[0])


def test_far_field_decreases():
    near = geertsma_influence(0, 0, 0, 0, 2000)
    far = geertsma_influence(100_000, 0, 0, 0, 2000)
    assert 0 < far < near


def test_superposition():
    x = np.array([0.0, 500.0, 1000.0])
    y = np.zeros_like(x)

    combined = calculate_subsidence(
        x, y,
        [-500.0, 500.0],
        [0.0, 0.0],
        [2000.0, 2000.0],
        [2200.0, 2200.0],
        [1000.0, 2000.0],
    )

    first = calculate_subsidence(
        x, y, [-500.0], [0.0], [2000.0], [2200.0], [1000.0]
    )
    second = calculate_subsidence(
        x, y, [500.0], [0.0], [2000.0], [2200.0], [2000.0]
    )

    assert np.allclose(combined, first + second)
