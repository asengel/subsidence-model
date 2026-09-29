import numpy as np
from reservoir_subsidence.compaction import calculate_compaction_volume


def test_zero_depletion_gives_zero_compaction():
    assert calculate_compaction_volume(40_000, 30, 1e-5, 0) == 0


def test_compaction_is_linear():
    v1 = calculate_compaction_volume(40_000, 30, 1e-5, 100)
    v2 = calculate_compaction_volume(40_000, 30, 1e-5, 200)
    assert np.isclose(v2, 2*v1)


def test_known_compaction_value():
    value = calculate_compaction_volume(200**2, 20, 1e-5, 200)
    assert np.isclose(value, 1600.0)
