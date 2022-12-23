import numpy as np

import rmsd as rmsdlib


def test_centroid() -> None:
    a1 = np.array([-19.658, 17.18, 25.163], dtype=float)
    a2 = np.array([-20.573, 18.059, 25.88], dtype=float)
    a3 = np.array([-22.018, 17.551, 26.0], dtype=float)

    atms = np.asarray([a1, a2, a3])
    centroid = rmsdlib.centroid(atms)

    assert 3 == len(centroid)

    np.testing.assert_array_almost_equal([-20.7496, 17.5966, 25.6810], centroid, decimal=3)
