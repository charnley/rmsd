import numpy as np
from conftest import RESOURCE_PATH  # type: ignore

import rmsd as rmsdlib


def test_kabash_fit_pdb() -> None:

    filename_p = RESOURCE_PATH / "ci2_1r+t.pdb"
    filename_q = RESOURCE_PATH / "ci2_1.pdb"

    _, p_coord = rmsdlib.get_coordinates_pdb(filename_p)
    _, q_coord = rmsdlib.get_coordinates_pdb(filename_q)

    new_p_coord = rmsdlib.kabsch_fit(p_coord, q_coord)

    np.testing.assert_array_almost_equal(q_coord[0], new_p_coord[0], decimal=2)


def test_kabash_weighted_fit_pdb() -> None:

    filename_1 = RESOURCE_PATH / "ci2_12.pdb"
    filename_2 = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsdlib.get_coordinates_pdb(filename_1)
    _, q_coord = rmsdlib.get_coordinates_pdb(filename_2)

    weights = np.zeros(len(p_coord))
    residue13_start = 200
    residue24_start = 383

    weights[residue13_start:residue24_start] = 1.0

    new_p_coord = rmsdlib.kabsch_fit(p_coord, q_coord, weights)

    np.testing.assert_array_almost_equal(q_coord[300], new_p_coord[300], decimal=2)
