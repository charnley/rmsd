import numpy as np
from conftest import RESOURCE_PATH  # type: ignore

import rmsd as rmsdlib


def test_kabash_algorith_rmsd() -> None:

    filename_1 = RESOURCE_PATH / "ci2_1.pdb"
    filename_2 = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsdlib.get_coordinates(filename_1, "pdb")
    _, q_coord = rmsdlib.get_coordinates(filename_2, "pdb")

    value = rmsdlib.kabsch_rmsd(p_coord, q_coord, translate=True)

    np.testing.assert_almost_equal(value, 11.7768, decimal=4)


def test_kabash_algorith_pdb() -> None:

    filename_1 = RESOURCE_PATH / "ci2_1.pdb"
    filename_2 = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsdlib.get_coordinates_pdb(filename_1)
    _, q_coord = rmsdlib.get_coordinates_pdb(filename_2)

    rotation_matrix = rmsdlib.kabsch(p_coord, q_coord)

    np.testing.assert_array_almost_equal([-0.5124, 0.8565, 0.0608], rotation_matrix[0], decimal=3)


def test_kabash_rotate_pdb() -> None:

    filename_1 = RESOURCE_PATH / "ci2_1.pdb"
    filename_2 = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsdlib.get_coordinates_pdb(filename_1)
    _, q_coord = rmsdlib.get_coordinates_pdb(filename_2)

    new_p_coord = rmsdlib.kabsch_rotate(p_coord, q_coord)

    np.testing.assert_array_almost_equal([10.6822, -2.8867, 12.6977], new_p_coord[0], decimal=3)
