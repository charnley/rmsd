import pathlib

import numpy as np
from context import RESOURCE_PATH

import rmsd


def test_kabash_algorith_pdb():

    filename_1 = pathlib.PurePath(RESOURCE_PATH, "ci2_1.pdb")
    filename_2 = pathlib.PurePath(RESOURCE_PATH, "ci2_2.pdb")

    p_atoms, p_coord = rmsd.get_coordinates_pdb(filename_1)
    q_atoms, q_coord = rmsd.get_coordinates_pdb(filename_2)

    rotation_matrix = rmsd.kabsch(p_coord, q_coord)

    np.testing.assert_array_almost_equal([-0.5124, 0.8565, 0.0608], rotation_matrix[0], decimal=3)


def test_kabash_rotate_pdb():

    filename_1 = pathlib.PurePath(RESOURCE_PATH, "ci2_1.pdb")
    filename_2 = pathlib.PurePath(RESOURCE_PATH, "ci2_2.pdb")

    p_atoms, p_coord = rmsd.get_coordinates_pdb(filename_1)
    q_atoms, q_coord = rmsd.get_coordinates_pdb(filename_2)

    new_p_coord = rmsd.kabsch_rotate(p_coord, q_coord)

    np.testing.assert_array_almost_equal([10.6822, -2.8867, 12.6977], new_p_coord[0], decimal=3)
