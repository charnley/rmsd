import pathlib

import numpy as np
from constants import RESOURCE_PATH

import rmsd


def test_kabash_fit_pdb():

    filename_p = pathlib.PurePath(RESOURCE_PATH, "ci2_1r+t.pdb")
    filename_q = pathlib.PurePath(RESOURCE_PATH, "ci2_1.pdb")

    p_atoms, p_coord = rmsd.get_coordinates_pdb(filename_p)
    q_atoms, q_coord = rmsd.get_coordinates_pdb(filename_q)

    new_p_coord = rmsd.kabsch_fit(p_coord, q_coord)

    np.testing.assert_array_almost_equal(q_coord[0], new_p_coord[0], decimal=2)


def test_kabash_weighted_fit_pdb():

    filename_1 = pathlib.PurePath(RESOURCE_PATH, "ci2_12.pdb")
    filename_2 = pathlib.PurePath(RESOURCE_PATH, "ci2_2.pdb")

    p_atoms, p_coord = rmsd.get_coordinates_pdb(filename_1)
    q_atoms, q_coord = rmsd.get_coordinates_pdb(filename_2)

    weights = np.zeros(len(p_coord))
    residue13_start = 200
    residue24_start = 383

    weights[residue13_start:residue24_start] = 1.0

    new_p_coord = rmsd.kabsch_fit(p_coord, q_coord, weights)

    np.testing.assert_array_almost_equal(
        q_coord[300], new_p_coord[300], decimal=2
    )
