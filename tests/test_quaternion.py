import numpy as np
from conftest import RESOURCE_PATH  # type: ignore

import rmsd as rmsdlib


def test_quaternion_rmsd_pdb() -> None:

    filename_p = RESOURCE_PATH / "ci2_1.pdb"
    filename_q = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsdlib.get_coordinates_pdb(filename_p)
    _, q_coord = rmsdlib.get_coordinates_pdb(filename_q)

    p_center = rmsdlib.centroid(p_coord)
    q_center = rmsdlib.centroid(q_coord)
    p_coord -= p_center
    q_coord -= q_center

    result = rmsdlib.quaternion_rmsd(p_coord, q_coord)

    np.testing.assert_almost_equal(11.7768, result, decimal=3)


def test_quaternion_rotate_pdb() -> None:

    filename_p = RESOURCE_PATH / "ci2_1.pdb"
    filename_q = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsdlib.get_coordinates_pdb(filename_p)
    _, q_coord = rmsdlib.get_coordinates_pdb(filename_q)

    new_p_coord = rmsdlib.quaternion_rotate(p_coord, q_coord)

    np.testing.assert_array_almost_equal([-0.5124, 0.8565, 0.0608], new_p_coord[0], decimal=3)


def test_quaternion_transform() -> None:

    r = np.array([-0.31019, -0.59291, 0.63612, -0.38415])
    U = rmsdlib.quaternion_transform(r)

    np.testing.assert_array_almost_equal([-0.5124, 0.8565, 0.0608], U[0], decimal=3)


def test_makeQ() -> None:

    r = [-0.31019, -0.59291, 0.63612, -0.38415]
    Q_r = rmsdlib.makeQ(*r)

    np.testing.assert_array_almost_equal([-0.3841, -0.6361, -0.5929, -0.3101], Q_r[0], decimal=3)


def test_makeW() -> None:
    r = [-0.31019, -0.59291, 0.63612, -0.38415]
    Wt_r = rmsdlib.makeW(*r)

    np.testing.assert_array_almost_equal([-0.3841, 0.6361, 0.5929, -0.3101], Wt_r[0], decimal=3)
