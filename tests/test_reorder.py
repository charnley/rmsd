import copy

import numpy as np
from numpy import ndarray

import rmsd as rmsdlib


def test_reorder_distance() -> None:

    N = 5
    atoms = np.array(["H"] * N)
    p_coord: ndarray = np.arange(N * 3)
    p_coord = p_coord.reshape((5, 3))
    q_coord = np.array(p_coord, copy=True)

    np.random.seed(6)
    np.random.shuffle(q_coord)

    review = rmsdlib.reorder_hungarian(atoms, atoms, p_coord, q_coord)

    assert p_coord.tolist() == q_coord[review].tolist()


def test_reorder_brute() -> None:
    N = 5
    atoms = np.array(["H"] * N)
    p_coord: ndarray = np.arange(N * 3)
    p_coord = p_coord.reshape((N, 3))
    q_coord = np.array(p_coord, copy=True)

    np.random.seed(6)
    np.random.shuffle(q_coord)

    review = rmsdlib.reorder_brute(atoms, atoms, p_coord, q_coord)
    assert p_coord.tolist() == q_coord[review].tolist()


def test_reorder_brute_ch() -> None:

    N = 6
    p_atoms_str = ["C"] * 3 + ["H"] * 3
    p_atoms_int = [rmsdlib.int_atom(atom) for atom in p_atoms_str]
    p_atoms = np.array(p_atoms_int)
    p_coord: ndarray = np.arange(N * 3)
    p_coord = p_coord.reshape((N, 3))

    # random index
    np.random.seed(6)
    idx = np.arange(N, dtype=int)
    np.random.shuffle(idx)

    q_coord = copy.deepcopy(p_coord)
    q_atoms = copy.deepcopy(p_atoms)

    q_coord = q_coord[idx]
    q_atoms = q_atoms[idx]

    review = rmsdlib.reorder_brute(p_atoms, q_atoms, p_coord, q_coord)

    assert p_coord.tolist() == q_coord[review].tolist()
    assert p_atoms.tolist() == q_atoms[review].tolist()


def test_reorder_hungarian() -> None:

    N = 5
    atoms = np.array(["H"] * N)
    p_coord: ndarray = np.arange(N * 3)
    p_coord = p_coord.reshape((5, 3))
    q_coord = copy.deepcopy(p_coord)

    np.random.seed(6)
    np.random.shuffle(q_coord)

    review = rmsdlib.reorder_distance(atoms, atoms, p_coord, q_coord)
    assert p_coord.tolist() == q_coord[review].tolist()


def test_reorder_inertia_hungarian() -> None:

    # coordinates of scrambled and rotated butane
    atoms = np.array(["C", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"])
    atoms_ = [rmsdlib.int_atom(atom) for atom in atoms]
    atoms = np.array(atoms_)

    p_coord = np.array(
        [
            [2.142e00, 1.395e00, -8.932e00],
            [3.631e00, 1.416e00, -8.537e00],
            [4.203e00, -1.200e-02, -8.612e00],
            [5.691e00, 9.000e-03, -8.218e00],
            [1.604e00, 7.600e-01, -8.260e00],
            [1.745e00, 2.388e00, -8.880e00],
            [2.043e00, 1.024e00, -9.930e00],
            [4.169e00, 2.051e00, -9.210e00],
            [3.731e00, 1.788e00, -7.539e00],
            [3.665e00, -6.470e-01, -7.940e00],
            [4.104e00, -3.840e-01, -9.610e00],
            [6.088e00, -9.830e-01, -8.270e00],
            [5.791e00, 3.810e-01, -7.220e00],
            [6.230e00, 6.440e-01, -8.890e00],
        ]
    )

    q_coord = np.array(
        [
            [6.71454, -5.53848, -3.50851],
            [6.95865, -6.22697, -2.15264],
            [8.16747, -5.57632, -1.45606],
            [5.50518, -6.19016, -4.20589],
            [5.33617, -5.71137, -5.14853],
            [7.58263, -5.64795, -4.12498],
            [6.51851, -4.49883, -3.35011],
            [6.09092, -6.11832, -1.53660],
            [5.70232, -7.22908, -4.36475],
            [7.15558, -7.26640, -2.31068],
            [8.33668, -6.05459, -0.51425],
            [7.97144, -4.53667, -1.29765],
            [4.63745, -6.08152, -3.58986],
            [9.03610, -5.68475, -2.07173],
        ]
    )

    p_coord -= rmsdlib.centroid(p_coord)
    q_coord -= rmsdlib.centroid(q_coord)

    review = rmsdlib.reorder_inertia_hungarian(atoms, atoms, p_coord, q_coord)

    result_rmsd = rmsdlib.kabsch_rmsd(p_coord, q_coord[review])

    np.testing.assert_almost_equal(0, result_rmsd, decimal=2)
