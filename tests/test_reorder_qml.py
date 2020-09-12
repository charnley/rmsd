
import pathlib
import copy
import numpy as np
import pytest

from constants import RESOURCE_PATH

qml = pytest.importorskip("qml")

import rmsd


def test_reorder_qml():

    filename_1 = pathlib.PurePath(
        RESOURCE_PATH,
        'CHEMBL3039407.xyz'
    )

    p_atoms, p_coord = rmsd.get_coordinates_xyz(filename_1)

    # Reorder atoms
    n_atoms = len(p_atoms)
    random_reorder = np.arange(n_atoms, dtype=int)
    np.random.seed(5)
    np.random.shuffle(random_reorder)

    q_atoms = copy.deepcopy(p_atoms)
    q_coord = copy.deepcopy(p_coord)
    q_atoms = q_atoms[random_reorder]
    q_coord = q_coord[random_reorder]

    # Mess up the distance matrix by rotating the molecule
    theta = 180.0
    rotation_y = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

    q_coord = np.dot(q_coord, rotation_y)

    # Reorder with standard hungarian, this will fail reorder and give large
    # RMSD
    view_dist = rmsd.reorder_hungarian(
        p_atoms,
        q_atoms,
        p_coord,
        q_coord
    )
    q_atoms_dist = q_atoms[view_dist]
    q_coord_dist = q_coord[view_dist]
    _rmsd_dist = rmsd.kabsch_rmsd(p_coord, q_coord_dist)
    assert q_atoms_dist.tolist() == p_atoms.tolist()
    assert _rmsd_dist > 3.0

    # Reorder based in chemical similarity
    view = rmsd.reorder_similarity(
        p_atoms,
        q_atoms,
        p_coord,
        q_coord
    )
    q_atoms = q_atoms[view]
    q_coord = q_coord[view]

    # Calculate new RMSD with correct atom order
    _rmsd = rmsd.kabsch_rmsd(p_coord, q_coord)

    # Assert correct atom order
    assert q_atoms.tolist() == p_atoms.tolist()

    # Assert this is the same molecule
    pytest.approx(0.0) == _rmsd


def test_reorder_qml_distmat():

    filename_1 = pathlib.PurePath(
        RESOURCE_PATH,
        'CHEMBL3039407.xyz'
    )

    p_atoms, p_coord = rmsd.get_coordinates_xyz(filename_1)

    # Reorder atoms
    n_atoms = len(p_atoms)
    random_reorder = np.arange(n_atoms, dtype=int)
    np.random.seed(5)
    np.random.shuffle(random_reorder)

    q_atoms = copy.deepcopy(p_atoms)
    q_coord = copy.deepcopy(p_coord)
    q_atoms = q_atoms[random_reorder]
    q_coord = q_coord[random_reorder]

    # Mess up the distance matrix by rotating the molecule
    theta = 180.0
    rotation_y = np.array([
        [np.cos(theta), 0, np.sin(theta)],
        [0, 1, 0],
        [-np.sin(theta), 0, np.cos(theta)]
    ])

    q_coord = np.dot(q_coord, rotation_y)

    # Reorder based in chemical similarity
    view = rmsd.reorder_similarity(
        p_atoms,
        q_atoms,
        p_coord,
        q_coord,
        use_kernel=False
    )
    q_atoms = q_atoms[view]
    q_coord = q_coord[view]

    # Calculate new RMSD with correct atom order
    _rmsd = rmsd.kabsch_rmsd(p_coord, q_coord)

    # Assert correct atom order
    assert q_atoms.tolist() == p_atoms.tolist()

    # Assert this is the same molecule
    pytest.approx(0.0) == _rmsd
