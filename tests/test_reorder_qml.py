import copy

import numpy as np
import pytest
from conftest import RESOURCE_PATH  # type: ignore

import rmsd as rmsdlib

qmllib = pytest.importorskip("qmllib")


def test_reorder_qml() -> None:

    filename_1 = RESOURCE_PATH / "CHEMBL3039407.xyz"

    p_atoms, p_coord = rmsdlib.get_coordinates_xyz(filename_1, return_atoms_as_int=True)

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
    rotation_y = np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)],
        ]
    )

    q_coord = np.dot(q_coord, rotation_y)

    # Reorder with standard hungarian, this will fail reorder and give large
    # RMSD
    view_dist = rmsdlib.reorder_hungarian(p_atoms, q_atoms, p_coord, q_coord)
    q_atoms_dist = q_atoms[view_dist]
    q_coord_dist = q_coord[view_dist]
    _rmsd_dist = rmsdlib.kabsch_rmsd(p_coord, q_coord_dist)
    assert q_atoms_dist.tolist() == p_atoms.tolist()
    assert _rmsd_dist > 3.0

    # Reorder based in chemical similarity
    view = rmsdlib.reorder_similarity(p_atoms, q_atoms, p_coord, q_coord)
    q_atoms = q_atoms[view]
    q_coord = q_coord[view]

    # Calculate new RMSD with correct atom order
    _rmsd = rmsdlib.kabsch_rmsd(p_coord, q_coord)

    # Assert correct atom order
    assert q_atoms.tolist() == p_atoms.tolist()

    # Assert this is the same molecule
    pytest.approx(0.0) == _rmsd


def test_reorder_qml_distmat() -> None:

    filename_1 = RESOURCE_PATH / "CHEMBL3039407.xyz"

    p_atoms, p_coord = rmsdlib.get_coordinates_xyz(filename_1, return_atoms_as_int=True)

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
    rotation_y = np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [0, 1, 0],
            [-np.sin(theta), 0, np.cos(theta)],
        ]
    )

    q_coord = np.dot(q_coord, rotation_y)

    # Reorder based in chemical similarity
    view = rmsdlib.reorder_similarity(p_atoms, q_atoms, p_coord, q_coord, use_kernel=False)
    q_atoms = q_atoms[view]
    q_coord = q_coord[view]

    # Calculate new RMSD with correct atom order
    _rmsd = rmsdlib.kabsch_rmsd(p_coord, q_coord)

    # Assert correct atom order
    assert q_atoms.tolist() == p_atoms.tolist()

    # Assert this is the same molecule
    pytest.approx(0.0) == _rmsd


def test_pdb_only_carbon_possible() -> None:

    filename_a = RESOURCE_PATH / "issue98" / "test1.pdb"
    filename_b = RESOURCE_PATH / "issue98" / "test2.pdb"

    cmd = f"--reorder --reorder-method qml --only-alpha-carbons {filename_a} {filename_b}"
    out = rmsdlib.main(cmd.split())

    print(out)
    rmsd = float(out)
    assert isinstance(rmsd, float)
