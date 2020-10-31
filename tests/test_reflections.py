import copy

import numpy as np

import rmsd


def test_reflections():

    atoms = np.array(["C", "H", "H", "H", "F"])

    p_coord = np.array(
        [
            [-0.000000, -0.000000, -0.000000],
            [1.109398, -0.000000, 0.000000],
            [-0.3697920, -0.7362220, -0.7429600],
            [-0.3698020, 1.011538, -0.2661100],
            [-0.3698020, -0.2753120, 1.009070],
        ]
    )

    q_coord = copy.deepcopy(p_coord)

    # insert reflection
    # TODO Insert a rotation on q
    q_coord[:, [0, 2]] = q_coord[:, [2, 0]]

    min_rmsd, min_swap, min_reflection, min_review = rmsd.check_reflections(
        atoms, atoms, p_coord, q_coord, reorder_method=None
    )

    assert np.isclose(min_rmsd, 0.0, atol=1e-6)


def test_reflections_norotation():

    atoms = np.array(["C", "H", "H", "H", "F"])

    p_coord = np.array(
        [
            [-0.000000, -0.000000, -0.000000],
            [1.109398, -0.000000, 0.000000],
            [-0.3697920, -0.7362220, -0.7429600],
            [-0.3698020, 1.011538, -0.2661100],
            [-0.3698020, -0.2753120, 1.009070],
        ]
    )

    q_coord = copy.deepcopy(p_coord)

    # Insert reflection
    q_coord[:, [0, 2]] = q_coord[:, [2, 0]]

    min_rmsd, min_swap, min_reflection, min_review = rmsd.check_reflections(
        atoms,
        atoms,
        p_coord,
        q_coord,
        reorder_method=None,
        rotation_method=None,
    )

    assert np.isclose(min_rmsd, 0.0, atol=1e-6)


def test_reflections_reorder():

    p_atoms = np.array(["C", "H", "H", "H", "F"])

    p_coord = np.array(
        [
            [-0.000000, -0.000000, -0.000000],
            [1.109398, -0.000000, 0.000000],
            [-0.3697920, -0.7362220, -0.7429600],
            [-0.3698020, 1.011538, -0.2661100],
            [-0.3698020, -0.2753120, 1.009070],
        ]
    )

    # random reflection
    q_coord = copy.deepcopy(p_coord)
    q_coord[:, [0, 2]] = q_coord[:, [2, 0]]

    # random order
    # review = list(range(len(atoms)))
    review = [1, 4, 2, 0, 3]
    q_coord = q_coord[review]
    q_atoms = p_atoms[review]

    min_rmsd, min_swap, min_reflection, min_review = rmsd.check_reflections(
        p_atoms, q_atoms, p_coord, q_coord
    )

    assert np.isclose(min_rmsd, 0.0, atol=1e-6)


def test_reflections_keep_stereo():

    atoms = np.array(["C", "H", "H", "H", "F"])

    p_coord = np.array(
        [
            [-0.000000, -0.000000, -0.000000],
            [1.109398, -0.000000, 0.000000],
            [-0.3697920, -0.7362220, -0.7429600],
            [-0.3698020, 1.011538, -0.2661100],
            [-0.3698020, -0.2753120, 1.009070],
        ]
    )

    q_coord = copy.deepcopy(p_coord)

    # swap [2,1,0]: prepare enantiomer coordinates (which is named as q_coord)
    # of p_coord
    q_coord[:, [0, 2]] = q_coord[:, [2, 0]]

    # If keep_stereo is off, enantiomer coordinates of q_coord are considered,
    # resulting into identical coordinates of p_coord.
    min_rmsd, min_swap, min_reflection, min_review = rmsd.check_reflections(
        atoms, atoms, p_coord, q_coord, reorder_method=None, keep_stereo=False
    )

    # the enantiomer of the enantiomer of the original molecule has zero RMSD
    # with the original molecule.
    assert np.isclose(min_rmsd, 0.0, atol=1e-6)

    # No enantiomer coordinates, non-zero RMSD.
    min_rmsd, min_swap, min_reflection, min_review = rmsd.check_reflections(
        atoms, atoms, p_coord, q_coord, reorder_method=None, keep_stereo=True
    )

    assert np.isclose(min_rmsd, 1.1457797, atol=1e-6)
