import numpy as np

import rmsd as rmsdlib
from tests.conftest import RESOURCE_PATH, rotate_coord  # type: ignore


def test_reorder_inertia_hungarian_butane() -> None:

    filename_a = RESOURCE_PATH / "butane.xyz"
    filename_b = RESOURCE_PATH / "butane_prime.xyz"

    atoms_a, coord_a = rmsdlib.get_coordinates_xyz(filename_a, return_atoms_as_int=True)
    atoms_b, coord_b = rmsdlib.get_coordinates_xyz(filename_b, return_atoms_as_int=True)

    coord_a -= rmsdlib.centroid(coord_a)
    coord_b -= rmsdlib.centroid(coord_b)

    review = rmsdlib.reorder_inertia_hungarian(atoms_a, atoms_b, coord_a, coord_b)
    print(review)

    result_rmsd = rmsdlib.kabsch_rmsd(coord_a, coord_b[review])

    np.testing.assert_almost_equal(0, result_rmsd, decimal=2)


def test_reorder_inertia_hungarian_complicated() -> None:

    filename_a = RESOURCE_PATH / "CHEMBL3039407.xyz"

    atoms_a, coord_a = rmsdlib.get_coordinates_xyz(filename_a, return_atoms_as_int=True)
    atoms_b, coord_b = rmsdlib.get_coordinates_xyz(filename_a, return_atoms_as_int=True)

    coord_a -= rmsdlib.centroid(coord_a)
    coord_b -= rmsdlib.centroid(coord_b)

    # Setup the problem
    # Rotate a 30 degrees in xy and 30 degrees in yz planes
    coord_a = rotate_coord(30, coord_a)
    coord_a = rotate_coord(30, coord_a, axis=[1, 2])

    # Randomize the atom order
    ind = np.array(list(range(len(atoms_b))))
    rng = np.random.default_rng(seed=56)
    rng.shuffle(ind)
    atoms_b = atoms_b[ind]
    coord_b = coord_b[ind]

    # Solve the problem
    review = rmsdlib.reorder_inertia_hungarian(atoms_a, atoms_b, coord_a, coord_b)
    print(review)

    result_rmsd = rmsdlib.kabsch_rmsd(coord_a, coord_b[review])
    print(result_rmsd)

    np.testing.assert_almost_equal(0.0, result_rmsd)


def test_reorder_and_distance() -> None:

    filename_a = RESOURCE_PATH / "CHEMBL3039407.xyz"

    atoms_a, coord_a = rmsdlib.get_coordinates_xyz(filename_a, return_atoms_as_int=True)
    atoms_b, coord_b = rmsdlib.get_coordinates_xyz(filename_a, return_atoms_as_int=True)

    coord_a -= rmsdlib.centroid(coord_a)
    coord_b -= rmsdlib.centroid(coord_b)

    # Setup the problem
    # Rotate a 30 degrees in xy and 30 degrees in yz planes
    coord_a = rotate_coord(30, coord_a)
    coord_a = rotate_coord(30, coord_a, axis=[1, 2])

    # Make coord_b have bigger bonds
    coord_b = coord_b * 1.2
    answer = 0.8035106793656143

    print(coord_b)

    # Randomize the atom order
    ind = np.array(list(range(len(atoms_b))))
    rng = np.random.default_rng(seed=56)
    rng.shuffle(ind)
    atoms_b = atoms_b[ind]
    coord_b = coord_b[ind]

    # Solve the problem
    review = rmsdlib.reorder_inertia_hungarian(atoms_a, atoms_b, coord_a, coord_b)
    print(review)

    result_rmsd = rmsdlib.kabsch_rmsd(coord_a, coord_b[review], translate=True)
    print(result_rmsd)

    np.testing.assert_almost_equal(answer, result_rmsd, decimal=3)
