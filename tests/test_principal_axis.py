import numpy as np

import rmsd as rmsdlib
from tests.conftest import RESOURCE_PATH  # type: ignore


def test_get_principal_axis_is_eigenvector() -> None:
    """Principal axis must satisfy I @ axis = lambda_max * axis."""
    filename = RESOURCE_PATH / "CHEMBL3039407.xyz"
    atoms, coord = rmsdlib.get_coordinates_xyz(filename, return_atoms_as_int=True)

    axis = rmsdlib.get_principal_axis(atoms, coord)

    inertia = rmsdlib.get_inertia_tensor(atoms, coord)
    eigval = np.linalg.eigvalsh(inertia)

    # Unit vector
    np.testing.assert_almost_equal(np.linalg.norm(axis), 1.0)

    # Eigen-equation residual (relative to lambda_max)
    residual = np.linalg.norm(inertia @ axis - eigval[-1] * axis) / eigval[-1]
    assert residual < 1e-10


def test_get_principal_axis_rotates_with_molecule() -> None:
    """Rotating the molecule must rotate the principal axis (up to sign)."""
    filename = RESOURCE_PATH / "CHEMBL3039407.xyz"
    atoms, coord = rmsdlib.get_coordinates_xyz(filename, return_atoms_as_int=True)

    axis_a = rmsdlib.get_principal_axis(atoms, coord)

    # Rotate 90 degrees about z: x -> y, y -> -x
    R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    axis_b = rmsdlib.get_principal_axis(atoms, coord @ R)

    # Eigenvector sign is arbitrary; compare up to sign
    np.testing.assert_almost_equal(np.abs(axis_a @ R @ axis_b), 1.0)
