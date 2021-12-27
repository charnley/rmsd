import gzip
from pathlib import Path

import numpy as np
import pytest
from context import RESOURCE_PATH

import rmsd


def test_get_coordinates_pdb() -> None:

    filename = RESOURCE_PATH / "ci2_1.pdb"
    atoms, coords = rmsd.get_coordinates_pdb(filename)
    assert "N" == atoms[0]
    assert [-7.173, -13.891, -6.266] == coords[0].tolist()

    filename = RESOURCE_PATH / "ci2_1.pdb"
    atoms, coords = rmsd.get_coordinates(filename, "pdb")
    assert "N" == atoms[0]
    assert [-7.173, -13.891, -6.266] == coords[0].tolist()


def test_get_coordinates_wrong() -> None:
    filename = RESOURCE_PATH / "ci2_1.pdb"
    with pytest.raises(ValueError):
        _, _ = rmsd.get_coordinates(filename, "qqq")


def test_get_coordinates_xyz() -> None:

    filename = RESOURCE_PATH / "ethane.xyz"
    atoms, coords = rmsd.get_coordinates_xyz(filename)

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_get_coordinates_xyz_bad() -> None:
    filename = RESOURCE_PATH / "bad_format.xyz"
    _, _ = rmsd.get_coordinates(filename, "xyz")


def test_get_coordinates() -> None:

    filename = RESOURCE_PATH / "ethane.xyz"
    atoms, coords = rmsd.get_coordinates(filename, "xyz")

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_get_coordinates_gzip(tmp_path: Path) -> None:

    filename = RESOURCE_PATH / "ethane.xyz"
    with open(filename, "r") as f:
        content = f.read()

    content_byte = content.encode()

    filename_gzip = tmp_path / "ethane.xyz.gz"

    with gzip.open(filename_gzip, "wb") as f:  # type: ignore
        f.write(content_byte)  # type: ignore

    atoms, coords = rmsd.get_coordinates(filename_gzip, "xyz", is_gzip=True)

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_rmsd_pdb() -> None:

    filename_1 = RESOURCE_PATH / "ci2_1.pdb"
    filename_2 = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsd.get_coordinates_pdb(filename_1)
    _, q_coord = rmsd.get_coordinates_pdb(filename_2)

    pure_rmsd = rmsd.rmsd(p_coord, q_coord)

    np.testing.assert_almost_equal(26.9750, pure_rmsd, decimal=3)


def test_rmsd_xyz() -> None:

    filename_1 = RESOURCE_PATH / "ethane.xyz"
    filename_2 = RESOURCE_PATH / "ethane_mini.xyz"

    _, p_coord = rmsd.get_coordinates_xyz(filename_1)
    _, q_coord = rmsd.get_coordinates_xyz(filename_2)

    pure_rmsd = rmsd.rmsd(p_coord, q_coord)

    np.testing.assert_almost_equal(0.33512, pure_rmsd, decimal=3)
