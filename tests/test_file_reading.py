import gzip
import pathlib

import numpy as np
from context import RESOURCE_PATH

import rmsd


def test_get_coordinates_pdb():

    filename = pathlib.PurePath(RESOURCE_PATH, "ci2_1.pdb")
    atoms, coords = rmsd.get_coordinates_pdb(filename)

    assert "N" == atoms[0]
    assert [-7.173, -13.891, -6.266] == coords[0].tolist()


def test_get_coordinates_xyz():

    filename = pathlib.PurePath(RESOURCE_PATH, "ethane.xyz")
    atoms, coords = rmsd.get_coordinates_xyz(filename)

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_get_coordinates():

    filename = pathlib.PurePath(RESOURCE_PATH, "ethane.xyz")
    atoms, coords = rmsd.get_coordinates(filename, "xyz")

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_get_coordinates_gzip(tmpdir):

    filename = pathlib.PurePath(RESOURCE_PATH, "ethane.xyz")
    with open(filename, "r") as f:
        content = f.read()

    content = content.encode()

    filename_gzip = tmpdir.join("ethane.xyz.gz")

    with gzip.open(filename_gzip, "wb") as f:
        f.write(content)

    atoms, coords = rmsd.get_coordinates(filename_gzip, "xyz", is_gzip=True)

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_rmsd_pdb():

    filename_1 = pathlib.PurePath(RESOURCE_PATH, "ci2_1.pdb")
    filename_2 = pathlib.PurePath(RESOURCE_PATH, "ci2_2.pdb")

    p_atoms, p_coord = rmsd.get_coordinates_pdb(filename_1)
    q_atoms, q_coord = rmsd.get_coordinates_pdb(filename_2)

    pure_rmsd = rmsd.rmsd(p_coord, q_coord)

    np.testing.assert_almost_equal(26.9750, pure_rmsd, decimal=3)


def test_rmsd_xyz():

    filename_1 = pathlib.PurePath(RESOURCE_PATH, "ethane.xyz")
    filename_2 = pathlib.PurePath(RESOURCE_PATH, "ethane_mini.xyz")

    p_atoms, p_coord = rmsd.get_coordinates_xyz(filename_1)
    q_atoms, q_coord = rmsd.get_coordinates_xyz(filename_2)

    pure_rmsd = rmsd.rmsd(p_coord, q_coord)

    np.testing.assert_almost_equal(0.33512, pure_rmsd, decimal=3)
