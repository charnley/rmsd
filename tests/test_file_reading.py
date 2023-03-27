import gzip
from pathlib import Path

import numpy as np
import pytest
from context import RESOURCE_PATH

import rmsd as rmsdlib


def test_get_coordinates_pdb_hetatm() -> None:

    filename = Path(RESOURCE_PATH) / "issue88" / "native.pdb"
    atoms, _ = rmsdlib.get_coordinates_pdb(filename)

    assert len(atoms)
    assert atoms[0] == "C"
    assert atoms[-1] == "C"


def test_get_coordinates_pdb() -> None:

    filename = RESOURCE_PATH / "ci2_1.pdb"
    atoms, coords = rmsdlib.get_coordinates_pdb(filename)
    assert "N" == atoms[0]
    assert [-7.173, -13.891, -6.266] == coords[0].tolist()


def test_get_coordinates_wrong() -> None:
    filename = RESOURCE_PATH / "ci2_1.pdb"
    with pytest.raises(ValueError):
        _, _ = rmsdlib.get_coordinates(filename, "qqq")


def test_get_coordinates_xyz() -> None:

    filename = RESOURCE_PATH / "ethane.xyz"
    atoms, coords = rmsdlib.get_coordinates_xyz(filename)

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_get_coordinates_xyz_bad() -> None:
    filename = RESOURCE_PATH / "bad_format.xyz"
    _, _ = rmsdlib.get_coordinates(filename, "xyz")


def test_get_coordinates() -> None:

    filename = RESOURCE_PATH / "ethane.xyz"
    atoms, coords = rmsdlib.get_coordinates(filename, "xyz")

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_get_coordinates_gzip(tmp_path: Path) -> None:

    filename = RESOURCE_PATH / "ci2_1.pdb"
    with open(filename, "r") as f:
        content = f.read()

    content_byte = content.encode()

    filename_gzip = tmp_path / "molecule.pdb.gz"

    with gzip.open(filename_gzip, "wb") as f:  # type: ignore
        f.write(content_byte)  # type: ignore

    atoms, coords = rmsdlib.get_coordinates(
        filename_gzip, "pdb", is_gzip=True, return_atoms_as_int=True
    )
    assert 7 == atoms[0]
    assert [-7.173, -13.891, -6.266] == coords[0].tolist()


def test_get_coordinates_gzip_pdb(tmp_path: Path) -> None:

    filename = RESOURCE_PATH / "ethane.xyz"
    with open(filename, "r") as f:
        content = f.read()

    content_byte = content.encode()

    filename_gzip = tmp_path / "ethane.xyz.gz"

    with gzip.open(filename_gzip, "wb") as f:  # type: ignore
        f.write(content_byte)  # type: ignore

    atoms, coords = rmsdlib.get_coordinates(filename_gzip, "xyz", is_gzip=True)

    assert "C" == atoms[0]
    assert [-0.98353, 1.81095, -0.0314] == coords[0].tolist()


def test_rmsd_pdb() -> None:

    filename_1 = RESOURCE_PATH / "ci2_1.pdb"
    filename_2 = RESOURCE_PATH / "ci2_2.pdb"

    _, p_coord = rmsdlib.get_coordinates_pdb(filename_1)
    _, q_coord = rmsdlib.get_coordinates_pdb(filename_2)

    pure_rmsd = rmsdlib.rmsd(p_coord, q_coord)

    np.testing.assert_almost_equal(26.9750, pure_rmsd, decimal=3)


def test_rmsd_xyz() -> None:

    filename_1 = RESOURCE_PATH / "ethane.xyz"
    filename_2 = RESOURCE_PATH / "ethane_mini.xyz"

    _, p_coord = rmsdlib.get_coordinates_xyz(filename_1)
    _, q_coord = rmsdlib.get_coordinates_xyz(filename_2)

    pure_rmsd = rmsdlib.rmsd(p_coord, q_coord)

    np.testing.assert_almost_equal(0.33512, pure_rmsd, decimal=3)


def test_rmsd_pdb_only_alpha() -> None:
    """Read only alpha-carbon in a PDB file

    Example, one two pdb files has same number of residues but one has
    hydrogens and side chains
    """

    filename_1 = RESOURCE_PATH / "issue98" / "test1.pdb"
    filename_2 = RESOURCE_PATH / "issue98" / "test2.pdb"

    kwargs = dict(
        return_atoms_as_int=True,
        only_alpha_carbon=True,
    )
    p_atoms, p_coord = rmsdlib.get_coordinates_pdb(filename_1, **kwargs)
    q_atoms, q_coord = rmsdlib.get_coordinates_pdb(filename_2, **kwargs)

    # only carbons
    assert len(np.unique(p_atoms)) == 1
    assert len(np.unique(q_atoms)) == 1
    assert q_atoms[0] == 6  # only Carbon

    # Always translate
    p_cent = rmsdlib.centroid(p_coord)
    q_cent = rmsdlib.centroid(q_coord)
    p_coord -= p_cent
    q_coord -= q_cent

    rmsd = rmsdlib.kabsch_rmsd(p_coord, q_coord)
    print(rmsd)

    np.testing.assert_almost_equal(1.94, rmsd, decimal=3)
