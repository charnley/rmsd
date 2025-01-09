import numpy as np
from conftest import RESOURCE_PATH  # type: ignore

import rmsd as rmsdlib
from rmsd import get_coordinates_xyz, get_coordinates_xyz_lines


def test_reorder_print_and_rmsd() -> None:

    # Issue 93, problem with printed structure after reorder.
    # - Calculate rmsd with --reorder (structure a and b)
    # - Calculate --print and --reorder to structure c.xyz
    # - Calculate normal rmsd a.xyz and c.xyz

    filename_a = RESOURCE_PATH / "issue93" / "a.xyz"
    filename_b = RESOURCE_PATH / "issue93" / "b.xyz"
    atoms_a, coord_a = get_coordinates_xyz(filename_a)
    atoms_b, coord_b = get_coordinates_xyz(filename_b)

    # Get reorder rmsd
    rmsd_ab = float(rmsdlib.main(f"--reorder {filename_a} {filename_b}".split()))
    print(rmsd_ab)
    assert isinstance(rmsd_ab, float)

    # Get printed structure
    stdout = rmsdlib.main(f"--reorder --print {filename_a} {filename_b}".split())
    print(stdout)
    atoms_c, coord_c = get_coordinates_xyz_lines(stdout.split("\n"))

    coord_c -= rmsdlib.centroid(coord_c)
    coord_a -= rmsdlib.centroid(coord_a)

    print(atoms_a)
    print(atoms_b)
    print(atoms_c)

    print(coord_a)
    print(coord_b)
    print(coord_c)

    assert (atoms_a == atoms_c).all()
    assert (atoms_a != atoms_b).any()

    rmsd_ac = rmsdlib.rmsd(coord_a, coord_c)
    print(rmsd_ac)

    np.testing.assert_almost_equal(rmsd_ab, rmsd_ac, decimal=8)
