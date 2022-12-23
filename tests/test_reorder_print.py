import numpy as np
from context import RESOURCE_PATH, call_main

import rmsd as rmsdlib
from rmsd import get_coordinates_xyz, get_coordinates_xyz_lines


def test_reorder_print_and_rmsd() -> None:

    # Issue 93, problem with printed structure after reorder.
    # - Calculate rmsd with --reorder (structure a and b)
    # - Calculate --print and --reorder to structure c.xyz
    # - Calculate normal rmsd a.xyz and c.xyz

    filename_a = RESOURCE_PATH / "issue93" / "a.xyz"
    filename_b = RESOURCE_PATH / "issue93" / "b.xyz"

    # Get reorder rmsd
    args = ["--reorder", f"{filename_a}", f"{filename_b}"]
    stdout = call_main(args)
    rmsd_ab = float(stdout[-1])
    print(rmsd_ab)

    # Get printed structure
    stdout = call_main(args + ["--print"])

    atoms_a, coord_a = get_coordinates_xyz(filename_a)
    atoms_c, coord_c = get_coordinates_xyz_lines(stdout)

    print(coord_a)
    print(atoms_a)

    print(coord_c)
    print(atoms_c)

    rmsd_ac = rmsdlib.rmsd(coord_a, coord_c)
    print(rmsd_ac)

    np.testing.assert_almost_equal(rmsd_ab, rmsd_ac, decimal=8)
