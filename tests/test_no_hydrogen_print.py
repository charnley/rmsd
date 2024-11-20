import numpy as np
from context import RESOURCE_PATH, call_main

from rmsd import get_coordinates_xyz, get_coordinates_xyz_lines


def get_printed_structure(filename_a, filename_b):
    args = ["--no-hydrogen", f"{filename_a}", f"{filename_b}"]
    stdout = call_main(args)
    print("RMSD between %s and %s" % (filename_a, filename_b))
    rmsd_ab = float(stdout[-1])
    print(rmsd_ab)

    stdout = call_main(args + ["--print"])

    atoms_a, coord_a = get_coordinates_xyz(filename_a)
    atoms_c, coord_c = get_coordinates_xyz_lines(stdout)

    print(filename_a)
    print(coord_a)
    print(atoms_a)

    print("Rotated", filename_b)
    print(coord_c)
    print(atoms_c)
    print()

    return rmsd_ab, coord_c, atoms_c


def test_reorder_print_and_rmsd() -> None:

    # Pull Request 111, problem with printed structure after selection of atoms
    # with --no-hydrogens.
    # Here a pair of ethane structures with hydrogens and withoud hydrogens are
    # used to test the problem.
    # The following steps are performed and results are compared:
    # - calculate_rmsd --no-hydrogen --print (structure a1 and b1 with hydrogens)
    # - calculate_rmsd --no-hydrogen --print (structure a2 and b2 without hydrogens)

    # Ethane with hydrogens
    filename_a1 = RESOURCE_PATH / "ethane.xyz"
    filename_b1 = RESOURCE_PATH / "ethane_translate.xyz"

    # Ethane does not contain hydrogen, but the coordinates of the carbon atoms
    # are the same as in the previous structure.
    filename_a2 = RESOURCE_PATH / "pr111" / "ethane_no_hydrogen.xyz"
    filename_b2 = RESOURCE_PATH / "pr111" / "ethane_translate_no_hydrogen.xyz"

    rmsd_a1b1, coord_c1, atoms_c1 = get_printed_structure(filename_a1, filename_b1)
    rmsd_a2b2, coord_c2, atoms_c2 = get_printed_structure(filename_a2, filename_b2)

    # Compare results
    # The RMSD between the structures with and without hydrogens should be the same.
    np.testing.assert_almost_equal(rmsd_a1b1, rmsd_a2b2, decimal=8)
    # The coordinates of the structures should be the same.
    np.testing.assert_almost_equal(coord_c1, coord_c2, decimal=8)
    # The atoms of the structures should be the same.
    np.testing.assert_array_equal(atoms_c1, atoms_c2, err_msg="Atoms are not the same!")

    print("*" * 10)
    print("Test passed!")
    print("*" * 10)


if __name__ == "__main__":
    test_reorder_print_and_rmsd()
