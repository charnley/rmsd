import copy

import numpy as np
import pytest
from conftest import RESOURCE_PATH  # type: ignore

import rmsd as rmsdlib
from rmsd.calculate_rmsd import get_coordinates_xyz_lines


def test_print_reflection_reorder() -> None:
    # Test issue 78

    filename_a = RESOURCE_PATH / "issue78" / "a.xyz"
    filename_b = RESOURCE_PATH / "issue78" / "b.xyz"

    # Function call
    atoms_a, coord_a = rmsdlib.get_coordinates_xyz(filename_a, return_atoms_as_int=True)
    atoms_b, coord_b = rmsdlib.get_coordinates_xyz(filename_b, return_atoms_as_int=True)

    # re-center molecules
    coord_a -= rmsdlib.centroid(coord_a)
    coord_b -= rmsdlib.centroid(coord_b)

    rmsd_method = rmsdlib.kabsch_rmsd
    reorder_method = rmsdlib.reorder_hungarian

    result_rmsd, q_swap, q_reflection, q_review = rmsdlib.check_reflections(
        atoms_a,
        atoms_b,
        coord_a,
        coord_b,
        reorder_method=reorder_method,
        rmsd_method=rmsd_method,
    )

    print(result_rmsd)
    print(q_swap)
    print(q_reflection)
    print(q_review)

    # Apply the swap, reflection and review
    tmp_coord = copy.deepcopy(coord_b)
    tmp_coord = tmp_coord[:, q_swap]
    tmp_coord = np.dot(tmp_coord, np.diag(q_reflection))
    tmp_coord -= rmsdlib.centroid(tmp_coord)
    tmp_coord = tmp_coord[q_review]
    tmp_coord = rmsdlib.kabsch_rotate(tmp_coord, coord_a)
    rmsd_ = rmsdlib.rmsd(coord_a, tmp_coord)
    print(rmsd_)
    print(tmp_coord)

    # Main call rmsd value
    args = f"--use-reflections --reorder {filename_a} {filename_b}"
    print(args.split())
    value = float(rmsdlib.main(args.split()))
    print(value)
    assert value is not None
    np.testing.assert_almost_equal(result_rmsd, value)

    # Main call print, check rmsd is still the same
    # Note, that --print is translating b to a center
    _args = f"--use-reflections --reorder --print {filename_a} {filename_b}"
    _stdout: str = rmsdlib.main(_args.split())
    atoms, coord = rmsdlib.get_coordinates_xyz_lines(_stdout.split("\n"), return_atoms_as_int=True)
    coord -= rmsdlib.centroid(coord)  # fix translation
    print(coord)
    print(atoms)
    print(atoms_b)

    rmsd_check1 = rmsdlib.kabsch_rmsd(coord, coord_a)
    rmsd_check2 = rmsdlib.rmsd(coord, coord_a)
    print(rmsd_check1)
    print(rmsd_check2)
    print(result_rmsd)
    np.testing.assert_almost_equal(rmsd_check2, rmsd_check1)
    np.testing.assert_almost_equal(rmsd_check2, result_rmsd)


def test_bad_different_molcules() -> None:

    filename_a = RESOURCE_PATH / "ethane.xyz"
    filename_b = RESOURCE_PATH / "water.xyz"

    args = f"{filename_a} {filename_b}"

    with pytest.raises(SystemExit):
        rmsdlib.main(args.split())


def test_bad_different_order() -> None:

    filename_a = RESOURCE_PATH / "CHEMBL3039407.xyz"
    filename_b = RESOURCE_PATH / "CHEMBL3039407_order.xyz"

    args = f"{filename_a} {filename_b}"

    with pytest.raises(SystemExit):
        rmsdlib.main(args.split())


def test_rotation_methods() -> None:

    filename_a = RESOURCE_PATH / "ethane_translate.xyz"
    filename_b = RESOURCE_PATH / "ethane.xyz"

    rmsdlib.main(f"{filename_a} {filename_b} --rotation quaternion".split())

    rmsdlib.main(f"{filename_a} {filename_b} --rotation none".split())


def test_reorder_methods() -> None:

    filename_a = RESOURCE_PATH / "CHEMBL3039407.xyz"
    filename_b = RESOURCE_PATH / "CHEMBL3039407_order.xyz"

    methods = ["hungarian", "inertia-hungarian", "distance"]

    for method in methods:
        rmsdlib.main(f"--reorder --reorder-method {method} {filename_a} {filename_b}".split())

    rmsdlib.main(
        f"--no-hydrogen --reorder --reorder-method hungarian {filename_a} {filename_b}".split()
    )


def test_reflections() -> None:

    filename_a = RESOURCE_PATH / "CHEMBL3039407.xyz"
    filename_b = RESOURCE_PATH / "CHEMBL3039407.xyz"

    rmsdlib.main(f"--use-reflections {filename_a} {filename_b}".split())

    rmsdlib.main(f"--use-reflections-keep-stereo {filename_a} {filename_b}".split())


def test_ignore() -> None:

    filename_a = RESOURCE_PATH / "CHEMBL3039407.xyz"
    filename_b = RESOURCE_PATH / "CHEMBL3039407.xyz"

    rmsdlib.main(f"--no-hydrogen {filename_a} {filename_b}".split())

    rmsdlib.main(f"{filename_a} {filename_b} --remove-idx 0 5".split())

    rmsdlib.main(f"{filename_a} {filename_b} --add-idx 0 1 2 3 4".split())


def test_print_match_no_hydrogen() -> None:

    filename_a = RESOURCE_PATH / "CHEMBL3039407_order.xyz"
    filename_b = RESOURCE_PATH / "CHEMBL3039407_order.xyz"

    cmd = f"--no-hydrogen --print {filename_a} {filename_b}"
    print(cmd)
    out = rmsdlib.main(cmd.split()).split("\n")
    atoms1, coord1 = get_coordinates_xyz_lines(out)

    print(atoms1)
    print(len(atoms1))

    assert len(atoms1) == 60
    assert coord1.shape
    assert "H" in atoms1

    cmd = f"--print {filename_a} {filename_b}"
    out = rmsdlib.main(cmd.split()).split("\n")
    atoms2, coord2 = get_coordinates_xyz_lines(out)

    print(atoms2)
    print(len(atoms2))

    assert len(atoms2) == 60
    assert coord2.shape
    assert "H" in atoms2

    out = rmsdlib.main(
        f"--no-hydrogen --print --print-only-rmsd-atoms {filename_a} {filename_b}".split()
    ).split("\n")
    atoms1, coord1 = get_coordinates_xyz_lines(out)

    print(atoms1)
    print(len(atoms1))

    assert len(atoms1) == 30
    assert coord1.shape
    assert "H" not in atoms1


def test_only_alpha_carbons() -> None:

    filename_a = RESOURCE_PATH / "ci2_1.pdb"
    filename_b = RESOURCE_PATH / "ci2_2.pdb"

    cmd = f"--only-alpha-carbons --print {filename_a} {filename_b}"
    print(cmd)
    out = rmsdlib.main(cmd.split()).split("\n")
    atoms1, _ = get_coordinates_xyz_lines(out)

    print(atoms1)
    print(len(atoms1))

    assert len(atoms1) == 64
    assert len(np.unique(atoms1)) == 1
    assert set(atoms1) == set(["C"])

    cmd2 = f"{filename_a} {filename_b}"
    print(cmd)
    out2 = rmsdlib.main(cmd2.split())
    rmsd2 = float(out2)
    assert isinstance(rmsd2, float)

    cmd3 = f"--only-alpha-carbons {filename_a} {filename_b}"
    print(cmd)
    out3 = rmsdlib.main(cmd3.split())
    rmsd3 = float(out3)
    assert isinstance(rmsd3, float)

    print(rmsd2)
    print(rmsd3)
    assert rmsd2 > rmsd3
