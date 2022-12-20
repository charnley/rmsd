import copy

import numpy as np
import pytest
from context import RESOURCE_PATH, call_main

import rmsd as rmsdlib


def test_print_reflection_reorder() -> None:
    # Test issue 78

    filename_a = RESOURCE_PATH / "issue78" / "a.xyz"
    filename_b = RESOURCE_PATH / "issue78" / "b.xyz"

    # Function call
    atoms_a, coord_a = rmsdlib.get_coordinates_xyz(filename_a, return_atoms_as_int=True)
    atoms_b, coord_b = rmsdlib.get_coordinates_xyz(filename_b, return_atoms_as_int=True)

    n_atoms = len(atoms_a)

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

    # tmp_atoms = copy.copy(atoms_b)
    tmp_coord = copy.deepcopy(coord_b)
    tmp_coord = tmp_coord[:, q_swap]
    tmp_coord = np.dot(tmp_coord, np.diag(q_reflection))
    tmp_coord -= rmsdlib.centroid(tmp_coord)
    tmp_coord = tmp_coord[q_review]
    tmp_coord = rmsdlib.kabsch_rotate(tmp_coord, coord_a)
    rmsd_ = rmsdlib.rmsd(coord_a, tmp_coord)
    print(rmsd_)

    # Main call rmsd value
    args = f"--use-reflections --reorder {filename_a} {filename_b}"
    stdout = call_main(args.split())
    value = float(stdout[-1])
    assert value is not None
    np.testing.assert_almost_equal(result_rmsd, value)

    # Main call print, check rmsd is still the same
    args = f"--use-reflections --reorder --print {filename_a} {filename_b}"
    stdout = call_main(args.split())
    _, coord = rmsdlib.get_coordinates_xyz_lines(stdout[-(n_atoms + 2) :])
    rmsd_check = rmsdlib.kabsch_rmsd(coord, coord_a, translate=True)
    np.testing.assert_almost_equal(result_rmsd, rmsd_check)


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
