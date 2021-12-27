import pytest

from rmsd import calculate_rmsd


def test_legal_arguments() -> None:

    args_ = "--rotation kabsch --ignore-hydrogen FILE_A FILE_B".split()
    args = calculate_rmsd.parse_arguments(args_)

    assert args.reorder is False
    assert args.ignore_hydrogen is True
    assert args.rotation == "kabsch"


def test_illegal_arguments() -> None:

    args = [
        "--rotation kabsch",
        "--reorder",
        "--print",
        "--ignore-hydrogen",
        "FILE_A",
        "FILE_B",
    ]

    with pytest.raises(SystemExit) as exception:
        _ = calculate_rmsd.parse_arguments(args)

    assert exception.type == SystemExit


def test_illegal_reflection() -> None:

    args = [
        "--rotation kabsch",
        "--use-reflections",
        "--print",
        "--ignore-hydrogen",
        "FILE_A",
        "FILE_B",
    ]

    with pytest.raises(SystemExit) as exception:
        _ = calculate_rmsd.parse_arguments(args)

    assert exception.type == SystemExit


def test_illegal_rotation_method() -> None:

    args = ["--rotation NeverHeardOfThisMethod", "FILE_A", "FILE_B"]

    with pytest.raises(SystemExit) as exception:
        _ = calculate_rmsd.parse_arguments(args)

    assert exception.type == SystemExit


def test_illegal_reorder_method() -> None:

    args = ["--reorder-method NotImplementedYet", "FILE_A", "FILE_B"]

    with pytest.raises(SystemExit) as exception:
        _ = calculate_rmsd.parse_arguments(args)

    assert exception.type == SystemExit
