import pytest

from rmsd import calculate_rmsd


def test_formats() -> None:

    args_ = "filename.xyz.gz filename2.xyz.gz".split()
    args = calculate_rmsd.parse_arguments(args_)

    assert args.format_is_gzip


def test_legal_arguments() -> None:

    args_ = "--rotation kabsch --ignore-hydrogen --format xyz FILE_A FILE_B".split()
    args = calculate_rmsd.parse_arguments(args_)

    assert args.reorder is False
    assert args.ignore_hydrogen is True
    assert args.rotation == "kabsch"


def test_illegal_arguments() -> None:

    with pytest.raises(SystemExit):
        args = calculate_rmsd.parse_arguments(
            "--reorder --ignore-hydrogen --print filea fileb".split()
        )
        print(args)

    with pytest.raises(SystemExit):
        args = calculate_rmsd.parse_arguments(
            "--print --ignore-hydrogen --use-reflections filea fileb".split()
        )
        print(args)

    with pytest.raises(SystemExit):
        args = calculate_rmsd.parse_arguments("--rotation do-not-exists filea fileb".split())
        print(args)

    with pytest.raises(SystemExit):
        args = calculate_rmsd.parse_arguments(
            "--reorder --reorder-method do-not-exists filea fileb".split()
        )
        print(args)


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
