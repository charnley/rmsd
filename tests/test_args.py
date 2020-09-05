
import pytest
import rmsd


def test_legal_arguments():

    args = "--rotation kabsch --ignore-hydrogen FILE_A FILE_B"
    args = args.split()

    args = rmsd.parse_arguments(args)

    assert args.reorder is False
    assert args.ignore_hydrogen is True
    assert args.rotation == 'kabsch'


def test_illegal_arguments():

    args = [
        "--rotation kabsch",
        "--reorder",
        "--print",
        "--ignore-hydrogen",
        "FILE_A",
        "FILE_B"
    ]

    with pytest.raises(SystemExit) as exception:
        args = rmsd.parse_arguments(args)

    assert exception.type == SystemExit


def test_illegal_reflection():

    args = [
        "--rotation kabsch",
        "--use-reflections",
        "--print",
        "--ignore-hydrogen",
        "FILE_A",
        "FILE_B"
    ]

    with pytest.raises(SystemExit) as exception:
        args = rmsd.parse_arguments(args)

    assert exception.type == SystemExit


def test_illegal_rotation_method():

    args = [
        "--rotation SuperKabsch",
        "FILE_A",
        "FILE_B"
    ]

    with pytest.raises(SystemExit) as exception:
        args = rmsd.parse_arguments(args)

    assert exception.type == SystemExit


def test_illegal_reorder_method():

    args = [
        "--reorder-method NotImplementedYet",
        "FILE_A",
        "FILE_B"
    ]

    with pytest.raises(SystemExit) as exception:
        args = rmsd.parse_arguments(args)

    assert exception.type == SystemExit
