import pytest
from conftest import RESOURCE_PATH  # type: ignore

import rmsd as rmsdlib


def test_pdb_only_carbon_possible() -> None:

    filename_a = RESOURCE_PATH / "issue98" / "test1.pdb"
    filename_b = RESOURCE_PATH / "issue98" / "test2.pdb"

    # Structure does not have the same size, so standard cannot work
    with pytest.raises(SystemExit):
        cmd = f"{filename_a} {filename_b}"
        _ = rmsdlib.main(cmd.split())

    cmd = f"--reorder --only-alpha-carbons {filename_a} {filename_b}"
    out = rmsdlib.main(cmd.split())

    print(out)
    rmsd = float(out)
    assert isinstance(rmsd, float)
