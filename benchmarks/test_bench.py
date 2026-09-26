"""Performance benchmarks for rmsd (100-1k atom use case).

Run with::

    pip install -e ".[bench]"
    pytest benchmarks --benchmark-only --benchmark-sort=mean

Each test benchmarks one hot operation in isolation. Datasets:

- ``ethane`` (8 atoms), ``chembl`` (60 atoms) and ``ci2`` (1064 atoms)
  are real files from ``tests/resources``.
- ``syn200`` / ``syn1000`` are seeded synthetic molecules with a
  realistic H/C/N/O mix, used so scaling in the 100-1k window is
  covered without depending on a single PDB.

Nothing here affects the library or the test suite; it is only
collected when ``pytest benchmarks`` is invoked explicitly.
"""

import numpy as np
import pytest

pytest.importorskip("pytest_benchmark", reason="benchmarks need the rmsd[bench] extra")

from rmsd import calculate_rmsd as R

RESOURCE = "tests/resources"

# ---------------------------------------------------------------------------
# Dataset helpers (module-scoped so file I/O is not part of the timing)
# ---------------------------------------------------------------------------


def _synthetic(n_atoms: int, seed: int = 7):
    """Seeded random molecule with realistic element mix + shuffled partner."""
    rng = np.random.default_rng(seed)
    choices = np.array([1, 1, 1, 1, 1, 6, 6, 6, 7, 8], dtype=int)
    atoms = rng.choice(choices, size=n_atoms)
    p_coord = rng.uniform(-10.0, 10.0, size=(n_atoms, 3))
    # Rigid rotation + translation + small noise, then shuffle atom order
    theta = np.deg2rad(37.0)
    rot = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0.0],
            [np.sin(theta), np.cos(theta), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    q_coord = p_coord @ rot + np.array([2.0, -1.0, 0.5])
    q_coord += rng.normal(0.0, 0.01, size=q_coord.shape)
    perm = rng.permutation(n_atoms)
    # Shuffle only within same element so reorder can recover it
    q_atoms = atoms.copy()
    q_shuffled = q_coord.copy()
    for el in np.unique(atoms):
        src = np.where(atoms == el)[0]
        dst = np.where(atoms == el)[0]
        q_shuffled[src] = q_coord[dst[rng.permutation(len(dst))]]
    _ = perm  # keep permutation intent explicit; per-element shuffle above
    return atoms, p_coord, q_atoms, q_shuffled


def _load_xyz(name):
    return R.get_coordinates_xyz(f"{RESOURCE}/{name}", return_atoms_as_int=True)


def _load_pdb(name):
    return R.get_coordinates_pdb(f"{RESOURCE}/{name}", return_atoms_as_int=True)


@pytest.fixture(scope="module")
def ds_ethane():
    p_atoms, p_coord = _load_xyz("ethane.xyz")
    q_atoms, q_coord = _load_xyz("ethane_translate.xyz")
    return p_atoms, p_coord, q_atoms, q_coord


@pytest.fixture(scope="module")
def ds_chembl():
    p_atoms, p_coord = _load_xyz("CHEMBL3039407.xyz")
    q_atoms, q_coord = _load_xyz("CHEMBL3039407_order.xyz")
    return p_atoms, p_coord, q_atoms, q_coord


@pytest.fixture(scope="module")
def ds_ci2():
    # 1064 atoms -- top of the 100-1k window
    p_atoms, p_coord = _load_pdb("ci2_1.pdb")
    q_atoms, q_coord = _load_pdb("ci2_2.pdb")
    return p_atoms, p_coord, q_atoms, q_coord


@pytest.fixture(scope="module")
def ds_syn200():
    return _synthetic(200)


@pytest.fixture(scope="module")
def ds_syn1000():
    return _synthetic(1000)


def _centered(ds):
    p_atoms, p_coord, q_atoms, q_coord = ds
    pc = p_coord - R.centroid(p_coord)
    qc = q_coord - R.centroid(q_coord)
    return p_atoms, pc, q_atoms, qc


# ---------------------------------------------------------------------------
# Rotation / RMSD core (cheap -- all N covered)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_rmsd(benchmark, ds, request):
    p_atoms, p_coord, q_atoms, q_coord = _centered(request.getfixturevalue(ds))
    benchmark(R.rmsd, p_coord, q_coord)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_kabsch(benchmark, ds, request):
    p_atoms, p_coord, q_atoms, q_coord = _centered(request.getfixturevalue(ds))
    benchmark(R.kabsch, p_coord, q_coord)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_kabsch_rmsd(benchmark, ds, request):
    p_atoms, p_coord, q_atoms, q_coord = _centered(request.getfixturevalue(ds))
    benchmark(R.kabsch_rmsd, p_coord, q_coord)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_kabsch_weighted(benchmark, ds, request):
    p_atoms, p_coord, q_atoms, q_coord = _centered(request.getfixturevalue(ds))
    w = np.ones(len(p_atoms))
    benchmark(R.kabsch_weighted, p_coord, q_coord, w)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_quaternion_rmsd(benchmark, ds, request):
    p_atoms, p_coord, q_atoms, q_coord = _centered(request.getfixturevalue(ds))
    benchmark(R.quaternion_rmsd, p_coord, q_coord)


# ---------------------------------------------------------------------------
# Mass / inertia helpers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_get_cm(benchmark, ds, request):
    p_atoms, p_coord, _, _ = request.getfixturevalue(ds)
    benchmark(R.get_cm, p_atoms, p_coord)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_get_inertia_tensor(benchmark, ds, request):
    p_atoms, p_coord, _, _ = request.getfixturevalue(ds)
    benchmark(R.get_inertia_tensor, p_atoms, p_coord)


# ---------------------------------------------------------------------------
# Reordering (Hungarian family -- the dominant cost at 100-1k atoms)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000", "ds_ci2"])
def test_reorder_hungarian(benchmark, ds, request):
    p_atoms, p_coord, q_atoms, q_coord = _centered(request.getfixturevalue(ds))
    benchmark(R.reorder_hungarian, p_atoms, q_atoms, p_coord, q_coord)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200"])
def test_reorder_inertia_hungarian(benchmark, ds, request):
    # 8x Hungarian + SVD; skip syn1000/ci2 here (covered separately below
    # with fewer rounds to keep suite runtime sane)
    p_atoms, p_coord, q_atoms, q_coord = request.getfixturevalue(ds)
    benchmark(R.reorder_inertia_hungarian, p_atoms, q_atoms, p_coord, q_coord)


@pytest.mark.benchmark(min_rounds=3)
def test_reorder_inertia_hungarian_ci2(benchmark, ds_ci2):
    p_atoms, p_coord, q_atoms, q_coord = ds_ci2
    benchmark(R.reorder_inertia_hungarian, p_atoms, q_atoms, p_coord, q_coord)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl", "ds_syn200", "ds_syn1000"])
def test_reorder_distance(benchmark, ds, request):
    p_atoms, p_coord, q_atoms, q_coord = request.getfixturevalue(ds)
    benchmark(R.reorder_distance, p_atoms, q_atoms, p_coord, q_coord)


@pytest.mark.parametrize("ds", ["ds_ethane", "ds_chembl"])
def test_check_reflections_hungarian(benchmark, ds, request):
    # 48x (reorder+rmsd); only small N by default
    p_atoms, p_coord, q_atoms, q_coord = _centered(request.getfixturevalue(ds))
    benchmark(
        R.check_reflections,
        p_atoms,
        q_atoms,
        p_coord,
        q_coord,
        reorder_method=R.reorder_hungarian,
        rmsd_method=R.kabsch_rmsd,
    )
