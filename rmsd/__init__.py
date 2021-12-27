from . import calculate_rmsd
from .calculate_rmsd import *
from .calculate_rmsd import __doc__, __version__

__all__ = [
    "brute_permutation",
    "centroid",
    "check_reflections",
    "generate_permutations",
    "get_coordinates",
    "get_coordinates_pdb",
    "get_coordinates_xyz",
    "hungarian",
    "kabsch",
    "kabsch_rmsd",
    "kabsch_rotate",
    "quaternion_rmsd",
    "quaternion_rotate",
    "quaternion_transform",
    "reorder_brute",
    "reorder_distance",
    "reorder_hungarian",
    "reorder_similarity",
    "rmsd",
    "set_coordinates",
]

if __name__ == "__main__":
    calculate_rmsd.main()
