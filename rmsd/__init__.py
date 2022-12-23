# flake8: noqa
from .calculate_rmsd import *
from .calculate_rmsd import __doc__, __version__

__all__ = [
    "str_atom",
    "int_atom",
    "rmsd",
    "kabsch_rmsd",
    "kabsch_rotate",
    "kabsch_fit",
    "kabsch",
    "kabsch_weighted",
    "kabsch_weighted_fit",
    "kabsch_weighted_rmsd",
    "quaternion_rmsd",
    "quaternion_transform",
    "makeW",
    "makeQ",
    "quaternion_rotate",
    "centroid",
    "hungarian_vectors",
    "reorder_similarity",
    "reorder_distance",
    "hungarian",
    "reorder_hungarian",
    "reorder_inertia_hungarian",
    "generate_permutations",
    "brute_permutation",
    "reorder_brute",
    "check_reflections",
    "rotation_matrix_vectors",
    "get_cm",
    "get_inertia_tensor",
    "get_principal_axis",
    "set_coordinates",
    "get_coordinates",
    "get_coordinates_pdb",
    "get_coordinates_xyz_lines",
    "get_coordinates_xyz",
    "main",
]

if __name__ == "__main__":
    main()  # pragma: no cover
