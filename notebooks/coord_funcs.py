import numpy as np
from numpy import ndarray


def rotation_matrix(sigma: float) -> ndarray:
    """https://en.wikipedia.org/wiki/Rotation_matrix"""

    radians = sigma * np.pi / 180.0

    r11 = np.cos(radians)
    r12 = -np.sin(radians)
    r21 = np.sin(radians)
    r22 = np.cos(radians)

    R = np.array([[r11, r12], [r21, r22]])

    return R


def degree2radiant(degrees):
    return degrees * np.pi / 180.0


# def find_best_inertia(
#     atoms_a,
#     coord_a,
#     atoms_b,
#     coord_b,
#     reorder_method=rmsdlib.hungarian,
#     rotation_method=rmsdlib.kabsch,
# ):
#     """Expect coord are centered"""

#     inertia_a = rmsdlib.get_inertia_tensor(atoms_a, coord_a)
#     eigval_a, eigvec_a = np.linalg.eig(inertia_a)
#     inertia_b = rmsdlib.get_inertia_tensor(atoms_b, coord_b)
#     eigval_b, eigvec_b = np.linalg.eig(inertia_b)

#     # Align by eignval

#     # Enumerate permutation
#     U = rmsdlib.kabsch(eigvec_a, eigvec_b)

#     # U1 = rotation_matrix_vectors(p_axis, q_axis)
#     # U2 = rotation_matrix_vectors(p_axis, -q_axis)

#     return
