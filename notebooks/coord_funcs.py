import numpy as np
from numpy import ndarray


def get_rotation_matrix(degrees: float) -> ndarray:
    """https://en.wikipedia.org/wiki/Rotation_matrix"""

    radians = degrees * np.pi / 180.0

    r11 = np.cos(radians)
    r12 = -np.sin(radians)
    r21 = np.sin(radians)
    r22 = np.cos(radians)

    R = np.array([[r11, r12], [r21, r22]])

    return R


def degree2radiant(degrees):
    return degrees * np.pi / 180.0


def rotate_coord(angle, coord, axis=[0, 1]):
    U = get_rotation_matrix(angle)
    _xy = np.dot(coord[:, axis], U)
    _coord = np.array(coord, copy=True)
    _coord[:, :2] = _xy
    return _coord
