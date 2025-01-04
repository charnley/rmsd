#!/usr/bin/env python3

__doc__ = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation.

Rotation methods:

    - Kabsch (Default)

    - Quaternion

Reorder methods:

    - Brute

        Brute-force enumeration of all atom-atom pair combinations

    - Distance

        Assign atom-atom pair based on the shortest distance

    - Hungarian

        Assign atom-atom pair based on linear-sum assignment of the distance combination

    - Inertia + Hungarian (Default)

        First, align the molecules with inertia moments. Secondly, use Hungarian from above.

For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

__version__ = "1.6.1"

import argparse
import copy
import gzip
import re
import sys
from functools import partial
from pathlib import Path
from typing import Any, Iterator, List, Optional, Protocol, Set, Tuple, Union

import numpy as np
from numpy import ndarray
from scipy.optimize import linear_sum_assignment  # type: ignore
from scipy.spatial import distance_matrix  # type: ignore
from scipy.spatial.distance import cdist  # type: ignore

try:
    import qmllib  # type: ignore
    from qmllib.kernels import laplacian_kernel  # type: ignore
    from qmllib.representations import generate_fchl19  # type: ignore
except ImportError:
    qmllib = None


METHOD_KABSCH = "kabsch"
METHOD_QUATERNION = "quaternion"
METHOD_NOROTATION = "none"
ROTATION_METHODS = [METHOD_KABSCH, METHOD_QUATERNION, METHOD_NOROTATION]

REORDER_NONE = "none"
REORDER_QML = "qml"
REORDER_HUNGARIAN = "hungarian"
REORDER_INERTIA_HUNGARIAN = "inertia-hungarian"
REORDER_BRUTE = "brute"
REORDER_DISTANCE = "distance"
REORDER_METHODS = [
    REORDER_NONE,
    REORDER_QML,
    REORDER_HUNGARIAN,
    REORDER_INERTIA_HUNGARIAN,
    REORDER_BRUTE,
    REORDER_DISTANCE,
]

FORMAT_XYZ = "xyz"
FORMAT_PDB = "pdb"
FORMATS = [FORMAT_XYZ, FORMAT_PDB]

AXIS_SWAPS = np.array([[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 1, 0], [2, 0, 1]])

AXIS_REFLECTIONS = np.array(
    [
        [1, 1, 1],
        [-1, 1, 1],
        [1, -1, 1],
        [1, 1, -1],
        [-1, -1, 1],
        [-1, 1, -1],
        [1, -1, -1],
        [-1, -1, -1],
    ]
)

ELEMENT_WEIGHTS = {
    1: 1.00797,
    2: 4.00260,
    3: 6.941,
    4: 9.01218,
    5: 10.81,
    6: 12.011,
    7: 14.0067,
    8: 15.9994,
    9: 18.998403,
    10: 20.179,
    11: 22.98977,
    12: 24.305,
    13: 26.98154,
    14: 28.0855,
    15: 30.97376,
    16: 32.06,
    17: 35.453,
    19: 39.0983,
    18: 39.948,
    20: 40.08,
    21: 44.9559,
    22: 47.90,
    23: 50.9415,
    24: 51.996,
    25: 54.9380,
    26: 55.847,
    28: 58.70,
    27: 58.9332,
    29: 63.546,
    30: 65.38,
    31: 69.72,
    32: 72.59,
    33: 74.9216,
    34: 78.96,
    35: 79.904,
    36: 83.80,
    37: 85.4678,
    38: 87.62,
    39: 88.9059,
    40: 91.22,
    41: 92.9064,
    42: 95.94,
    43: 98,
    44: 101.07,
    45: 102.9055,
    46: 106.4,
    47: 107.868,
    48: 112.41,
    49: 114.82,
    50: 118.69,
    51: 121.75,
    53: 126.9045,
    52: 127.60,
    54: 131.30,
    55: 132.9054,
    56: 137.33,
    57: 138.9055,
    58: 140.12,
    59: 140.9077,
    60: 144.24,
    61: 145,
    62: 150.4,
    63: 151.96,
    64: 157.25,
    65: 158.9254,
    66: 162.50,
    67: 164.9304,
    68: 167.26,
    69: 168.9342,
    70: 173.04,
    71: 174.967,
    72: 178.49,
    73: 180.9479,
    74: 183.85,
    75: 186.207,
    76: 190.2,
    77: 192.22,
    78: 195.09,
    79: 196.9665,
    80: 200.59,
    81: 204.37,
    82: 207.2,
    83: 208.9804,
    84: 209,
    85: 210,
    86: 222,
    87: 223,
    88: 226.0254,
    89: 227.0278,
    91: 231.0359,
    90: 232.0381,
    93: 237.0482,
    92: 238.029,
    94: 242,
    95: 243,
    97: 247,
    96: 247,
    102: 250,
    98: 251,
    99: 252,
    108: 255,
    109: 256,
    100: 257,
    101: 258,
    103: 260,
    104: 261,
    107: 262,
    105: 262,
    106: 263,
    110: 269,
    111: 272,
    112: 277,
}

ELEMENT_NAMES = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    114: "Uuq",
    116: "Uuh",
}

NAMES_ELEMENT = {value: key for key, value in ELEMENT_NAMES.items()}


class ReorderCallable(Protocol):
    def __call__(
        self,
        p_atoms: ndarray,
        q_atoms: ndarray,
        p_coord: ndarray,
        q_coord: ndarray,
        **kwargs: Any,
    ) -> ndarray:
        """
        Protocol for a reorder callable function

        Return:
            ndarray dtype=int  # Array of indices
        """
        ...  # pragma: no cover


class RmsdCallable(Protocol):
    def __call__(
        self,
        P: ndarray,
        Q: ndarray,
        **kwargs: Any,
    ) -> float:
        """
        Protocol for a rotation callable function

        return:
            RMSD after rotation
        """
        ...  # pragma: no cover


def str_atom(atom: int) -> str:
    """
    Convert atom type from integer to string

    Parameters
    ----------
    atoms : string

    Returns
    -------
    atoms : integer

    """
    return ELEMENT_NAMES[atom]


def int_atom(atom: str) -> int:
    """
    Convert atom type from string to integer

    Parameters
    ----------
    atoms : string

    Returns
    -------
    atoms : integer
    """

    atom = atom.capitalize().strip()
    return NAMES_ELEMENT[atom]


def rmsd(P: ndarray, Q: ndarray, **kwargs) -> float:
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : array
        (N,D) matrix, where N is points and D is dimension.
    W : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
        Root-mean-square deviation between the two vectors
    """
    diff = P - Q
    return np.sqrt((diff * diff).sum() / P.shape[0])


def kabsch_rmsd(
    P: ndarray,
    Q: ndarray,
    W: Optional[ndarray] = None,
    translate: bool = False,
    **kwargs: Any,
) -> float:
    """
    Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.
    An optional vector of weights W may be provided.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    W : array or None
        (N) vector, where N is points.
    translate : bool
        Use centroids to translate vector P and Q unto each other.

    Returns
    -------
    rmsd : float
        root-mean squared deviation
    """

    if translate:
        Q = Q - centroid(Q)
        P = P - centroid(P)

    if W is not None:
        return kabsch_weighted_rmsd(P, Q, W)

    P = kabsch_rotate(P, Q)
    return rmsd(P, Q)


def kabsch_rotate(P: ndarray, Q: ndarray) -> ndarray:
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    P : array
        (N,D) matrix, where N is points and D is dimension,
        rotated

    """
    U = kabsch(P, Q)

    # Rotate P
    P = np.dot(P, U)
    return P


def kabsch_fit(P: ndarray, Q: ndarray, W: Optional[ndarray] = None) -> ndarray:
    """
    Rotate and translate matrix P unto matrix Q using Kabsch algorithm.
    An optional vector of weights W may be provided.

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    W : array or None
        (N) vector, where N is points.

    Returns
    -------
    P : array
        (N,D) matrix, where N is points and D is dimension,
        rotated and translated.

    """
    if W is not None:
        P, _ = kabsch_weighted_fit(P, Q, W, return_rmsd=False)
    else:
        QC = centroid(Q)
        Q = Q - QC
        P = P - centroid(P)
        P = kabsch_rotate(P, Q) + QC
    return P


def kabsch(P: ndarray, Q: ndarray) -> ndarray:
    """
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    The algorithm works in three steps:
    - a centroid translation of P and Q (assumed done before this function
      call)
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    U : matrix
        Rotation matrix (D,D)
    """

    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U
    U: ndarray = np.dot(V, W)

    return U


def kabsch_weighted(
    P: ndarray, Q: ndarray, W: Optional[ndarray] = None
) -> Tuple[ndarray, ndarray, float]:
    """
    Using the Kabsch algorithm with two sets of paired point P and Q.
    Each vector set is represented as an NxD matrix, where D is the
    dimension of the space.
    An optional vector of weights W may be provided.

    Note that this algorithm does not require that P and Q have already
    been overlayed by a centroid translation.

    The function returns the rotation matrix U, translation vector V,
    and RMS deviation between Q and P', where P' is:

        P' = P * U + V

    For more info see http://en.wikipedia.org/wiki/Kabsch_algorithm

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    W : array or None
        (N) vector, where N is points.

    Returns
    -------
    U    : matrix
           Rotation matrix (D,D)
    V    : vector
           Translation vector (D)
    RMSD : float
           Root mean squared deviation between P and Q
    """
    # Computation of the weighted covariance matrix
    CMP = np.zeros(3)
    CMQ = np.zeros(3)
    C = np.zeros((3, 3))
    if W is None:
        W = np.ones(len(P)) / len(P)
    W = np.array([W, W, W]).T
    # NOTE UNUSED psq = 0.0
    # NOTE UNUSED qsq = 0.0
    iw = 3.0 / W.sum()
    n = len(P)
    for i in range(3):
        for j in range(n):
            for k in range(3):
                C[i, k] += P[j, i] * Q[j, k] * W[j, i]
    CMP = (P * W).sum(axis=0)
    CMQ = (Q * W).sum(axis=0)
    PSQ = (P * P * W).sum() - (CMP * CMP).sum() * iw
    QSQ = (Q * Q * W).sum() - (CMQ * CMQ).sum() * iw
    C = (C - np.outer(CMP, CMQ) * iw) * iw

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # Create Rotation matrix U, translation vector V, and calculate RMSD:
    U = np.dot(V, W)
    msd = (PSQ + QSQ) * iw - 2.0 * S.sum()
    if msd < 0.0:
        msd = 0.0
    rmsd_ = np.sqrt(msd)
    V = np.zeros(3)
    for i in range(3):
        t = (U[i, :] * CMQ).sum()
        V[i] = CMP[i] - t
    V = V * iw
    return U, V, rmsd_


def kabsch_weighted_fit(
    P: ndarray,
    Q: ndarray,
    W: Optional[ndarray] = None,
    return_rmsd: bool = False,
) -> Tuple[ndarray, Optional[float]]:
    """
    Fit P to Q with optional weights W.
    Also returns the RMSD of the fit if return_rmsd=True.

    Parameters
    ----------
    P    : array
           (N,D) matrix, where N is points and D is dimension.
    Q    : array
           (N,D) matrix, where N is points and D is dimension.
    W    : vector
           (N) vector, where N is points
    rmsd : Bool
           If True, rmsd is returned as well as the fitted coordinates.

    Returns
    -------
    P'   : array
           (N,D) matrix, where N is points and D is dimension.
    RMSD : float
           if the function is called with rmsd=True
    """
    rmsd_: float
    R, T, rmsd_ = kabsch_weighted(Q, P, W)
    PNEW: ndarray = np.dot(P, R.T) + T
    if return_rmsd:
        return (PNEW, rmsd_)

    return (PNEW, None)


def kabsch_weighted_rmsd(P: ndarray, Q: ndarray, W: Optional[ndarray] = None) -> float:
    """
    Calculate the RMSD between P and Q with optional weights W

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    W : vector
        (N) vector, where N is points

    Returns
    -------
    RMSD : float
    """
    _, _, w_rmsd = kabsch_weighted(P, Q, W)
    return w_rmsd


def quaternion_rmsd(P: ndarray, Q: ndarray, **kwargs: Any) -> float:
    """
    Rotate matrix P unto Q and calculate the RMSD
    based on doi:10.1016/1049-9660(91)90036-O

    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rmsd : float
    """
    rot = quaternion_rotate(P, Q)
    P = np.dot(P, rot)
    return rmsd(P, Q)


def quaternion_transform(r: ndarray) -> ndarray:
    """
    Get optimal rotation
    note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot: ndarray = Wt_r.dot(Q_r)[:3, :3]
    return rot


def makeW(r1: float, r2: float, r3: float, r4: float = 0) -> ndarray:
    """
    matrix involved in quaternion rotation
    """
    W = np.asarray(
        [
            [r4, r3, -r2, r1],
            [-r3, r4, r1, r2],
            [r2, -r1, r4, r3],
            [-r1, -r2, -r3, r4],
        ]
    )
    return W


def makeQ(r1: float, r2: float, r3: float, r4: float = 0) -> ndarray:
    """
    matrix involved in quaternion rotation
    """
    Q = np.asarray(
        [
            [r4, -r3, r2, r1],
            [r3, r4, -r1, r2],
            [-r2, r1, r4, r3],
            [-r1, -r2, -r3, r4],
        ]
    )
    return Q


def quaternion_rotate(X: ndarray, Y: ndarray) -> ndarray:
    """
    Calculate the rotation

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.
    Y: array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    rot : matrix
        Rotation matrix (D,D)
    """
    N = X.shape[0]
    W = np.asarray([makeW(*Y[k]) for k in range(N)])
    Q = np.asarray([makeQ(*X[k]) for k in range(N)])
    Qt_dot_W = np.asarray([np.dot(Q[k].T, W[k]) for k in range(N)])
    # NOTE UNUSED W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    A = np.sum(Qt_dot_W, axis=0)
    eigen = np.linalg.eigh(A)
    r = eigen[1][:, eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def centroid(X: ndarray) -> ndarray:
    """
    Centroid is the mean position of all the points in all of the coordinate
    directions, from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid

    C = sum(X)/len(X)

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : ndarray
        centroid
    """
    C: ndarray = X.mean(axis=0)
    return C


def hungarian_vectors(
    p_vecs: ndarray, q_vecs: ndarray, sigma: float = 1e-0, use_kernel: bool = True
) -> ndarray:
    """

    Hungarian cost assignment of a similiarty molecule kernel.

    Note: Assumes p and q are atoms of same type

    Parameters
    ----------
    p_vecs : array
        (N,L) matrix, where N is no. of atoms and L is representation length
    q_vecs : array
        (N,L) matrix, where N is no. of atoms and L is representation length

    Returns
    -------
    indices_b : array
        (N) view vector of reordered assignment

    """

    if use_kernel:
        # Calculate cost matrix from similarity kernel
        kernel = laplacian_kernel(p_vecs, q_vecs, sigma)
        kernel *= -1.0
        kernel += 1.0

    else:
        kernel = distance_matrix(p_vecs, q_vecs)

    _, indices_q = linear_sum_assignment(kernel)

    return indices_q


def reorder_similarity(
    p_atoms: ndarray,
    q_atoms: ndarray,
    p_coord: ndarray,
    q_coord: ndarray,
    use_kernel: bool = True,
    **kwargs: Any,
) -> ndarray:
    """
    Re-orders the input atom list and xyz coordinates using QML similarity
    the Hungarian method for assignment.

    Parameters
    ----------
    p_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    p_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    p_coord : array
        (N,D) matrix, where N is points and D is dimension
    q_coord : array
        (N,D) matrix, where N is points and D is dimension

    Returns
    -------
    view_reorder : array
             (N,1) matrix, reordered indexes of atom alignment based on the
             coordinates of the atoms
    """

    elements = np.unique(p_atoms)
    n_atoms = p_atoms.shape[0]
    distance_cut = 20.0

    parameters = {
        "elements": elements,
        "pad": n_atoms,
        "rcut": distance_cut,
        "acut": distance_cut,
    }

    p_vecs = generate_fchl19(p_atoms, p_coord, **parameters)
    q_vecs = generate_fchl19(q_atoms, q_coord, **parameters)

    # generate full view from q shape to fill in atom view on the fly
    view_reorder = np.zeros(q_atoms.shape, dtype=int)

    for atom in elements:

        (p_atom_idx,) = np.where(p_atoms == atom)
        (q_atom_idx,) = np.where(q_atoms == atom)

        p_vecs_atom = p_vecs[p_atom_idx]
        q_vecs_atom = q_vecs[q_atom_idx]

        view = hungarian_vectors(p_vecs_atom, q_vecs_atom, use_kernel=use_kernel)
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


def reorder_distance(
    p_atoms: ndarray,
    q_atoms: ndarray,
    p_coord: ndarray,
    q_coord: ndarray,
    **kwargs: Any,
) -> ndarray:
    """
    Re-orders the input atom list and xyz coordinates by atom type and then by
    distance of each atom from the centroid.

    Parameters
    ----------
    atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    coord : array
        (N,D) matrix, where N is points and D is dimension

    Returns
    -------
    atoms_reordered : array
        (N,1) matrix, where N is points holding the ordered atoms' names
    coords_reordered : array
        (N,D) matrix, where N is points and D is dimension (rows re-ordered)
    """

    # Find unique atoms
    unique_atoms = np.unique(p_atoms)

    # generate full view from q shape to fill in atom view on the fly
    view_reorder = np.zeros(q_atoms.shape, dtype=int)

    for atom in unique_atoms:

        (p_atom_idx,) = np.where(p_atoms == atom)
        (q_atom_idx,) = np.where(q_atoms == atom)

        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]

        # Calculate distance from each atom to centroid
        A_norms = np.linalg.norm(A_coord, axis=1)
        B_norms = np.linalg.norm(B_coord, axis=1)

        reorder_indices_A = np.argsort(A_norms)
        reorder_indices_B = np.argsort(B_norms)

        # Project the order of P onto Q
        translator = np.argsort(reorder_indices_A)
        view = reorder_indices_B[translator]
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


def _hungarian(A: ndarray, B: ndarray) -> ndarray:
    """
    Hungarian reordering, minimizing the cost distances[winner_rows, winner_cols].sum())

    Assume A and B are coordinates for atoms of SAME type only.
    """

    distances = cdist(A, B, "euclidean")

    # Perform Hungarian analysis on distance matrix between atoms of 1st
    # structure and trial structure
    winner_rows, winner_cols = linear_sum_assignment(distances)

    return winner_cols


def reorder_hungarian(
    p_atoms: ndarray,
    q_atoms: ndarray,
    p_coord: ndarray,
    q_coord: ndarray,
    **kwargs: Any,
) -> ndarray:
    """
    Re-orders the input atom array and coordinates using the Hungarian-Distance
    sum assignment method. Returns view of q-atoms.

    Parameters
    ----------
    p_atoms : array
        (N,1) matrix, where N is points holding the atom types
    p_atoms : array
        (N,1) matrix, where N is points holding the atom types
    p_coord : array
        (N,D) matrix, where N is points and D is dimension
    q_coord : array
        (N,D) matrix, where N is points and D is dimension

    Returns
    -------
    view_reorder : array
             (N,1) matrix, reordered indexes of atom alignment based on the
             coordinates of the atoms

    """

    # Find unique atoms
    unique_atoms = np.unique(p_atoms)

    # generate full view from q shape to fill in atom view on the fly
    view_reorder = np.zeros(q_atoms.shape, dtype=int)
    view_reorder -= 1

    for atom in unique_atoms:
        (p_atom_idx,) = np.where(p_atoms == atom)
        (q_atom_idx,) = np.where(q_atoms == atom)

        _p_coord = p_coord[p_atom_idx]
        _q_coord = q_coord[q_atom_idx]

        view = _hungarian(_p_coord, _q_coord)
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


def reorder_inertia_hungarian(
    p_atoms: ndarray,
    q_atoms: ndarray,
    p_coord: ndarray,
    q_coord: ndarray,
    **kwargs: Any,
) -> ndarray:
    """
    First, align structures with the intertia moment eignvectors, then using
    distance hungarian, assign the best possible atom pair combinations. While
    also checking all possible reflections of intertia moments, selecting the
    one with minimal RMSD.

    Parameters
    ----------
    p_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    p_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    p_coord : array
        (N,D) matrix, where N is points and D is dimension
    q_coord : array
        (N,D) matrix, where N is points and D is dimension

    Returns
    -------
    view_reorder : array
             (N,1) matrix, reordered indexes of `Q` atom alignment based on the
             coordinates of the atoms

    """

    p_coord = np.array(p_coord, copy=True)
    q_coord = np.array(q_coord, copy=True)

    p_coord -= get_cm(p_atoms, p_coord)
    q_coord -= get_cm(q_atoms, p_coord)

    # Calculate inertia vectors for both structures
    inertia_p = get_inertia_tensor(p_atoms, p_coord)
    eigval_p, eigvec_p = np.linalg.eig(inertia_p)

    eigvec_p = eigvec_p.T
    eigvec_p = eigvec_p[np.argsort(eigval_p)]
    eigvec_p = eigvec_p.T

    inertia_q = get_inertia_tensor(q_atoms, q_coord)
    eigval_q, eigvec_q = np.linalg.eig(inertia_q)

    eigvec_q = eigvec_q.T
    eigvec_q = eigvec_q[np.argsort(eigval_q)]
    eigvec_q = eigvec_q.T

    # Reset the p coords, so the inertia vectors align with axis
    p_coord = np.dot(p_coord, eigvec_p)

    best_rmsd = np.inf
    best_review = np.arange(len(p_atoms))

    for mirror in AXIS_REFLECTIONS:

        tmp_eigvec = eigvec_q * mirror.T
        tmp_coord = np.dot(q_coord, tmp_eigvec)

        test_review = reorder_hungarian(p_atoms, q_atoms, p_coord, tmp_coord)
        test_rmsd = kabsch_rmsd(tmp_coord[test_review], p_coord)

        if test_rmsd < best_rmsd:
            best_rmsd = test_rmsd
            best_review = test_review

    return best_review


def generate_permutations(elements: List[int], n: int) -> Iterator[List[int]]:
    """
    Heap's algorithm for generating all n! permutations in a list
    https://en.wikipedia.org/wiki/Heap%27s_algorithm

    """
    c = [0] * n
    yield elements
    i = 0
    while i < n:
        if c[i] < i:
            if i % 2 == 0:
                elements[0], elements[i] = elements[i], elements[0]
            else:
                elements[c[i]], elements[i] = elements[i], elements[c[i]]
            yield elements
            c[i] += 1
            i = 0
        else:
            c[i] = 0
            i += 1


def brute_permutation(A: ndarray, B: ndarray) -> ndarray:
    """
    Re-orders the input atom list and xyz coordinates using the brute force
    method of permuting all rows of the input coordinates

    Parameters
    ----------
    A : array
        (N,D) matrix, where N is points and D is dimension
    B : array
        (N,D) matrix, where N is points and D is dimension

    Returns
    -------
    view : array
        (N,1) matrix, reordered view of B projected to A
    """

    rmsd_min = np.inf
    view_min: ndarray

    # Sets initial ordering for row indices to [0, 1, 2, ..., len(A)], used in
    # brute-force method

    num_atoms = A.shape[0]
    initial_order = list(range(num_atoms))

    for reorder_indices in generate_permutations(initial_order, num_atoms):

        # Re-order the atom array and coordinate matrix
        coords_ordered = B[reorder_indices]

        # Calculate the RMSD between structure 1 and the Hungarian re-ordered
        # structure 2
        rmsd_temp = kabsch_rmsd(A, coords_ordered)

        # Replaces the atoms and coordinates with the current structure if the
        # RMSD is lower
        if rmsd_temp < rmsd_min:
            rmsd_min = rmsd_temp
            view_min = np.asarray(copy.deepcopy(reorder_indices))

    return view_min


def reorder_brute(
    p_atoms: ndarray,
    q_atoms: ndarray,
    p_coord: ndarray,
    q_coord: ndarray,
    **kwargs: Any,
) -> ndarray:
    """
    Re-orders the input atom list and xyz coordinates using all permutation of
    rows (using optimized column results)

    Parameters
    ----------
    p_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    q_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    p_coord : array
        (N,D) matrix, where N is points and D is dimension
    q_coord : array
        (N,D) matrix, where N is points and D is dimension

    Returns
    -------
    view_reorder : array
        (N,1) matrix, reordered indexes of atom alignment based on the
        coordinates of the atoms

    """

    # Find unique atoms
    unique_atoms = np.unique(p_atoms)

    # generate full view from q shape to fill in atom view on the fly
    view_reorder = np.zeros(q_atoms.shape, dtype=int)
    view_reorder -= 1

    for atom in unique_atoms:
        (p_atom_idx,) = np.where(p_atoms == atom)
        (q_atom_idx,) = np.where(q_atoms == atom)

        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]

        view = brute_permutation(A_coord, B_coord)
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


def check_reflections(
    p_atoms: ndarray,
    q_atoms: ndarray,
    p_coord: ndarray,
    q_coord: ndarray,
    reorder_method: Optional[ReorderCallable] = None,
    rmsd_method: RmsdCallable = kabsch_rmsd,
    keep_stereo: bool = False,
) -> Tuple[float, ndarray, ndarray, ndarray]:
    """
    Minimize RMSD using reflection planes for molecule P and Q

    Warning: This will affect stereo-chemistry

    Parameters
    ----------
    p_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    q_atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    p_coord : array
        (N,D) matrix, where N is points and D is dimension
    q_coord : array
        (N,D) matrix, where N is points and D is dimension

    Returns
    -------
    min_rmsd
    min_swap
    min_reflection
    min_review

    """

    if reorder_method is None:
        assert (p_atoms == q_atoms).all(), "No reorder method selected, but atoms are not ordered"

    min_rmsd = np.inf
    min_swap: ndarray
    min_reflection: ndarray
    min_review: ndarray = np.array(range(len(p_atoms)))
    tmp_review: ndarray = min_review
    swap_mask = [1, -1, -1, 1, -1, 1]
    reflection_mask = [1, -1, -1, -1, 1, 1, 1, -1]

    for swap, i in zip(AXIS_SWAPS, swap_mask):
        for reflection, j in zip(AXIS_REFLECTIONS, reflection_mask):

            # skip enantiomers
            if keep_stereo and i * j == -1:
                continue

            tmp_atoms = copy.copy(q_atoms)
            tmp_coord = copy.deepcopy(q_coord)
            tmp_coord = tmp_coord[:, swap]
            tmp_coord = np.dot(tmp_coord, np.diag(reflection))
            tmp_coord -= centroid(tmp_coord)

            # Reorder
            if reorder_method is not None:
                tmp_review = reorder_method(p_atoms, tmp_atoms, p_coord, tmp_coord)
                tmp_coord = tmp_coord[tmp_review]
                tmp_atoms = tmp_atoms[tmp_review]

            # Rotation
            this_rmsd = rmsd_method(p_coord, tmp_coord)

            if this_rmsd < min_rmsd:
                min_rmsd = this_rmsd
                min_swap = swap
                min_reflection = reflection
                min_review = tmp_review

    assert (p_atoms == q_atoms[min_review]).all(), "error: Not aligned"

    return min_rmsd, min_swap, min_reflection, min_review


def rotation_matrix_vectors(v1: ndarray, v2: ndarray) -> ndarray:
    """
    Returns the rotation matrix that rotates v1 onto v2
    using Rodrigues' rotation formula.
    (see https://math.stackexchange.com/a/476311)
    ----------
    v1 : array
        Dim 3 float array
    v2 : array
        Dim 3 float array

    Return
    ------
    output : 3x3 matrix
        Rotation matrix
    """

    rot: ndarray

    if (v1 == v2).all():
        rot = np.eye(3)

    # return a rotation of pi around the y-axis
    elif (v1 == -v2).all():
        rot = np.array([[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])

    else:
        v = np.cross(v1, v2)
        s = np.linalg.norm(v)
        c = np.vdot(v1, v2)

        vx = np.array([[0.0, -v[2], v[1]], [v[2], 0.0, -v[0]], [-v[1], v[0], 0.0]])

        rot = np.eye(3) + vx + np.dot(vx, vx) * ((1.0 - c) / (s * s))

    return rot


def get_cm(atoms: ndarray, V: ndarray) -> ndarray:
    """
    Get the center of mass of V.
    ----------
    atoms : list
        List of atomic types
    V : array
        (N,3) matrix of atomic coordinates

    Return
    ------
    output : (3) array
        The CM vector
    """

    weights: Union[List[float], ndarray] = [ELEMENT_WEIGHTS[x] for x in atoms]
    weights = np.asarray(weights)
    center_of_mass: ndarray = np.average(V, axis=0, weights=weights)

    return center_of_mass


def get_inertia_tensor(atoms: ndarray, V: ndarray) -> ndarray:
    """
    Get the tensor of intertia of V.
    ----------
    atoms : list
        List of atomic types
    V : array
        (N,3) matrix of atomic coordinates

    Return
    ------
    output : 3x3 float matrix
        The tensor of inertia
    """

    CV = V - get_cm(atoms, V)

    Ixx = 0.0
    Iyy = 0.0
    Izz = 0.0
    Ixy = 0.0
    Ixz = 0.0
    Iyz = 0.0

    for sp, acoord in zip(atoms, CV):
        amass = ELEMENT_WEIGHTS[sp]
        Ixx += amass * (acoord[1] * acoord[1] + acoord[2] * acoord[2])
        Iyy += amass * (acoord[0] * acoord[0] + acoord[2] * acoord[2])
        Izz += amass * (acoord[0] * acoord[0] + acoord[1] * acoord[1])
        Ixy += -amass * acoord[0] * acoord[1]
        Ixz += -amass * acoord[0] * acoord[2]
        Iyz += -amass * acoord[1] * acoord[2]

    coordinates = V
    com = get_cm(atoms, V)

    atomic_masses = np.asarray([ELEMENT_WEIGHTS[a] for a in atoms])
    coordinates -= com

    mass_matrix = np.diag(atomic_masses)
    helper = coordinates.T.dot(mass_matrix).dot(coordinates)
    inertia_tensor: np.ndarray = np.diag(np.ones(3)) * helper.trace() - helper
    return inertia_tensor

    return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])


def get_principal_axis(atoms: ndarray, V: ndarray) -> ndarray:
    """
    Get the molecule's principal axis.
    ----------
    atoms : list
        List of atomic types
    V : array
        (N,3) matrix of atomic coordinates

    Return
    ------
    output : array
        Array of dim 3 containing the principal axis
    """
    inertia = get_inertia_tensor(atoms, V)

    eigval, eigvec = np.linalg.eig(inertia)

    principal_axis: ndarray = eigvec[np.argmax(eigval)]

    return principal_axis


def set_coordinates(
    atoms: ndarray, V: ndarray, title: str = "", decimals: int = 8, set_atoms_as_symbols=True
) -> str:
    """
    Print coordinates V with corresponding atoms to stdout in XYZ format.
    Parameters
    ----------
    atoms : list
        List of atomic types
    V : array
        (N,3) matrix of atomic coordinates
    title : string (optional)
        Title of molecule
    decimals : int (optional)
        number of decimals for the coordinates

    Return
    ------
    output : str
        Molecule in XYZ format

    """
    N, D = V.shape

    if N != len(atoms):
        raise ValueError("Mismatch between expected atoms and coordinate size")

    if not isinstance(atoms[0], str) and set_atoms_as_symbols:
        atoms = np.array([str_atom(atom) for atom in atoms])

    fmt = "{:<2}" + (" {:15." + str(decimals) + "f}") * 3

    out = list()
    out += [str(N)]
    out += [title]

    for i in range(N):
        atom = atoms[i]
        out += [fmt.format(atom, V[i, 0], V[i, 1], V[i, 2])]

    return "\n".join(out)


def get_coordinates(
    filename: Path, fmt: str, is_gzip: bool = False, return_atoms_as_int: bool = False
) -> Tuple[ndarray, ndarray]:
    """
    Get coordinates from filename in format fmt. Supports XYZ and PDB.
    Parameters
    ----------
    filename : string
        Filename to read
    fmt : string
        Format of filename. Either xyz or pdb.
    Returns
    -------
    atoms : list
        List of atomic types
    V : array
        (N,3) where N is number of atoms
    """
    if fmt == "xyz":
        get_func = get_coordinates_xyz

    elif fmt == "pdb":
        get_func = get_coordinates_pdb

    else:
        raise ValueError("Could not recognize file format: {:s}".format(fmt))

    val = get_func(filename, is_gzip=is_gzip, return_atoms_as_int=return_atoms_as_int)

    return val


def _parse_pdb_alphacarbon_line(line: str) -> bool:
    """Try to read Alpha carbons based on PDB column-based format"""

    atom_col = line[12:16]
    atom = atom_col[0:2]
    atom = re.sub(r"\d", " ", atom)
    atom = atom.strip()
    atom = atom.capitalize()
    location = atom_col[2]

    if atom == "C" and location == "A":
        return True

    return False


def _parse_pdb_atom_line(line: str) -> Optional[str]:
    """
    Will try it best to find atom from an atom-line. The standard of PDB
    *should* be column based, however, there are many examples of non-standard
    files. We try our best to have a minimal reader.

    From PDB Format 1992 pdf:

        ATOM Atomic coordinate records for "standard" groups
        HETATM Atomic coordinate records for "non-standard" groups

        Cols. 1 - 4  ATOM
           or 1 - 6  HETATM

              7 - 11 Atom serial number(i)
             13 - 16 Atom name(ii)
             17      Alternate location indicator(iii)
             18 - 20 Residue name(iv,v)
             22      Chain identifier, e.g., A for hemoglobin α chain
             23 - 26 Residue seq. no.
             27      Code for insertions of residues, e.g., 66A, 66B, etc.
             31 - 38 X
             39 - 46 Y Orthogonal Å coordinates
             47 - 54 Z
             55 - 60 Occupancy
             61 - 66 Temperature factor(vi)
             68 - 70 Footnote number

        For (II)

        Within each residue the atoms occur in the order specified by the
        superscripts. The extra oxygen atom of the carboxy terminal amino acid
        is designated OXT

        Four characters are reserved for these atom names. They are assigned as
        follows:

        1-2 Chemical symbol - right justified
        3 Remoteness indicator (alphabetic)
        4 Branch designator (numeric)

        For protein coordinate sets containing hydrogen atoms, the IUPAC-IUB
        rules1 have been followed. Recommendation rule number 4.4 has been
        modified as follows: When more than one hydrogen atom is bonded to a
        single non-hydrogen atom, the hydrogen atom number designation is given
        as the first character of the atom name rather than as the last
        character (e.g. Hβ1 is denoted as 1HB). Exceptions to these rules may
        occur in certain data sets at the depositors’ request. Any such
        exceptions will be delineated clearly in FTNOTE and REMARK records

    but, from [PDB Format Version 2.1]

        In large het groups it sometimes is not possible to follow the
        convention of having the first two characters be the chemical symbol
        and still use atom names that are meaningful to users. A example is
        nicotinamide adenine dinucleotide, atom names begin with an A or N,
        depending on which portion of the molecule they appear in, e.g., AC6 or
        NC6, AN1 or NN1.

        Hydrogen naming sometimes conflicts with IUPAC conventions. For
        example, a hydrogen named HG11 in columns 13 - 16 is differentiated
        from a mercury atom by the element symbol in columns 77 - 78. Columns
        13 - 16 present a unique name for each atom.

    """

    atom_col = line[12:16]
    atom = atom_col[0:2]
    atom = re.sub(r"\d", " ", atom)
    atom = atom.strip()
    atom = atom.capitalize()

    # Highly unlikely that it is Mercury, Helium, Hafnium etc. See comment in
    # function description. [PDB Format v2.1]
    if len(atom) == 2 and atom[0] == "H":
        atom = "H"

    if atom in NAMES_ELEMENT.keys():
        return atom

    tokens = line.split()
    atom = tokens[2][0]
    if atom in NAMES_ELEMENT.keys():
        return atom

    # e.g. 1HD1
    atom = tokens[2][1]
    if atom.upper() == "H":
        return atom

    return None


def _parse_pdb_coord_line(line: str) -> Optional[ndarray]:
    """
    Try my best to coordinates from a PDB ATOM or HETATOM line

    The coordinates should be located in
        31 - 38 X
        39 - 46 Y
        47 - 54 Z

    as defined in PDB, ATOMIC COORDINATE AND BIBLIOGRAPHIC ENTRY FORMAT DESCRIPTION, Feb, 1992
    """

    # If that doesn't work, use hardcoded indices
    try:
        x = line[30:38]
        y = line[38:46]
        z = line[46:54]
        coord = np.asarray([x, y, z], dtype=float)
        return coord

    except ValueError:
        coord = None

    tokens = line.split()

    x_column: Optional[int] = None

    # look for x column
    for i, x in enumerate(tokens):
        if "." in x and "." in tokens[i + 1] and "." in tokens[i + 2]:
            x_column = i
            break

    if x_column is None:
        return None

    # Try to read the coordinates
    try:
        coord = np.asarray(tokens[x_column : x_column + 3], dtype=float)
        return coord
    except ValueError:
        coord = None

    return None


def get_coordinates_pdb(
    filename: Path,
    is_gzip: bool = False,
    return_atoms_as_int: bool = False,
    only_alpha_carbon: bool = False,
) -> Tuple[ndarray, ndarray]:
    """
    Get coordinates from the first chain in a pdb file
    and return a vectorset with all the coordinates.

    Parameters
    ----------
    filename : string
        Filename to read

    Returns
    -------
    atoms : list
        List of atomic types
    V : array
        (N,3) where N is number of atoms
    """

    # PDB files tend to be a bit of a mess. The x, y and z coordinates
    # are supposed to be in column 31-38, 39-46 and 47-54, but this is
    # not always the case.
    # Because of this the three first columns containing a decimal is used.
    # Since the format doesn't require a space between columns, we use the
    # above column indices as a fallback.

    V: Union[List[ndarray], ndarray] = list()
    assert isinstance(V, list)

    # Same with atoms and atom naming.
    # The most robust way to do this is probably
    # to assume that the atomtype is given in column 3.

    atoms: List[str] = list()
    alpha_carbons: List[bool] = list()
    assert isinstance(atoms, list)
    openfunc: Any

    if is_gzip:
        openfunc = gzip.open
        openarg = "rt"
    else:
        openfunc = open
        openarg = "r"

    with openfunc(filename, openarg) as f:
        lines = f.readlines()
        for line in lines:

            if line.startswith("TER") or line.startswith("END"):
                break

            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom = _parse_pdb_atom_line(line)
            if atom is None:
                raise ValueError(f"error: Parsing for atom line: {line}")
            atoms.append(atom)

            coord = _parse_pdb_coord_line(line)
            if coord is None:
                raise ValueError(f"error: Parsing coordinates for line: {line}")
            V.append(coord)

            # Check if alpha-carbon
            is_alpha = _parse_pdb_alphacarbon_line(line)
            alpha_carbons.append(is_alpha)

    if return_atoms_as_int:
        _atoms = np.asarray([int_atom(str(atom)) for atom in atoms])
    else:
        _atoms = np.asarray(atoms)

    V = np.asarray(V)
    assert isinstance(V, ndarray)

    assert isinstance(_atoms, ndarray)
    assert V.shape[0] == _atoms.size

    if only_alpha_carbon:
        # Check that any alpha carbons were found
        if not sum(alpha_carbons):
            raise ValueError("Trying to filter for alpha carbons, but couldn't find any")

        _alpha_carbons = np.asarray(alpha_carbons, dtype=bool)

        V = V[_alpha_carbons, :]
        _atoms = _atoms[_alpha_carbons]

    return _atoms, V


def get_coordinates_xyz_lines(
    lines: List[str], return_atoms_as_int: bool = False
) -> Tuple[ndarray, ndarray]:

    V: Union[List[ndarray], ndarray] = list()
    atoms: Union[List[str], ndarray] = list()
    n_atoms = 0

    assert isinstance(V, list)
    assert isinstance(atoms, list)

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(lines[0])
    except ValueError:
        exit("error: Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    # Use the number of atoms to not read beyond the end of a file
    for lines_read, line in enumerate(lines[2:]):

        line = line.strip()

        if lines_read == n_atoms:
            break

        values = line.split()

        if len(values) < 4:
            atom = re.findall(r"[a-zA-Z]+", line)[0]
            atom = atom.upper()
            numbers = re.findall(r"[-]?\d+\.\d*(?:[Ee][-\+]\d+)?", line)
            numbers = [float(number) for number in numbers]
        else:
            atom = values[0]
            numbers = [float(number) for number in values[1:]]

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) >= 3:
            V.append(np.array(numbers)[:3])
            atoms.append(atom)
        else:
            msg = (
                f"Reading the .xyz file failed in line {lines_read + 2}."
                "Please check the format."
            )
            exit(msg)

    try:
        # I've seen examples where XYZ are written with integer atoms types
        atoms_ = [int(atom) for atom in atoms]
        atoms = [str_atom(atom) for atom in atoms_]

    except ValueError:
        # Correct atom spelling
        atoms = [atom.capitalize() for atom in atoms]

    if return_atoms_as_int:
        atoms_ = [int_atom(atom) for atom in atoms]
        atoms = np.array(atoms_)
    else:
        atoms = np.array(atoms)

    V = np.array(V)

    return atoms, V


def get_coordinates_xyz(
    filename: Path,
    is_gzip: bool = False,
    return_atoms_as_int: bool = False,
) -> Tuple[ndarray, ndarray]:
    """
    Get coordinates from filename and return a vectorset with all the
    coordinates, in XYZ format.

    Parameters
    ----------
    filename : string
        Filename to read

    Returns
    -------
    atoms : list
        List of atomic types
    V : array
        (N,3) where N is number of atoms
    """

    openfunc: Any

    if is_gzip:
        openfunc = gzip.open
        openarg = "rt"
    else:
        openfunc = open
        openarg = "r"

    with openfunc(filename, openarg) as f:
        lines = f.readlines()

    atoms, V = get_coordinates_xyz_lines(lines, return_atoms_as_int=return_atoms_as_int)

    return atoms, V


def parse_arguments(arguments: Optional[Union[str, List[str]]] = None) -> argparse.Namespace:

    description = __doc__

    version_msg = f"""
rmsd {__version__}

See https://github.com/charnley/rmsd for citation information

"""

    epilog = """
"""

    valid_reorder_methods = ", ".join(REORDER_METHODS)
    valid_rotation_methods = ", ".join(ROTATION_METHODS)

    parser = argparse.ArgumentParser(
        usage="calculate_rmsd [options] FILE_A FILE_B",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog,
    )

    # Input structures
    parser.add_argument(
        "structure_a",
        metavar="FILE_A",
        type=str,
        help="structures in .xyz or .pdb format",
    )
    parser.add_argument("structure_b", metavar="FILE_B", type=str)

    # Admin
    parser.add_argument("-v", "--version", action="version", version=version_msg)

    # Rotation
    parser.add_argument(
        "-r",
        "--rotation",
        action="store",
        default="kabsch",
        help=(
            "select rotation method. Valid methods are "
            f"{valid_rotation_methods}. "
            "Default is Kabsch."
        ),
        metavar="METHOD",
        choices=ROTATION_METHODS,
    )

    # Reorder arguments
    parser.add_argument(
        "-e",
        "--reorder",
        action="store_true",
        help="align the atoms of molecules",
    )
    parser.add_argument(
        "--reorder-method",
        action="store",
        default="inertia-hungarian",
        metavar="METHOD",
        help=(
            "select reorder method. Valid method are "
            f"{valid_reorder_methods}. "
            "Default is Inertia-Hungarian."
        ),
        choices=REORDER_METHODS,
    )
    parser.add_argument(
        "-ur",
        "--use-reflections",
        action="store_true",
        help=(
            "scan through reflections in planes "
            "(eg Y transformed to -Y -> X, -Y, Z) "
            "and axis changes, (eg X and Z coords exchanged -> Z, Y, X). "
            "This will affect stereo-chemistry."
        ),
    )
    parser.add_argument(
        "-urks",
        "--use-reflections-keep-stereo",
        action="store_true",
        help=(
            "scan through reflections in planes "
            "(eg Y transformed to -Y -> X, -Y, Z) "
            "and axis changes, (eg X and Z coords exchanged -> Z, Y, X). "
            "Stereo-chemistry will be kept."
        ),
    )

    # Filter
    index_group = parser.add_mutually_exclusive_group()
    index_group.add_argument(
        "--only-alpha-carbons",
        action="store_true",
        help="use only alpha carbons (only for pdb)",
    )
    index_group.add_argument(
        "-nh",
        "--ignore-hydrogen",
        "--no-hydrogen",
        action="store_true",
        help="ignore hydrogens when calculating RMSD",
    )
    index_group.add_argument(
        "--remove-idx",
        nargs="+",
        type=int,
        help="index list of atoms NOT to consider",
        metavar="IDX",
    )
    index_group.add_argument(
        "--add-idx",
        nargs="+",
        type=int,
        help="index list of atoms to consider",
        metavar="IDX",
    )

    parser.add_argument(
        "--format",
        action="store",
        help=f"format of input files. valid format are {', '.join(FORMATS)}.",
        metavar="FMT",
    )
    parser.add_argument(
        "--format-is-gzip",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "-p",
        "--output",
        "--print",
        action="store_true",
        help=(
            "print out structure B, centered and rotated unto structure A's coordinates in XYZ format"
        ),
    )

    parser.add_argument(
        "--print-only-rmsd-atoms",
        action="store_true",
        help=(
            "Print only atoms used in finding optimal RMSD calculation (relevant if filtering e.g. Hydrogens)"
        ),
    )

    args = parser.parse_args(arguments)

    # Check illegal combinations
    if args.output and args.reorder and (args.ignore_hydrogen or args.add_idx or args.remove_idx):
        print(
            "error: Cannot reorder atoms and print structure, when excluding atoms (such as --ignore-hydrogen)"
        )
        sys.exit(5)

    if (
        args.use_reflections
        and args.output
        and (args.ignore_hydrogen or args.add_idx or args.remove_idx)
    ):
        print(
            "error: Cannot use reflections on atoms and print, "
            "when excluding atoms (such as --ignore-hydrogen)"
        )
        sys.exit(5)

    # Check methods
    args.rotation = args.rotation.lower()
    if args.rotation not in ROTATION_METHODS:
        print(
            f"error: Unknown rotation method: '{args.rotation}'. "
            f"Please use {valid_rotation_methods}"
        )
        sys.exit(5)

    # Check reorder methods
    args.reorder_method = args.reorder_method.lower()
    if args.reorder_method not in REORDER_METHODS:
        print(
            f'error: Unknown reorder method: "{args.reorder_method}". '
            f"Please use {valid_reorder_methods}"
        )
        sys.exit(5)

    # Check fileformat
    if args.format is None:
        filename = args.structure_a
        suffixes = Path(filename).suffixes

        if len(suffixes) == 0:
            ext = None

        elif suffixes[-1] == ".gz":
            args.format_is_gzip = True
            ext = suffixes[-2].strip(".")

        else:
            ext = suffixes[-1].strip(".")

        args.format = ext

    # Check if format exist
    if args.format not in FORMATS:
        print(f"error: Format not supported {args.format}")
        sys.exit(5)

    # Check illegal argument
    if args.format != FORMAT_PDB and args.only_alpha_carbons:
        print("Alpha carbons only exist in pdb files")
        sys.exit(5)

    # Check QML is installed
    if args.reorder_method == REORDER_QML and qmllib is None:
        print(
            "'qmllib' is not installed. Package is avaliable from: github.com/qmlcode/qmllib or pip install qmllib."
        )
        sys.exit(1)

    return args


def main(args: Optional[List[str]] = None) -> str:

    # Parse arguments
    settings = parse_arguments(args)

    # Define the read function
    if settings.format == FORMAT_XYZ:
        get_coordinates = partial(
            get_coordinates_xyz, is_gzip=settings.format_is_gzip, return_atoms_as_int=True
        )

    elif settings.format == FORMAT_PDB:
        get_coordinates = partial(
            get_coordinates_pdb,
            is_gzip=settings.format_is_gzip,
            return_atoms_as_int=True,
            only_alpha_carbon=settings.only_alpha_carbons,
        )
    else:
        print(f"Unknown format: {settings.format}")
        sys.exit(1)

    # As default, load the extension as format
    # Parse pdb.gz and xyz.gz as pdb and xyz formats
    p_atoms, p_coord = get_coordinates(
        settings.structure_a,
    )

    q_atoms, q_coord = get_coordinates(
        settings.structure_b,
    )

    p_size = p_coord.shape[0]
    q_size = q_coord.shape[0]

    if not p_size == q_size:
        print("error: Structures not same size")
        sys.exit()

    if np.count_nonzero(p_atoms != q_atoms) and not settings.reorder:
        msg = """
error: Atoms are not in the same order.

Use --reorder to align the atoms (can be expensive for large structures).

Please see --help or documentation for more information or
https://github.com/charnley/rmsd for further examples.
"""
        print(msg)
        sys.exit()

    # Typing
    index: Union[Set[int], List[int], ndarray]

    # Set local view
    p_view: Optional[ndarray] = None
    q_view: Optional[ndarray] = None
    use_view: bool = True

    if settings.ignore_hydrogen:
        (p_view,) = np.where(p_atoms != 1)
        (q_view,) = np.where(q_atoms != 1)

    elif settings.remove_idx:
        index = np.array(list(set(range(p_size)) - set(settings.remove_idx)))
        p_view = index
        q_view = index

    elif settings.add_idx:
        p_view = settings.add_idx
        q_view = settings.add_idx

    else:
        use_view = False

    # Set local view
    if use_view:
        p_coord_sub = copy.deepcopy(p_coord[p_view])
        q_coord_sub = copy.deepcopy(q_coord[q_view])
        p_atoms_sub = copy.deepcopy(p_atoms[p_view])
        q_atoms_sub = copy.deepcopy(q_atoms[q_view])

    else:
        p_coord_sub = copy.deepcopy(p_coord)
        q_coord_sub = copy.deepcopy(q_coord)
        p_atoms_sub = copy.deepcopy(p_atoms)
        q_atoms_sub = copy.deepcopy(q_atoms)

    # Recenter to centroid
    p_cent_sub = centroid(p_coord_sub)
    q_cent_sub = centroid(q_coord_sub)
    p_coord_sub -= p_cent_sub
    q_coord_sub -= q_cent_sub

    rmsd_method: RmsdCallable
    reorder_method: Optional[ReorderCallable]

    # set rotation method
    if settings.rotation == METHOD_KABSCH:
        rmsd_method = kabsch_rmsd
    elif settings.rotation == METHOD_QUATERNION:
        rmsd_method = quaternion_rmsd
    else:
        rmsd_method = rmsd

    # set reorder method
    reorder_method = None
    if settings.reorder_method == REORDER_QML:
        reorder_method = reorder_similarity
    elif settings.reorder_method == REORDER_HUNGARIAN:
        reorder_method = reorder_hungarian
    elif settings.reorder_method == REORDER_INERTIA_HUNGARIAN:
        reorder_method = reorder_inertia_hungarian
    elif settings.reorder_method == REORDER_BRUTE:
        reorder_method = reorder_brute  # pragma: no cover
    elif settings.reorder_method == REORDER_DISTANCE:
        reorder_method = reorder_distance

    # Save the resulting RMSD
    result_rmsd: Optional[float] = None

    # Collect changes to be done on q coords
    q_swap = None
    q_reflection = None
    q_review = None

    if settings.use_reflections:

        result_rmsd, q_swap, q_reflection, q_review = check_reflections(
            p_atoms_sub,
            q_atoms_sub,
            p_coord_sub,
            q_coord_sub,
            reorder_method=reorder_method,
            rmsd_method=rmsd_method,
        )

    elif settings.use_reflections_keep_stereo:

        result_rmsd, q_swap, q_reflection, q_review = check_reflections(
            p_atoms_sub,
            q_atoms_sub,
            p_coord_sub,
            q_coord_sub,
            reorder_method=reorder_method,
            rmsd_method=rmsd_method,
            keep_stereo=True,
        )

    elif settings.reorder:

        assert reorder_method is not None, "Cannot reorder without selecting --reorder method"
        q_review = reorder_method(p_atoms_sub, q_atoms_sub, p_coord_sub, q_coord_sub)

    # If there is a reorder, then apply before print
    if q_review is not None:

        q_atoms_sub = q_atoms_sub[q_review]
        q_coord_sub = q_coord_sub[q_review]

        assert all(
            p_atoms_sub == q_atoms_sub
        ), "error: Structure not aligned. Please submit bug report at http://github.com/charnley/rmsd"

    # Calculate the RMSD value
    if result_rmsd is None:
        result_rmsd = rmsd_method(p_coord_sub, q_coord_sub)

    # print result
    if settings.output:

        if q_swap is not None:
            q_coord_sub = q_coord_sub[:, q_swap]

        if q_reflection is not None:
            q_coord_sub = np.dot(q_coord_sub, np.diag(q_reflection))

        U = kabsch(q_coord_sub, p_coord_sub)

        if settings.print_only_rmsd_atoms or not use_view:
            q_coord_sub = np.dot(q_coord_sub, U)
            q_coord_sub += p_cent_sub
            return set_coordinates(
                q_atoms_sub,
                q_coord_sub,
                title=f"Rotated '{settings.structure_b}' to match '{settings.structure_a}', with a RMSD of {result_rmsd:.8f}",
            )

        # Swap, reflect, rotate and re-center on the full atom and coordinate set
        q_coord -= q_cent_sub

        if q_swap is not None:
            q_coord = q_coord[:, q_swap]

        if q_reflection is not None:
            q_coord = np.dot(q_coord, np.diag(q_reflection))

        q_coord = np.dot(q_coord, U)
        q_coord += p_cent_sub
        return set_coordinates(
            q_atoms,
            q_coord,
            title=f"Rotated {settings.structure_b} to match {settings.structure_a}, with RMSD of {result_rmsd:.8f}",
        )

    return str(result_rmsd)


if __name__ == "__main__":
    result = main()
    print(result)
