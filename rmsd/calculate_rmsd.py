#!/usr/bin/env python
__doc__ = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation.

For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

__version__ = "1.4"

import argparse
import copy
import gzip
import pathlib
import re
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist

try:
    import qml
except ImportError:
    qml = None


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


def str_atom(atom):
    """
    Convert atom type from integer to string

    Parameters
    ----------
    atoms : string

    Returns
    -------
    atoms : integer

    """
    atom = ELEMENT_NAMES[atom]
    return atom


def int_atom(atom):
    """
    Convert atom type from string to integer

    Parameters
    ----------
    atoms : string

    Returns
    -------
    atoms : integer

    """

    atom = atom.capitalize()
    return NAMES_ELEMENT[atom]


def rmsd(V, W):
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
    diff = np.array(V) - np.array(W)
    N = len(V)
    return np.sqrt((diff * diff).sum() / N)


def kabsch_rmsd(P, Q, W=None, translate=False):
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


def kabsch_rotate(P, Q):
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


def kabsch_fit(P, Q, W=None):
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
        P = kabsch_weighted_fit(P, Q, W, return_rmsd=False)
    else:
        QC = centroid(Q)
        Q = Q - QC
        P = P - centroid(P)
        P = kabsch_rotate(P, Q) + QC
    return P


def kabsch(P, Q):
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
    U = np.dot(V, W)

    return U


def kabsch_weighted(P, Q, W=None):
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


def kabsch_weighted_fit(P, Q, W=None, return_rmsd=False):
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
    R, T, rmsd_ = kabsch_weighted(Q, P, W)
    PNEW = np.dot(P, R.T) + T
    if return_rmsd:
        return PNEW, rmsd_
    else:
        return PNEW


def kabsch_weighted_rmsd(P, Q, W=None):
    """
    Calculate the RMSD between P and Q with optional weighhts W

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
    R, T, w_rmsd = kabsch_weighted(P, Q, W)
    return w_rmsd


def quaternion_rmsd(P, Q):
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


def quaternion_transform(r):
    """
    Get optimal rotation
    note: translation will be zero when the centroids of each molecule are the
    same
    """
    Wt_r = makeW(*r).T
    Q_r = makeQ(*r)
    rot = Wt_r.dot(Q_r)[:3, :3]
    return rot


def makeW(r1, r2, r3, r4=0):
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


def makeQ(r1, r2, r3, r4=0):
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


def quaternion_rotate(X, Y):
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


def centroid(X):
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
    C : float
        centroid
    """
    C = X.mean(axis=0)
    return C


def hungarian_vectors(p_vecs, q_vecs, sigma=1e-0, use_kernel=True):
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
        K = qml.kernels.laplacian_kernel(p_vecs, q_vecs, sigma)
        K *= -1.0
        K += 1.0

    else:
        K = distance_matrix(p_vecs, q_vecs)

    # Perform Hungarian analysis on distance matrix between atoms of 1st
    # structure and trial structure
    indices_a, indices_b = linear_sum_assignment(K)

    return indices_b


def reorder_similarity(p_atoms, q_atoms, p_coord, q_coord, use_kernel=True):
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

    if qml is None:
        raise ImportError(
            "QML is not installed. Package is avaliable from"
            "\n github.com/qmlcode/qml"
            "\n pip install qml"
        )

    if isinstance(p_atoms[0], str):
        p_atoms = [int_atom(atom) for atom in p_atoms]
        q_atoms = [int_atom(atom) for atom in q_atoms]
        p_atoms = np.array(p_atoms)
        q_atoms = np.array(q_atoms)

    elements = np.unique(p_atoms)
    n_atoms = p_atoms.shape[0]
    distance_cut = 20.0

    parameters = {
        "elements": elements,
        "pad": n_atoms,
        "rcut": distance_cut,
        "acut": distance_cut,
    }

    p_vecs = qml.representations.generate_fchl_acsf(p_atoms, p_coord, **parameters)

    q_vecs = qml.representations.generate_fchl_acsf(q_atoms, q_coord, **parameters)

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


def reorder_distance(p_atoms, q_atoms, p_coord, q_coord):
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


def hungarian(A, B):
    """
    Hungarian reordering.

    Assume A and B are coordinates for atoms of SAME type only
    """

    # should be kabasch here i think
    distances = cdist(A, B, "euclidean")

    # Perform Hungarian analysis on distance matrix between atoms of 1st
    # structure and trial structure
    indices_a, indices_b = linear_sum_assignment(distances)

    return indices_b


def reorder_hungarian(p_atoms, q_atoms, p_coord, q_coord):
    """
    Re-orders the input atom list and xyz coordinates using the Hungarian
    method (using optimized column results)

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

        view = hungarian(A_coord, B_coord)
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


def reorder_inertia_hungarian(p_atoms, q_atoms, p_coord, q_coord):
    """
    Align the principal intertia axis and then re-orders the input atom list
    and xyz coordinates using the Hungarian method (using optimized column
    results)

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
    # get the principal axis of P and Q
    p_axis = get_principal_axis(p_atoms, p_coord)
    q_axis = get_principal_axis(q_atoms, q_coord)

    # rotate Q onto P considering that the axis are parallel and antiparallel
    U1 = rotation_matrix_vectors(p_axis, q_axis)
    U2 = rotation_matrix_vectors(p_axis, -q_axis)
    q_coord1 = np.dot(q_coord, U1)
    q_coord2 = np.dot(q_coord, U2)

    q_review1 = reorder_hungarian(p_atoms, q_atoms, p_coord, q_coord1)
    q_review2 = reorder_hungarian(p_atoms, q_atoms, p_coord, q_coord2)
    q_coord1 = q_coord1[q_review1]
    q_coord2 = q_coord2[q_review2]

    rmsd1 = kabsch_rmsd(p_coord, q_coord1)
    rmsd2 = kabsch_rmsd(p_coord, q_coord2)

    if rmsd1 < rmsd2:
        return q_review1
    else:
        return q_review2


def generate_permutations(elements, n):
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


def brute_permutation(A, B):
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
    view_min = None

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
            view_min = copy.deepcopy(reorder_indices)

    return view_min


def reorder_brute(p_atoms, q_atoms, p_coord, q_coord):
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
    p_atoms,
    q_atoms,
    p_coord,
    q_coord,
    reorder_method=reorder_hungarian,
    rotation_method=kabsch_rmsd,
    keep_stereo=False,
):
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

    min_rmsd = np.inf
    min_swap = None
    min_reflection = None
    min_review = None
    tmp_review = None
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
            if rotation_method is None:
                this_rmsd = rmsd(p_coord, tmp_coord)
            else:
                this_rmsd = rotation_method(p_coord, tmp_coord)

            if this_rmsd < min_rmsd:
                min_rmsd = this_rmsd
                min_swap = swap
                min_reflection = reflection
                min_review = tmp_review

    if not (p_atoms == q_atoms[min_review]).all():
        print("error: Not aligned")
        quit()

    return min_rmsd, min_swap, min_reflection, min_review


def rotation_matrix_vectors(v1, v2):
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


def get_cm(atoms, V):
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

    if isinstance(atoms[0], str):
        atoms = [int_atom(atom) for atom in atoms]

    weights = [ELEMENT_WEIGHTS[x] for x in atoms]

    center_of_mass = np.average(V, axis=0, weights=weights)

    return center_of_mass


def get_inertia_tensor(atoms, V):
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

    if isinstance(atoms[0], str):
        atoms = [int_atom(atom) for atom in atoms]

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

    return np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])


def get_principal_axis(atoms, V):
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

    return eigvec[np.argmax(eigval)]


def set_coordinates(atoms, V, title="", decimals=8):
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

    if not isinstance(atoms[0], str):
        atoms = [str_atom(atom) for atom in atoms]

    fmt = "{:<2}" + (" {:15." + str(decimals) + "f}") * 3

    out = list()
    out += [str(N)]
    out += [title]

    for i in range(N):
        atom = atoms[i]
        out += [fmt.format(atom, V[i, 0], V[i, 1], V[i, 2])]

    return "\n".join(out)


def print_coordinates(atoms, V, title=""):
    """
    Print coordinates V with corresponding atoms to stdout in XYZ format.

    Parameters
    ----------
    atoms : list
        List of element types
    V : array
        (N,3) matrix of atomic coordinates
    title : string (optional)
        Title of molecule

    """
    print(set_coordinates(atoms, V, title=title))


def get_coordinates(filename, fmt, is_gzip=False, return_atoms_as_int=False):
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
        exit("Could not recognize file format: {:s}".format(fmt))

    val = get_func(filename, is_gzip=is_gzip, return_atoms_as_int=return_atoms_as_int)

    return val


def get_coordinates_pdb(filename, is_gzip=False, return_atoms_as_int=False):
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

    x_column = None
    V = list()

    # Same with atoms and atom naming.
    # The most robust way to do this is probably
    # to assume that the atomtype is given in column 3.

    atoms = list()

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
            if line.startswith("ATOM") or line.startswith("HETATM"):
                tokens = line.split()
                # Try to get the atomtype
                try:
                    atom = tokens[2][0]
                    if atom in ("H", "C", "N", "O", "S", "P"):
                        atoms.append(atom)
                    else:
                        # e.g. 1HD1
                        atom = tokens[2][1]
                        if atom == "H":
                            atoms.append(atom)
                        else:
                            raise Exception

                except ValueError:
                    msg = f"error: Parsing atomtype for the following line:" f" \n{line}"
                    exit(msg)

                if x_column is None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i + 1] and "." in tokens[i + 2]:
                                x_column = i
                                break

                    except IndexError:
                        msg = "error: Parsing coordinates " "for the following line:" f"\n{line}"
                        exit(msg)

                # Try to read the coordinates
                try:
                    V.append(np.asarray(tokens[x_column : x_column + 3], dtype=float))

                except ValueError:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        V.append(np.asarray([x, y, z], dtype=float))
                    except ValueError:
                        msg = "error: Parsing input " "for the following line:" f"\n{line}"
                        exit(msg)

    if return_atoms_as_int:
        atoms = [int_atom(atom) for atom in atoms]

    V = np.asarray(V)
    atoms = np.asarray(atoms)

    assert V.shape[0] == atoms.size

    return atoms, V


def get_coordinates_xyz(filename, is_gzip=False, return_atoms_as_int=False):
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

    if is_gzip:
        openfunc = gzip.open
        openarg = "rt"
    else:
        openfunc = open
        openarg = "r"

    f = openfunc(filename, openarg)
    V = list()
    atoms = list()
    n_atoms = 0

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(f.readline())
    except ValueError:
        exit("error: Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    f.readline()

    # Use the number of atoms to not read beyond the end of a file
    for lines_read, line in enumerate(f):

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

    f.close()

    try:
        # I've seen examples where XYZ are written with integer atoms types
        atoms = [int(atom) for atom in atoms]
        atoms = [str_atom(atom) for atom in atoms]

    except ValueError:
        # Correct atom spelling
        atoms = [atom.capitalize() for atom in atoms]

    if return_atoms_as_int:
        atoms = [int_atom(atom) for atom in atoms]

    atoms = np.array(atoms)
    V = np.array(V)
    return atoms, V


def parse_arguments(args=None):

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
        default="hungarian",
        metavar="METHOD",
        help=(
            "select reorder method. Valid method are "
            f"{valid_reorder_methods}. "
            "Default is Hungarian."
        ),
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
        help="format of input files. valid format are xyz and pdb",
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
            "print out structure B, "
            "centered and rotated unto structure A's coordinates "
            "in XYZ format"
        ),
    )

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    # Check illegal combinations
    if args.output and args.reorder and (args.ignore_hydrogen or args.add_idx or args.remove_idx):
        print(
            "error: Cannot reorder atoms and print structure, "
            "when excluding atoms (such as --ignore-hydrogen)"
        )
        sys.exit()

    if (
        args.use_reflections
        and args.output
        and (args.ignore_hydrogen or args.add_idx or args.remove_idx)
    ):
        print(
            "error: Cannot use reflections on atoms and print, "
            "when excluding atoms (such as --ignore-hydrogen)"
        )
        sys.exit()

    # Check methods
    args.rotation = args.rotation.lower()
    if args.rotation not in ROTATION_METHODS:
        print(
            f"error: Unknown rotation method: '{args.rotation}'. "
            f"Please use {valid_rotation_methods}"
        )
        sys.exit()

    # Check reorder methods
    args.reorder_method = args.reorder_method.lower()
    if args.reorder_method not in REORDER_METHODS:
        print(
            f'error: Unknown reorder method: "{args.reorder_method}". '
            f"Please use {valid_reorder_methods}"
        )
        sys.exit()

    # Check fileformat
    if args.format is None:
        filename = args.structure_a
        suffixes = pathlib.Path(filename).suffixes

        if len(suffixes) == 0:
            ext = None

        elif suffixes[-1] == ".gz":
            args.format_is_gzip = True
            ext = suffixes[-2].strip(".")

        else:
            ext = suffixes[-1].strip(".")

        args.format = ext

    return args


def main(args=None):

    # Parse arguments
    args = parse_arguments(args)

    # As default, load the extension as format
    # Parse pdb.gz and xyz.gz as pdb and xyz formats
    p_all_atoms, p_all = get_coordinates(
        args.structure_a,
        args.format,
        is_gzip=args.format_is_gzip,
        return_atoms_as_int=True,
    )

    q_all_atoms, q_all = get_coordinates(
        args.structure_b,
        args.format,
        is_gzip=args.format_is_gzip,
        return_atoms_as_int=True,
    )

    p_size = p_all.shape[0]
    q_size = q_all.shape[0]

    if not p_size == q_size:
        print("error: Structures not same size")
        sys.exit()

    if np.count_nonzero(p_all_atoms != q_all_atoms) and not args.reorder:
        msg = """
error: Atoms are not in the same order.

Use --reorder to align the atoms (can be expensive for large structures).

Please see --help or documentation for more information or
https://github.com/charnley/rmsd for further examples.
"""
        print(msg)
        sys.exit()

    # Set local view
    p_view = None
    q_view = None

    if args.ignore_hydrogen:
        assert type(p_all_atoms[0]) != str
        assert type(q_all_atoms[0]) != str
        p_view = np.where(p_all_atoms != 1)
        q_view = np.where(q_all_atoms != 1)

    elif args.remove_idx:
        index = range(p_size)
        index = set(index) - set(args.remove_idx)
        index = list(index)
        p_view = index
        q_view = index

    elif args.add_idx:
        p_view = args.add_idx
        q_view = args.add_idx

    # Set local view
    if p_view is None:
        p_coord = copy.deepcopy(p_all)
        q_coord = copy.deepcopy(q_all)
        p_atoms = copy.deepcopy(p_all_atoms)
        q_atoms = copy.deepcopy(q_all_atoms)

    else:
        p_coord = copy.deepcopy(p_all[p_view])
        q_coord = copy.deepcopy(q_all[q_view])
        p_atoms = copy.deepcopy(p_all_atoms[p_view])
        q_atoms = copy.deepcopy(q_all_atoms[q_view])

    # Recenter to centroid
    p_cent = centroid(p_coord)
    q_cent = centroid(q_coord)
    p_coord -= p_cent
    q_coord -= q_cent

    # set rotation method
    if args.rotation.lower() == METHOD_KABSCH:
        rotation_method = kabsch_rmsd
    elif args.rotation.lower() == METHOD_QUATERNION:
        rotation_method = quaternion_rmsd
    else:
        rotation_method = None

    # set reorder method
    if not args.reorder:
        reorder_method = None
    elif args.reorder_method == REORDER_QML:
        reorder_method = reorder_similarity
    elif args.reorder_method == REORDER_HUNGARIAN:
        reorder_method = reorder_hungarian
    elif args.reorder_method == REORDER_INERTIA_HUNGARIAN:
        reorder_method = reorder_inertia_hungarian
    elif args.reorder_method == REORDER_BRUTE:
        reorder_method = reorder_brute
    elif args.reorder_method == REORDER_DISTANCE:
        reorder_method = reorder_distance

    # Save the resulting RMSD
    result_rmsd = None

    if args.use_reflections:

        result_rmsd, _, _, q_review = check_reflections(
            p_atoms,
            q_atoms,
            p_coord,
            q_coord,
            reorder_method=reorder_method,
            rotation_method=rotation_method,
        )

    elif args.use_reflections_keep_stereo:

        result_rmsd, _, _, q_review = check_reflections(
            p_atoms,
            q_atoms,
            p_coord,
            q_coord,
            reorder_method=reorder_method,
            rotation_method=rotation_method,
            keep_stereo=True,
        )

    elif args.reorder:

        q_review = reorder_method(p_atoms, q_atoms, p_coord, q_coord)
        q_coord = q_coord[q_review]
        q_atoms = q_atoms[q_review]

        if not all(p_atoms == q_atoms):
            print(
                "error: Structure not aligned. "
                "Please submit bug report at "
                "http://github.com/charnley/rmsd"
            )
            sys.exit()

    # print result
    if args.output:

        if args.reorder:

            if q_review.shape[0] != q_all.shape[0]:
                print("error: Reorder length error. " "Full atom list needed for --print")
                quit()

            q_all = q_all[q_review]
            q_all_atoms = q_all_atoms[q_review]

        # Get rotation matrix
        U = kabsch(q_coord, p_coord)

        # recenter all atoms and rotate all atoms
        q_all -= q_cent
        q_all = np.dot(q_all, U)

        # center q on p's original coordinates
        q_all += p_cent

        # done and done
        xyz = set_coordinates(q_all_atoms, q_all, title=f"{args.structure_b} - modified")
        print(xyz)

    else:

        if result_rmsd:
            pass

        elif rotation_method is None:
            result_rmsd = rmsd(p_coord, q_coord)

        else:
            result_rmsd = rotation_method(p_coord, q_coord)

        print("{0}".format(result_rmsd))


if __name__ == "__main__":
    main()
