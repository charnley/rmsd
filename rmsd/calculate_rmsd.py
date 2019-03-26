#!/usr/bin/env python
__doc__ = \
"""
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation.

For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

__version__ = '1.3.2'

import copy
import re

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


AXIS_SWAPS = np.array([
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 1, 0],
    [2, 0, 1]])

AXIS_REFLECTIONS = np.array([
    [1, 1, 1],
    [-1, 1, 1],
    [1, -1, 1],
    [1, 1, -1],
    [-1, -1, 1],
    [-1, 1, -1],
    [1, -1, -1],
    [-1, -1, -1]])


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

    P = kabsch_rotate(P, Q, W)
    return rmsd(P, Q)


def kabsch_rotate(P, Q, W=None):
    """
    Rotate matrix P unto matrix Q using Kabsch algorithm.
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
        rotated

    """
    U = kabsch(P, Q, W)

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
    QC = centroid(Q, W)
    Q = Q - QC
    P = P - centroid(P, W)

    P = kabsch_rotate(P, Q, W) + QC
    return P

def kabsch(P, Q, W=None):
    """
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    An optional vector of weights W may be provided.

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
    W : array or None
        (N) vector, where N is points.

    Returns
    -------
    U : matrix
        Rotation matrix (D,D)
    """
    if W is None:
        C = np.dot(P.T, Q)
    else:
        # Computation of the weighted covariance matrix
        v = np.sum(W)
        try:
            P -= np.sum((P.T * W).T, axis=1, keepdims=True) / v
        except TypeError:
            P -= np.sum((P.T * W).T, axis=1, keepdims=True) // v
        try:
            Q -= np.sum((Q.T * W).T, axis=1, keepdims=True) / v
        except TypeError:
            Q -= np.sum((Q.T * W).T, axis=1, keepdims=True) // v
        C = np.dot(P.T * W, Q)

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
    W = np.asarray([
        [r4, r3, -r2, r1],
        [-r3, r4, r1, r2],
        [r2, -r1, r4, r3],
        [-r1, -r2, -r3, r4]])
    return W


def makeQ(r1, r2, r3, r4=0):
    """
    matrix involved in quaternion rotation
    """
    Q = np.asarray([
        [r4, -r3, r2, r1],
        [r3, r4, -r1, r2],
        [-r2, r1, r4, r3],
        [-r1, -r2, -r3, r4]])
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
    W_minus_Q = np.asarray([W[k] - Q[k] for k in range(N)])
    A = np.sum(Qt_dot_W, axis=0)
    eigen = np.linalg.eigh(A)
    r = eigen[1][:, eigen[0].argmax()]
    rot = quaternion_transform(r)
    return rot


def centroid(X, W=None):
    """
    Centroid is the mean position of all the points in all of the coordinate
    directions, from a vectorset X.
    An optional vector of weights W may be provided.

    https://en.wikipedia.org/wiki/Centroid

    C = sum(X)/len(X)

    Parameters
    ----------
    X : array
        (N,D) matrix, where N is points and D is dimension.
    W : array or None
        (N) vector, where N is points.

    Returns
    -------
    C : float
        centroid
    """
    if W is None:
        W = np.ones(len(X), dtype=X.dtype)
    v = np.sum(W)
    X = (X.T * W).T
    C = X.sum(axis=0) / v
    return C


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

        p_atom_idx, = np.where(p_atoms == atom)
        q_atom_idx, = np.where(q_atoms == atom)

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
    distances = cdist(A, B, 'euclidean')

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
        p_atom_idx, = np.where(p_atoms == atom)
        q_atom_idx, = np.where(q_atoms == atom)

        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]

        view = hungarian(A_coord, B_coord)
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


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
        p_atom_idx, = np.where(p_atoms == atom)
        q_atom_idx, = np.where(q_atoms == atom)

        A_coord = p_coord[p_atom_idx]
        B_coord = q_coord[q_atom_idx]

        view = brute_permutation(A_coord, B_coord)
        view_reorder[p_atom_idx] = q_atom_idx[view]

    return view_reorder


def check_reflections(p_atoms, q_atoms, p_coord, q_coord,
                      reorder_method=reorder_hungarian,
                      rotation_method=kabsch_rmsd,
                      keep_stereo=False):
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
    swap_mask = [1,-1,-1,1,-1,1]
    reflection_mask = [1,-1,-1,-1,1,1,1,-1]

    for swap, i in zip(AXIS_SWAPS, swap_mask):
        for reflection, j in zip(AXIS_REFLECTIONS, reflection_mask):
            if keep_stereo and  i * j == -1: continue # skip enantiomers

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

    fmt = "{:2s}" + (" {:15."+str(decimals)+"f}")*3

    out = list()
    out += [str(N)]
    out += [title]

    for i in range(N):
        atom = atoms[i]
        atom = atom[0].upper() + atom[1:]
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

    return


def get_coordinates(filename, fmt):
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

    return get_func(filename)


def get_coordinates_pdb(filename):
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

    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("TER") or line.startswith("END"):
                break
            if line.startswith("ATOM"):
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
                except:
                    exit("error: Parsing atomtype for the following line: \n{0:s}".format(line))

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i + 1] and "." in tokens[i + 2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("error: Parsing coordinates for the following line: \n{0:s}".format(line))
                # Try to read the coordinates
                try:
                    V.append(np.asarray(tokens[x_column:x_column + 3], dtype=float))
                except:
                    # If that doesn't work, use hardcoded indices
                    try:
                        x = line[30:38]
                        y = line[38:46]
                        z = line[46:54]
                        V.append(np.asarray([x, y ,z], dtype=float))
                    except:
                        exit("error: Parsing input for the following line: \n{0:s}".format(line))


    V = np.asarray(V)
    atoms = np.asarray(atoms)

    assert V.shape[0] == atoms.size

    return atoms, V


def get_coordinates_xyz(filename):
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

    f = open(filename, 'r')
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

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        atom = atom.upper()

        numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
        numbers = [float(number) for number in numbers]

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) >= 3:
            V.append(np.array(numbers)[:3])
            atoms.append(atom)
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

    f.close()
    atoms = np.array(atoms)
    V = np.array(V)
    return atoms, V


def main():

    import argparse
    import sys

    description = __doc__

    version_msg = """
rmsd {}

See https://github.com/charnley/rmsd for citation information

"""
    version_msg = version_msg.format(__version__)

    epilog = """
"""

    parser = argparse.ArgumentParser(
        usage='calculate_rmsd [options] FILE_A FILE_B',
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog)


    # Input structures
    parser.add_argument('structure_a', metavar='FILE_A', type=str, help='structures in .xyz or .pdb format')
    parser.add_argument('structure_b', metavar='FILE_B', type=str)

    # Admin
    parser.add_argument('-v', '--version', action='version', version=version_msg)

    # Rotation
    parser.add_argument('-r', '--rotation', action='store', default="kabsch", help='select rotation method. "kabsch" (default), "quaternion" or "none"', metavar="METHOD")

    # Reorder arguments
    parser.add_argument('-e', '--reorder', action='store_true', help='align the atoms of molecules (default: Hungarian)')
    parser.add_argument('--reorder-method', action='store', default="hungarian", metavar="METHOD", help='select which reorder method to use; hungarian (default), brute, distance')
    parser.add_argument('--use-reflections', action='store_true', help='scan through reflections in planes (eg Y transformed to -Y -> X, -Y, Z) and axis changes, (eg X and Z coords exchanged -> Z, Y, X). This will affect stereo-chemistry.')
    parser.add_argument('--use-reflections-keep-stereo', action='store_true', help='scan through reflections in planes (eg Y transformed to -Y -> X, -Y, Z) and axis changes, (eg X and Z coords exchanged -> Z, Y, X). Stereo-chemistry will be kept.')

    # Filter
    index_group = parser.add_mutually_exclusive_group()
    index_group.add_argument('-nh', '--no-hydrogen', action='store_true', help='ignore hydrogens when calculating RMSD')
    index_group.add_argument('--remove-idx', nargs='+', type=int, help='index list of atoms NOT to consider', metavar='IDX')
    index_group.add_argument('--add-idx', nargs='+', type=int, help='index list of atoms to consider', metavar='IDX')

    # format and print
    parser.add_argument('--format', action='store', help='format of input files. valid format are xyz and pdb', metavar='FMT')
    parser.add_argument('-p', '--output', '--print', action='store_true', help='print out structure B, centered and rotated unto structure A\'s coordinates in XYZ format')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # As default, load the extension as format
    if args.format is None:
        args.format = args.structure_a.split('.')[-1]

    p_all_atoms, p_all = get_coordinates(args.structure_a, args.format)
    q_all_atoms, q_all = get_coordinates(args.structure_b, args.format)

    p_size = p_all.shape[0]
    q_size = q_all.shape[0]

    if not p_size == q_size:
        print("error: Structures not same size")
        quit()

    if np.count_nonzero(p_all_atoms != q_all_atoms) and not args.reorder:
        msg = """
error: Atoms are not in the same order.

Use --reorder to align the atoms (can be expensive for large structures).

Please see --help or documentation for more information or
https://github.com/charnley/rmsd for further examples.
"""
        print(msg)
        exit()


    # Set local view
    p_view = None
    q_view = None


    if args.no_hydrogen:
        p_view = np.where(p_all_atoms != 'H')
        q_view = np.where(q_all_atoms != 'H')

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

        if args.reorder and args.output:
            print("error: Cannot reorder atoms and print structure, when excluding atoms (such as --no-hydrogen)")
            quit()

        if args.use_reflections and args.output:
            print("error: Cannot use reflections on atoms and print, when excluding atoms (such as --no-hydrogen)")
            quit()

        p_coord = copy.deepcopy(p_all[p_view])
        q_coord = copy.deepcopy(q_all[q_view])
        p_atoms = copy.deepcopy(p_all_atoms[p_view])
        q_atoms = copy.deepcopy(q_all_atoms[q_view])


    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    p_cent = centroid(p_coord)
    q_cent = centroid(q_coord)
    p_coord -= p_cent
    q_coord -= q_cent


    # set rotation method
    if args.rotation.lower() == "kabsch":
        rotation_method = kabsch_rmsd

    elif args.rotation.lower() == "quaternion":
        rotation_method = quaternion_rmsd

    elif args.rotation.lower() == "none":
        rotation_method = None

    else:
        print("error: Unknown rotation method:", args.rotation)
        quit()


    # set reorder method
    if not args.reorder:
        reorder_method = None

    if args.reorder_method == "hungarian":
        reorder_method = reorder_hungarian

    elif args.reorder_method == "brute":
        reorder_method = reorder_brute

    elif args.reorder_method == "distance":
        reorder_method = reorder_distance

    else:
        print("error: Unknown reorder method:", args.reorder_method)
        quit()


    # Save the resulting RMSD
    result_rmsd = None


    if args.use_reflections:

        result_rmsd, q_swap, q_reflection, q_review = check_reflections(
            p_atoms,
            q_atoms,
            p_coord,
            q_coord,
            reorder_method=reorder_method,
            rotation_method=rotation_method)

    elif args.use_reflections_keep_stereo:

        result_rmsd, q_swap, q_reflection, q_review = check_reflections(
            p_atoms,
            q_atoms,
            p_coord,
            q_coord,
            reorder_method=reorder_method,
            rotation_method=rotation_method,
            keep_stereo=True)

    elif args.reorder:

        q_review = reorder_method(p_atoms, q_atoms, p_coord, q_coord)
        q_coord = q_coord[q_review]
        q_atoms = q_atoms[q_review]

        if not all(p_atoms == q_atoms):
            print("error: Structure not aligned")
            quit()


    # print result
    if args.output:

        if args.reorder:

            if q_review.shape[0] != q_all.shape[0]:
                print("error: Reorder length error. Full atom list needed for --print")
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
        xyz = set_coordinates(q_all_atoms, q_all, title="{} - modified".format(args.structure_b))
        print(xyz)

    else:

        if result_rmsd:
            pass

        elif rotation_method is None:
            result_rmsd = rmsd(p_coord, q_coord)

        else:
            result_rmsd = rotation_method(p_coord, q_coord)

        print("{0}".format(result_rmsd))


    return

if __name__ == "__main__":
    main()
