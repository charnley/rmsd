#!/usr/bin/env python
__doc__ = \
"""
Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation
===========================================================================
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation. The order of the atoms *must*
be the same for both structures.
For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

__version__ = '1.2.7'

import copy
import numpy as np
import re
import scipy as scipy

# Needed for creating the distance ("cost") matrix in the Hungarian method
from scipy.spatial import distance as spatial_distance

# Needed for finding the re-ordering in the Hungarian method
from scipy.optimize import linear_sum_assignment as optimize_lsm


# Python 2/3 compatibility
# Make range a iterator in Python 2
try:
    range = xrange
except NameError:
    pass


def kabsch_rmsd(P, Q):
    """
    Rotate matrix P unto Q using Kabsch algorithm and calculate the RMSD.
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    Q : array
        (N,D) matrix, where N is points and D is dimension.
    Returns
    -------
    rmsd : float
        root-mean squared deviation
    """
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


def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix
    P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
    calculated.
    Using the Kabsch algorithm with two sets of paired point P and Q, centered
    around the centroid. Each vector set is represented as an NxD
    matrix, where D is the the dimension of the space.
    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U
    http://en.wikipedia.org/wiki/Kabsch_algorithm
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
    Example
    -----
    TODO
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


def quaternion_rmsd(P, Q):
    """
    Rotate matrix P unto Q and calculate the RMSD
    based on doi:10.1016/1049-9660(91)90036-O
    Parameters
    ----------
    P : array
        (N,D) matrix, where N is points and D is dimension.
    P : array
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


def centroid(X):
    """
    Calculate the centroid from a vectorset X.
    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.
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


def reorder_distance(atoms, X, debug = False):
    """
    Re-orders the input atom list and xyz coordinates by atom type and then 
    by distance of each atom from the centroid
    Parameters
    ----------
    atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    X : array
        (N,D) matrix, where N is points and D is dimension
    Returns
    -------
    atoms_reordered : array
             (N,1) matrix, where N is points holding the ordered atoms' names
    coords_reordered : array
             (N,D) matrix, where N is points and D is dimension (rows re-ordered)
    """

    # Calculate distance from each atom to centroid
    norms = np.linalg.norm(X, axis=1)

    if debug:
        # unsorted array print
        print('------------------')
        print('Original atoms: {}'.format(atoms))
        print('Original norms: {}'.format(norms))
        print('Original array: {}'.format(X))

    # Array indices for re-ordering atoms
    reorder_indices = np.lexsort((norms, atoms))

    if debug:
        print('Sorted indices of original atoms, then norms: {}'.format(reorder_indices))

    # Calculates length (# of atoms) for new atom and coord matrics
    num_atoms = reorder_indices.shape[0]

    # Create new arrays for the atoms, norms (not really necessary), and coordinate matrix
    atoms_ordered = np.zeros(num_atoms,  'U2')
    norms_ordered = np.zeros(num_atoms, dtype = float)
    coords_ordered = np.zeros((num_atoms, 3), dtype = float)

    # Re-order the atom array, norms (not really necessary), and coordinate matrix
    atoms_ordered = atoms[reorder_indices]
    norms_ordere = norms[reorder_indices]
    coords_ordered = X[reorder_indices]


    if debug:
        for (i, j, k) in zip(atoms_ordered, norms_ordered, coords_ordered):
            print('%3s %8.5f %s' % (i, j, k))
        print('\n')

    return (atoms_ordered, coords_ordered)


def generate_permutations(elements, n):
    # Heap's algorithm for generating all n! permutations in a list
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


def reorder_hungarian(atoms, A, B, debug = False, debug_abridged = False):
    """
    Re-orders the input atom list and xyz coordinates using the Hungarian method (using optimized column results)
    Parameters
    ----------
    atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    A : array
        (N,D) matrix, where N is points and D is dimension
    B : array
        (N,D) matrix, where N is points and D is dimension
    reorder_type : string
        Contains the method for re-ordering, either "brute" or "hungarian"
    Returns
    -------
    atoms_min : array
             (N,1) matrix, where N is points holding the ordered atoms' names
    coords_min : array
             (N,D) matrix, where N is points and D is dimension (rows re-ordered)
    """

    swaps = np.array([[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 1, 0], [2, 0, 1]])
    reflections = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])

    rmsd_min = 100000000

    # Calculates length (# of atoms) for new atom and coord matrics
    num_atoms = A.shape[0]
    num_swaps = swaps.shape[0]
    num_reflections = reflections.shape[0]

    atoms_min = np.zeros(num_atoms,  'U2')
    coords_min = np.zeros((num_atoms, 3), dtype = float)

    # Sets initial ordering for row indices to [0, 1, 2, ..., len(A)], used in brute-force method
    initial_order = list(range(num_atoms))


    # Iterate over all swaps (eg X and Z coords exchanged -> Z, Y, X)
    for r in range(0, num_swaps):

        # Iterate over all reflections (eg Y transformed to -Y -> X, -Y, Z)
        for s in range(0, num_reflections):

            # Create blank arrays for the re-ordered atoms and coordinates
            atoms_ordered = np.zeros(num_atoms,  'U2')
            coords_ordered = np.zeros((num_atoms, 3), dtype = float)

            # Create blank array for trial coordiantes
            coords_temp = np.zeros((num_atoms, 3), dtype = float)

            # Apply coordinate swap r and coordinate reflection s to 2nd structure
            coords_temp[:, range(0, 3)] = np.dot(B[:, swaps[r, range(0, 3)]], np.diag(reflections[s])) 

            # Create distance matrix between atoms in structure 1 and 2
            distances = spatial_distance.cdist(A, coords_temp, 'euclidean')

            # unsorted array print
            if debug:
                print('------------------')
                print('Swaps: {}'.format(swaps[r]))
                print('Reflections: {}'.format(reflections[s]))
                print('Original atoms: {}'.format(atoms))
                print('Original array: {}'.format(B))
                print('New array: {}'.format(coords_temp))
                print('Distance matrix: {}'.format(distances))

            # Perform Hungarian analysis on distance matrix between atoms of 1st structure and trial structure
            row_indices, reorder_indices = optimize_lsm(distances)

            # Re-order the atom array and coordinate matrix
            atoms_ordered = atoms[reorder_indices]
            coords_ordered = coords_temp[reorder_indices]

            if debug:
                print('New assigned row: {}'.format(reorder_indices)) 
                for (i, j) in zip(atoms_ordered, coords_ordered):
                    print("%3s  %s" % (i, j))

            # Calculate the RMSD between structure 1 and the Hungarian re-ordered structure 2
            rmsd_temp = quaternion_rmsd(A, coords_ordered)

            if debug:
                print('temp RMSD = %10.5e, min RMSD = %10.5e\n' % (rmsd_temp, rmsd_min))

            # Replaces the atoms and coordinates with the current structure if the RMSD is lower
            if rmsd_temp < rmsd_min:
                swap_min = swaps[r]
                reflection_min = reflections[s]
                atoms_min = atoms_ordered
                rmsd_min = rmsd_temp
                coords_min = coords_ordered

                if debug_abridged:
                    print('------------------')
                    print('Swaps: {}'.format(swaps[r]))
                    print('Reflections: {}'.format(reflections[s]))
                    print('New row indices: {}'.format(reorder_indices))
                    print('temp RMSD = %10.5e, min RMSD = %10.5e\n' % (rmsd_temp, rmsd_min))

    return (atoms_min, coords_min)



def reorder_brute(atoms, A, B, debug = False, debug_abridged = False):
    """
    Re-orders the input atom list and xyz coordinates using the brute force method of permuting all rows of the input coordinates
    Parameters
    ----------
    atoms : array
        (N,1) matrix, where N is points holding the atoms' names
    A : array
        (N,D) matrix, where N is points and D is dimension
    B : array
        (N,D) matrix, where N is points and D is dimension
    reorder_type : string
        Contains the method for re-ordering, either "brute" or "hungarian"
    Returns
    -------
    atoms_min : array
             (N,1) matrix, where N is points holding the ordered atoms' names
    coords_min : array
             (N,D) matrix, where N is points and D is dimension (rows re-ordered)
    """

    swaps = np.array([[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 1, 0], [2, 0, 1]])
    reflections = np.array([[1, 1, 1], [-1, 1, 1], [1, -1, 1], [1, 1, -1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1], [-1, -1, -1]])

    rmsd_min = 100000000

    # Calculates length (# of atoms) for new atom and coord matrics
    num_atoms = A.shape[0]
    num_swaps = swaps.shape[0]
    num_reflections = reflections.shape[0]

    atoms_min = np.zeros(num_atoms,  'U2')
    coords_min = np.zeros((num_atoms, 3), dtype = float)

    # Sets initial ordering for row indices to [0, 1, 2, ..., len(A)], used in brute-force method
    initial_order = list(range(num_atoms))


    # Iterate over all swaps (eg X and Z coords exchanged -> Z, Y, X)
    for r in range(0, num_swaps):

        # Iterate over all reflections (eg Y transformed to -Y -> X, -Y, Z)
        for s in range(0, num_reflections):

            # Create blank arrays for the re-ordered atoms and coordinates
            atoms_ordered = np.zeros(num_atoms,  'U2')
            coords_ordered = np.zeros((num_atoms, 3), dtype = float)

            # Create blank array for trial coordiantes
            coords_temp = np.zeros((num_atoms, 3), dtype = float)

            # Apply coordinate swap r and coordinate reflection s to 2nd structure
            coords_temp[:, range(0, 3)] = np.dot(B[:, swaps[r, range(0, 3)]], np.diag(reflections[s])) 

            for reorder_indices in generate_permutations(initial_order, num_atoms):

                # unsorted array print
                if debug:
                    print('Swaps: {}'.format(swaps[r]))
                    print('Reflections: {}'.format(reflections[s]))
                    print('Original atoms: {}'.format(atoms))
                    print('Original array: {}'.format(B))
                    print('New array: {}'.format(coords_temp))


                # Re-order the atom array and coordinate matrix
                atoms_ordered = atoms[reorder_indices]
                coords_ordered = coords_temp[reorder_indices]

                if debug:
                    print('New assigned row: {}'.format(reorder_indices)) 
                    for (i, j) in zip(atoms_ordered, coords_ordered):
                        print("%3s  %s" % (i, j))

                # Calculate the RMSD between structure 1 and the Hungarian re-ordered structure 2
                rmsd_temp = quaternion_rmsd(A, coords_ordered)

                if debug:
                    print('temp RMSD = %10.5f, min RMSD = %10.5e\n' % (rmsd_temp, rmsd_min))

                # Replaces the atoms and coordinates with the current structure if the RMSD is lower
                if rmsd_temp < rmsd_min:
                    swap_min = swaps[r]
                    reflection_min = reflections[s]
                    atoms_min = atoms_ordered
                    rmsd_min = rmsd_temp
                    coords_min = coords_ordered

                    if debug_abridged:
                        print('------------------')
                        print('Swaps: {}'.format(swaps[r]))
                        print('Reflections: {}'.format(reflections[s]))
                        print('New row indices: {}'.format(reorder_indices))
                        print('temp RMSD = %10.5e, min RMSD = %10.5e\n' % (rmsd_temp, rmsd_min))

    return (atoms_min, coords_min)


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
        Root-mean-square deviation
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)


def write_coordinates(atoms, V, title=""):
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
    """
    N, D = V.shape

    print(str(N))
    print(title)

    for i in range(N):
        atom = atoms[i]
        atom = atom[0].upper() + atom[1:]
        print("{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}".format(
                atom, V[i, 0], V[i, 1], V[i, 2]))


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
        return get_coordinates_xyz(filename)
    elif fmt == "pdb":
        return get_coordinates_pdb(filename)
    exit("Could not recognize file format: {:s}".format(fmt))


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
                        exit("Error parsing atomtype for the following line: \n{0:s}".format(line))

                if x_column == None:
                    try:
                        # look for x column
                        for i, x in enumerate(tokens):
                            if "." in x and "." in tokens[i + 1] and "." in tokens[i + 2]:
                                x_column = i
                                break
                    except IndexError:
                        exit("Error parsing coordinates for the following line: \n{0:s}".format(line))
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
                        exit("Error parsing input for the following line: \n{0:s}".format(line))


    V = np.asarray(V)
    atoms = np.asarray(atoms)
    assert(V.shape[0] == atoms.size)
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
        exit("Could not obtain the number of atoms in the .xyz file.")

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
        if len(numbers) == 3:
            V.append(np.array(numbers))
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

    description = """
Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
or PDB format, using transformation and rotation. The order of the atoms *must*
be the same for both structures.
For more information, usage, example and citation read more at
https://github.com/charnley/rmsd
"""

    epilog = """output:
  Normal - RMSD calculated the straight-forward way, no translation or rotation.
  Kabsch - RMSD after coordinates are translated and rotated using Kabsch.
  Quater - RMSD after coordinates are translated and rotated using quaternions.
"""

    parser = argparse.ArgumentParser(
                    usage='%(prog)s [options] structure_a structure_b',
                    description=description,
                    formatter_class=argparse.RawDescriptionHelpFormatter,
                    epilog=epilog)

    parser.add_argument('-v', '--version', action='version', version='rmsd ' + __version__ + "\nhttps://github.com/charnley/rmsd")

    parser.add_argument('structure_a', metavar='structure_a', type=str, help='Structure in .xyz or .pdb format')
    parser.add_argument('structure_b', metavar='structure_b', type=str)
    parser.add_argument('-o', '--output', action='store_true', help='print out structure A, centered and rotated unto structure B\'s coordinates in XYZ format')
    parser.add_argument('-f', '--format', action='store', help='Format of input files. Valid format are XYZ and PDB', metavar='fmt')
    parser.add_argument('-s', '--reorder', action='store_true', help='Uses Hungarian method for re-ordering atoms for improved alignment')
    parser.add_argument('-sb', '--reorder_brute', action='store_true', help='Uses brute-force method for re-ordering atoms for optimal alignment  (only use if # atoms < 10 or you enjoy twiddling your thumbs)')
    parser.add_argument('-sd', '--reorder_distance', action='store_true', help='Uses distance from centroid to re-order atoms for improved alignment')

    parser.add_argument('-m', '--normal', action='store_true', help='Use no transformation')
    parser.add_argument('-k', '--kabsch', action='store_true', help='Use Kabsch algorithm for transformation')
    parser.add_argument('-q', '--quater', action='store_true', help='Use Quaternion algorithm for transformation')

    index_group = parser.add_mutually_exclusive_group()
    index_group.add_argument('-n', '--no-hydrogen', action='store_true', help='ignore hydrogens when calculating RMSD')
    index_group.add_argument('-r', '--remove-idx', nargs='+', type=int, help='index list of atoms NOT to consider', metavar='idx')
    index_group.add_argument('-a', '--add-idx', nargs='+', type=int, help='index list of atoms to consider', metavar='idx')


    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # As default use all three methods
    if args.output:
        args.normal = False
        args.kabsch = False
        args.quater = False
    elif not args.normal and not args.kabsch and not args.quater:
        args.normal = True
        args.kabsch = True
        args.quater = True

    # As default, load the extension as format
    if args.format == None:
        args.format = args.structure_a.split('.')[-1]

    p_atoms, p_all = get_coordinates(args.structure_a, args.format)
    q_atoms, q_all = get_coordinates(args.structure_b, args.format)


    if np.count_nonzero(p_atoms != q_atoms) and not args.reorder and not args.reorder_brute and not args.order_distance:
        exit("Atoms not in the same order")

    if args.no_hydrogen:
        not_hydrogens = np.where(p_atoms != 'H')
        P = copy.deepcopy(p_all[not_hydrogens])
        Q = copy.deepcopy(q_all[not_hydrogens])

    elif args.remove_idx:
        N, = p_atoms.shape
        index = range(N)
        index = set(index) - set(args.remove_idx)
        index = list(index)
        P = copy.deepcopy(p_all[index])
        Q = copy.deepcopy(q_all[index])

    elif args.add_idx:
        P = copy.deepcopy(p_all[args.add_idx])
        Q = copy.deepcopy(q_all[args.add_idx])

    else:
        P = copy.deepcopy(p_all)
        Q = copy.deepcopy(q_all)


    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    
    # Calculate 'dumb' RMSD
    if args.normal and not args.output:
        normal_rmsd = rmsd(P, Q)
        print("Normal RMSD: {0}".format(normal_rmsd))

    if args.reorder and args.format == 'xyz':
       (q_atoms, Q) = reorder_hungarian(q_atoms, P, Q)

    elif args.reorder_brute and args.format == 'xyz':
       (q_atoms, Q) = reorder_brute(q_atoms, P, Q)

    elif args.reorder_distance and args.format == 'xyz':
       (p_atoms, P) = reorder_distance(p_atoms, P)
       (q_atoms, Q) = reorder_distance(q_atoms, Q)


    if args.kabsch:
        print("Kabsch RMSD: {0}".format(kabsch_rmsd(P, Q)))

    if args.quater:
        print("Quater RMSD: {0}".format(quaternion_rmsd(P, Q)))

    if args.output:
        U = kabsch(P, Q)
        p_all -= Pc
        p_all = np.dot(p_all, U)
        p_all += Qc
        write_coordinates(p_atoms, p_all, title="{} translated".format(args.structure_a))

    return

if __name__ == "__main__":
    main()
