#!/usr/bin/env python

""" Calculate RMSD between two XYZ files

by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Bratholm
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE

"""

import numpy
import sys
import re

def fit(P, Q):
    """ Varies the distance between P and Q, and optimizes rotation for each step
    until a minimum is found.
    """
    step_size = P.max(0)
    threshold = step_size*1e-9
    rmsd_best = kabsch(P, Q)
    while True:
        for i in range(3):
            temp = numpy.zeros(3)
            temp[i] = step_size[i]
            rmsd_new = kabsch(P+temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P[:,i] += step_size[i]
            else:
                rmsd_new = kabsch(P-temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P[:,i] -= step_size[i]
                else:
                    step_size[i] /= 2
        if (step_size<threshold).all():
            break
    return rmsd_best


def kabsch(P, Q, output=False):
    """ The Kabsch algorithm

    http://en.wikipedia.org/wiki/Kabsch_algorithm

    The algorithm starts with two sets of paired points P and Q.
    P and Q should already be centered on top of each other.

    Each vector set is represented as an NxD matrix, where D is the
    the dimension of the space.

    The algorithm works in three steps:
    - a translation of P and Q
    - the computation of a covariance matrix C
    - computation of the optimal rotation matrix U

    The optimal rotation matrix U is then used to
    rotate P unto Q so the RMSD can be caculated
    from a straight forward fashion.

    """

    # Computation of the covariance matrix
    C = numpy.dot(numpy.transpose(P), Q)

    # Computation of the optimal rotation matrix
    # This can be done using singular value decomposition (SVD)
    # Getting the sign of the det(V)*(W) to decide
    # whether we need to correct our rotation matrix to ensure a
    # right-handed coordinate system.
    # And finally calculating the optimal rotation matrix U
    # see http://en.wikipedia.org/wiki/Kabsch_algorithm
    V, S, W = numpy.linalg.svd(C)
    d = (numpy.linalg.det(V) * numpy.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]

    # Create Rotation matrix U
    U = numpy.dot(V, W)

    # Rotate P
    P = numpy.dot(P, U)

    if output:
        return P, rmsd(P, Q)

    return rmsd(P, Q)


def centroid(X):
    """ Calculate the centroid from a vectorset X """
    C = sum(X)/len(X)
    return C


def rmsd(V, W):
    """ Calculate Root-mean-square deviation from two sets of vectors V and W.
    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
    return numpy.sqrt(rmsd/N)


def write_coordinates(atoms, V, filename):
    """ Write coordinates V to filename """

    f = open(filename, 'w')

    N, D = V.shape

    f.write(str(N))
    f.write('\n')

    for i in xrange(N):
        line = "\n{0:2s} {1:10.5f} {2:10.5f} {3:10.5f}".format(atoms[i], V[i,0], V[i,1], V[i,2])
        f.write(line)

    f.close()


def get_coordinates(filename):
    """ Get coordinates from filename.

    Get coordinates from a filename.xyz and return a vectorset with all the
    coordinates.

    This function has been written to parse XYZ files, but can easily be
    written to parse others.

    """
    f = open(filename, 'r')
    V = []
    atoms = []
    n_atoms = 0
    lines_read = 0

    # Read the first line to obtain the number of atoms to read
    try:
        n_atoms = int(f.next())
    except ValueError:
        exit("Could not obtain the number of atoms in the .xyz file.")

    # Skip the title line
    f.next()

    # Use the number of atoms to not read beyond the end of a file
    for line in f:
        if lines_read == n_atoms:
            break

        atom = re.findall(r'[a-zA-Z]+', line)[0]
        numbers = re.findall(r'[-]?\d+\.\d+', line)
        numbers = [float(number) for number in numbers]

        # The numbers are not valid unless we obtain exacly three
        if len(numbers) == 3:
            V.append(numpy.array(numbers))
            atoms.append(atom)
        else:
            exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read +2))

        lines_read += 1

    f.close()
    V = numpy.array(V)
    return atoms, V


if __name__ == "__main__":

    args = sys.argv[1:]

    usage = """
usage:
calculate_rmsd [--output] <structure_a.xyz> <structure_b.xyz>

Options:
    --output  Print out the structure a, centered and transformed to structure
              b's coordinates.

Calculate Root-mean-square deviation (RMSD) between two molecules.
The two sets of atoms must be in the same order.

The script will return three RMSD values;

1) Normal: The RMSD calculated the straight-forward way.
2) Kabsch: The RMSD after the two coordinate sets are translated and rotated onto eachother.
3) Fitted: The RMSD after a fitting function has optimized the centers of the two coordinat sets.

"""

    if len(args) < 2:
        print usage
        sys.exit(0)

    output = False

    if len(args) == 3:
        output = True
        output_mol = args[2]

    mol1 = args[0]
    mol2 = args[1]

    atomsP, P = get_coordinates(mol1)
    atomsQ, Q = get_coordinates(mol2)

    print "Normal RMSD:", rmsd(P, Q)

    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    if output:
        V, r = kabsch(P, Q, output)
        print "Kabsch RMSD:", r
        V += Qc
        write_coordinates(atomsP, V, output_mol)
    else:
        print "Kabsch RMSD:", kabsch(P, Q)

    print "Fitted RMSD:", fit(P, Q)



