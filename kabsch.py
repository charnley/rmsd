#!/usr/bin/env python

README = """ Calculate RMSD between two XYZ files

by: Jimmy Charnley Kromann <jimmy@charnley.dk> and Lars Bratholm
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE

"""
LICENSE = """
Copyright (c) 2013, Jimmy Charnley Kromann <jimmy@charnley.dk> & Lars Bratholm
All rights reserved.

Modifications Copyright (c) 2014, Ryan G. Coleman  <ryan.g.coleman@gmail.com>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import numpy
import sys
import re

def fit(P, Q):
    """ Varies the distance between P and Q, and optimizes rotation for each 
    step until a minimum is found.
    """
    step_size = P.max(0)
    threshold = step_size*1e-9
    rmsd_best = kabsch_centered(P, Q)
    while True:
        for i in range(3):
            temp = numpy.zeros(3)
            temp[i] = step_size[i]
            rmsd_new = kabsch_centered(P+temp, Q)
            if rmsd_new < rmsd_best:
                rmsd_best = rmsd_new
                P[:,i] += step_size[i]
            else:
                rmsd_new = kabsch_centered(P-temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P[:,i] -= step_size[i]
                else:
                    step_size[i] /= 2
        if (step_size<threshold).all():
            break
    return rmsd_best

def kabsch_align_other_notnumpy(Pin, Qin, Pother):
  '''turns lists into numpy arrays, returns list as well. basically lets
  calling code ignore numpy things.'''
  PinNumpy = numpy.array(Pin)
  QinNumpy = numpy.array(Qin)
  PotherNumpy = numpy.array(Pother)
  PotherRotNumpy = kabsch_align_other(PinNumpy, QinNumpy, PotherNumpy)
  return PotherRotNumpy.tolist()

def kabsch_align_other(Pin, Qin, Pother):
  '''calls kabsch. returns aligned version of Pother  that is mapped for best 
  rmsd of Pin onto Qin.'''
  Pc = centroid(Pin)
  Qc = centroid(Qin)
  Pin -= Pc
  Qin -= Qc
  rotMatrix = kabsch_centered_find_rot(Pin, Qin)
  Pother -= Pc
  Prot = numpy.dot(Pother, rotMatrix)
  Prot += Qc   #retranslate back to orig coordinates on top of Q
  return Prot

def kabsch_align(P, Q):
  '''calls kabsch. returns aligned version of P that is mapped onto Q.'''
  Pc = centroid(P)
  Qc = centroid(Q)
  P -= Pc
  Q -= Qc
  rotMatrix = kabsch_centered_find_rot(P, Q)
  Prot = numpy.dot(P, rotMatrix)
  Prot += Qc #retranslate back to orig coordinates on top of Q
  return Prot

def kabsch(P, Q):
  '''calls kabsch_centered after moving both to the same centroid'''
  Pc = centroid(P)
  Qc = centroid(Q)
  P -= Pc
  Q -= Qc
  return kabsch_centered(P, Q)

def kabsch_centered(P, Q):
  ''' calls kabsch_centered_find_rot and gets rotation matrix U. 
  does rotation, returns rmsd'''
  Urot = kabsch_centered_find_rot(P, Q)
  # Rotate P
  P = numpy.dot(P, Urot)
  return rmsd(P, Q)

def kabsch_centered_find_rot(P, Q):
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
    if(d):
        S[-1] = -S[-1]
        V[:,-1] = -V[:,-1]
    # Create Rotation matrix U
    return numpy.dot(V, W)

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

def get_coordinates(filename):
    """ Get coordinates from filename.

    Get coordinates from a filename.xyz and return a vectorset with all the
    coordinates.

    This function has been written to parse XYZ files, but can easily be
    written to parse others.

    """
    f = open(filename, 'r')
    V = []
    # Skip the first two lines
    for _ in xrange(2):
        f.next()
    for line in f:
        numbers = re.findall(r'[-]?\d+\.\d+', line)
        numbers = [float(number) for number in numbers]
        V.append(numpy.array(numbers))
    f.close()
    V = numpy.array(V)
    return V

if __name__ == "__main__":
  args = sys.argv[1:]
  usage = """
Usage:
python calculate_rmsd.py <mol1.xyz> <mol2.xyz>

Calculate Root-mean-square deviation (RMSD) between two molecules, where the
two sets of xyz atoms are in the same order.

The script will return three RMSD values;

1) Normal: The RMSD calculated the straight-forward way (no alignment)
2) Kabsch: The RMSD after the two coordinate sets are translated and rotated onto eachother.
3) Fitted: The RMSD after a fitting function has optimized the centers of the two coordinat sets.
"""
  if len(args) < 2 and (len(args) > 1 and args[0] != 'test'):
    print usage
    sys.exit(0)
  elif args[0] == 'test':
    #run this simple test, including a bunch of output files using pdb
    import pdb
    data1 = [(59.085, 5.081, -16.853), (59.044, 4.51, -18.197),\
               (58.189, 3.235, -18.307), (58.162, 2.379, -17.426)]      
    data2 = [(59.21356, 4.898, -16.90747), (59.18413, 4.5536, -18.31091), 
               (57.86664, 3.85246, -18.59904), (56.86584, 4.09684, -17.92705)] 
    data3 = [(59.085, 5.081, -16.853), (59.044, 4.51, -18.197), \
               (58.189, 3.235, -18.307), (58.162, 2.379, -17.426), \
               (60.468, 4.19, -18.613), (61.31, 5.458, -18.547), \
               (61.052, 3.137, -17.665), (59.932, 5.289, -16.404), \
               (58.655, 5.252, -18.873), (60.462, 3.818, -19.614), \
               (60.895, 6.191, -19.211), (62.326, 5.236, -18.842), \
               (61.297, 5.839, -17.541), (61.216, 2.217, -18.206), \
               (60.363, 2.958, -16.854), (61.986, 3.494, -17.264)]
    outdata = kabsch_align_other_notnumpy(data1, data2, data3)
    pdb.debugCoords(data1, 'test.1.pdb')
    pdb.debugCoords(data2, 'test.2.pdb')
    pdb.debugCoords(data3, 'test.3.pdb')
    pdb.debugCoords(outdata, 'test.4.pdb')
    print outdata
    sys.exit(1)
  else:
    mol1 = args[0]
    mol2 = args[1]
    P = get_coordinates(mol1)
    Q = get_coordinates(mol2)
    print "Normal RMSD:", rmsd(P, Q)
    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    # http://en.wikipedia.org/wiki/Centroid
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc
    print "Kabsch RMSD:", kabsch_centered(P, Q)
    print "Fitted RMSD:", fit(P, Q)


