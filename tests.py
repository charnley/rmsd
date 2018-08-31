#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE
"""

import os
import sys
import unittest
import numpy as np
import copy
from contextlib import contextmanager

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from rmsd import \
    get_coordinates_pdb, \
    get_coordinates_xyz, \
    get_coordinates, \
    rmsd, \
    centroid, \
    kabsch_rmsd, \
    kabsch_rotate, \
    kabsch, \
    quaternion_rmsd, \
    quaternion_rotate, \
    quaternion_transform, \
    makeQ, \
    makeW, \
    write_coordinates, \
    reorder_distance, \
    reorder_hungarian, \
    reorder_brute

@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestRMSD(unittest.TestCase):
    """Test the DSSP parser methods."""

    def setUp(self):
        """Initialize the framework for testing."""

        self.centroid = centroid
        self.rmsd = rmsd
        abs_path = os.path.abspath(os.path.dirname(__file__))
        self.xyzpath = abs_path + "/tests/"
        self.get_coordinates = get_coordinates
        self.get_coordinates_pdb = get_coordinates_pdb
        self.get_coordinates_xyz = get_coordinates_xyz

        self.kabsch_rmsd = kabsch_rmsd
        self.kabsch_rotate = kabsch_rotate
        self.kabsch_algo = kabsch

        self.quaternion_rmsd = quaternion_rmsd
        self.quaternion_rotate = quaternion_rotate
        self.quaternion_transform = quaternion_transform
        self.makeQ = makeQ
        self.makeW = makeW
        self.write_coordinates = write_coordinates

    def tearDown(self):
        """Clear the testing framework."""

        self.centroid = None
        self.rmsd = None
        self.xyzpath = None
        self.get_coordinates = None
        self.get_coordinates_pdb = None
        self.get_coordinates_xyz = None

        self.kabsch_rmsd = None
        self.kabsch_rotate = None
        self.kabsch_algo = None

        self.quaternion_rmsd = None
        self.quaternion_rotate = None
        self.quaternion_transform = None
        self.makeQ = None
        self.makeW = None
        self.write_coordinates = None

    def assertListAlmostEqual(self, list1, list2, places):
        self.assertEqual(len(list1), len(list2))
        for a, b in zip(list1, list2):
            self.assertAlmostEqual(a, b, places=places)

    def test_get_coordinates_pdb(self):
        infile = self.xyzpath + 'ci2_1.pdb'
        coords = self.get_coordinates_pdb(infile)
        self.assertEqual('N', coords[0][0])
        self.assertEqual([-7.173, -13.891, -6.266], coords[1][0].tolist())

    def test_get_coordinates_xyz(self):
        infile = self.xyzpath + 'ethane.xyz'
        coords = self.get_coordinates_xyz(infile)
        self.assertEqual('C', coords[0][0])
        self.assertEqual([-0.98353, 1.81095, -0.0314], coords[1][0].tolist())

    def test_get_coordinates(self):
        infile = self.xyzpath + 'ci2_1.pdb'
        coords = self.get_coordinates(infile, 'pdb')
        self.assertEqual('N', coords[0][0])
        self.assertEqual([-7.173, -13.891, -6.266], coords[1][0].tolist())

        infile = self.xyzpath + 'ethane.xyz'
        coords = self.get_coordinates(infile, 'xyz')
        self.assertEqual('C', coords[0][0])
        self.assertEqual([-0.98353, 1.81095, -0.0314], coords[1][0].tolist())

    def test_centroid(self):
        a1 = np.array([-19.658, 17.18, 25.163], dtype=float)
        a2 = np.array([-20.573, 18.059, 25.88], dtype=float)
        a3 = np.array([-22.018, 17.551, 26.0], dtype=float)
        atms = np.asarray([a1, a2, a3])
        centroid = self.centroid(atms)
        self.assertEqual(3, len(centroid))
        self.assertListAlmostEqual([-20.7496, 17.5966, 25.6810],
                                   centroid, places=3)

    def test_rmsd_pdb(self):
        infile1 = self.xyzpath + 'ci2_1.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        rmsd = self.rmsd(P, Q)
        self.assertAlmostEqual(26.9750, rmsd, places=3)

    def test_rmsd_xyz(self):
        infile1 = self.xyzpath + 'ethane.xyz'
        infile2 = self.xyzpath + 'ethane_mini.xyz'
        p_atoms, P = self.get_coordinates_xyz(infile1)
        q_atoms, Q = self.get_coordinates_xyz(infile2)
        rmsd = self.rmsd(P, Q)
        self.assertAlmostEqual(0.33512, rmsd, places=3)

    def test_kabash_algorith_pdb(self):
        infile1 = self.xyzpath + 'ci2_1.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        U = self.kabsch_algo(P, Q)
        self.assertListAlmostEqual([-0.5124, 0.8565, 0.0608],
                                   U[0].tolist(), places=3)

    def test_kabash_rotate_pdb(self):
        infile1 = self.xyzpath + 'ci2_1.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        nP = self.kabsch_rotate(P, Q)
        self.assertListAlmostEqual([10.6822, -2.8867, 12.6977],
                                   nP[0].tolist(), places=3)

    def test_kabash_rmsd_pdb(self):
        infile1 = self.xyzpath + 'ci2_1.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        Pc = self.centroid(P)
        Qc = self.centroid(Q)
        P -= Pc
        Q -= Qc
        rmsd = self.kabsch_rmsd(P, Q)
        self.assertAlmostEqual(11.7768, rmsd, places=3)

    def test_quaternion_rmsd_pdb(self):
        infile1 = self.xyzpath + 'ci2_1.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        Pc = self.centroid(P)
        Qc = self.centroid(Q)
        P -= Pc
        Q -= Qc
        rmsd = self.quaternion_rmsd(P, Q)
        self.assertAlmostEqual(11.7768, rmsd, places=3)

    def test_quaternion_rotate_pdb(self):
        infile1 = self.xyzpath + 'ci2_1.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        nP = self.quaternion_rotate(P, Q)
        self.assertListAlmostEqual([-0.5124, 0.8565, 0.0608],
                                   nP[0].tolist(), places=3)

    def test_quaternion_transform(self):
        r = [-0.31019, -0.59291, 0.63612, -0.38415]
        U = self.quaternion_transform(r)
        self.assertListAlmostEqual([-0.5124, 0.8565, 0.0608],
                                   U[0].tolist(), places=3)

    def test_makeQ(self):
        r = [-0.31019, -0.59291, 0.63612, -0.38415]
        Q_r = self.makeQ(*r)
        self.assertListAlmostEqual([-0.3841, -0.6361, -0.5929, -0.3101],
                                   Q_r[0].tolist(), places=3)

    def test_makeW(self):
        r = [-0.31019, -0.59291, 0.63612, -0.38415]
        Wt_r = self.makeW(*r)
        self.assertListAlmostEqual([-0.3841,  0.6361, 0.5929, -0.3101],
                                   Wt_r[0].tolist(), places=3)

    def test_write_coordinates(self):
        infile1 = self.xyzpath + 'ci2_1.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        Pc = self.centroid(P)
        Qc = self.centroid(Q)
        P -= Pc
        Q -= Qc
        U = self.kabsch_algo(P, Q)
        P -= Pc
        p_all = np.dot(P, U)

        with captured_output() as (out, err):
            self.write_coordinates(p_atoms, p_all, title="ci2_1.pdb")
        output = out.getvalue().strip().split('\n')
        self.assertEqual(output[0:4],
                         ['1064', 'ci2_1.pdb',
                          'N      10.57220149     -0.21712538     12.41498910',
                          'C       9.34616675     -0.15741025     11.60646766'])


    def test_reorder_distance(self):
        N = 5
        atoms = np.array(["H"]*N)
        p_coord = np.arange(N*3)
        p_coord = p_coord.reshape((5,3))
        q_coord = copy.deepcopy(p_coord)

        np.random.seed(6)
        np.random.shuffle(q_coord)

        review = reorder_hungarian(atoms, atoms, p_coord, q_coord)
        self.assertEqual(p_coord.tolist(), q_coord[review].tolist())

        return


    def test_reorder_brute(self):
        N = 5
        atoms = np.array(["H"]*N)
        p_coord = np.arange(N*3)
        p_coord = p_coord.reshape((5,3))
        q_coord = copy.deepcopy(p_coord)

        np.random.seed(6)
        np.random.shuffle(q_coord)

        review = reorder_brute(atoms, atoms, p_coord, q_coord)
        self.assertEqual(p_coord.tolist(), q_coord[review].tolist())

    def test_reorder_hungarian(self):
        N = 5
        atoms = np.array(["H"]*N)
        p_coord = np.arange(N*3)
        p_coord = p_coord.reshape((5,3))
        q_coord = copy.deepcopy(p_coord)

        np.random.seed(6)
        np.random.shuffle(q_coord)

        review = reorder_distance(atoms, atoms, p_coord, q_coord)
        self.assertEqual(p_coord.tolist(), q_coord[review].tolist())



if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRMSD)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(not result.wasSuccessful())
