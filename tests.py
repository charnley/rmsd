#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
project: https://github.com/charnley/rmsd
license: https://github.com/charnley/rmsd/blob/master/LICENSE
"""

import copy
import os
import sys
import unittest
from contextlib import contextmanager

import numpy as np
import rmsd

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


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

        abs_path = os.path.abspath(os.path.dirname(__file__))
        self.xyzpath = abs_path + "/tests/"

        self.centroid = rmsd.centroid
        self.rmsd = rmsd.rmsd

        self.get_coordinates = rmsd.get_coordinates
        self.get_coordinates_pdb = rmsd.get_coordinates_pdb
        self.get_coordinates_xyz = rmsd.get_coordinates_xyz

        self.kabsch_rmsd = rmsd.kabsch_rmsd
        self.kabsch_rotate = rmsd.kabsch_rotate
        self.kabsch_fit = rmsd.kabsch_fit
        self.kabsch_algo = rmsd.kabsch

        self.quaternion_rmsd = rmsd.quaternion_rmsd
        self.quaternion_rotate = rmsd.quaternion_rotate
        self.quaternion_transform = rmsd.quaternion_transform
        self.makeQ = rmsd.makeQ
        self.makeW = rmsd.makeW
        self.print_coordinates = rmsd.print_coordinates

        self.reorder_brute = rmsd.reorder_brute
        self.reorder_hungarian = rmsd.reorder_hungarian
        self.reorder_distance = rmsd.reorder_distance

        self.check_reflections = rmsd.check_reflections

    def tearDown(self):
        """Clear the testing framework."""

        self.xyzpath = None

        self.centroid = None
        self.rmsd = None
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
        self.print_coordinates = None

        self.reorder_brute = None
        self.reorder_hungarian = None
        self.reorder_distance = None

        self.check_reflections = None


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

    def test_kabash_fit_pdb(self):
        infile1 = self.xyzpath + 'ci2_1r+t.pdb'
        infile2 = self.xyzpath + 'ci2_1.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        nP = self.kabsch_fit(P, Q)
        self.assertListAlmostEqual(Q[0].tolist(),
                                   nP[0].tolist(), places=2)

    def test_kabash_weighted_fit_pdb(self):
        infile1 = self.xyzpath + 'ci2_12.pdb'
        infile2 = self.xyzpath + 'ci2_2.pdb'
        p_atoms, P = self.get_coordinates_pdb(infile1)
        q_atoms, Q = self.get_coordinates_pdb(infile2)
        weights = np.zeros(len(P))
        residue13_start = 200 
        residue24_start = 383
        weights[residue13_start:residue24_start] = 1.0
        nP = self.kabsch_fit(P, Q, weights)
        self.assertListAlmostEqual(Q[300].tolist(),
                                   nP[300].tolist(), places=2)

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

    def test_print_coordinates(self):
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
            self.print_coordinates(p_atoms, p_all, title="ci2_1.pdb")
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

        review = self.reorder_hungarian(atoms, atoms, p_coord, q_coord)
        self.assertEqual(p_coord.tolist(), q_coord[review].tolist())

        return


    def test_reorder_brute(self):
        N = 5
        atoms = np.array(["H"]*N)
        p_coord = np.arange(N*3)
        p_coord = p_coord.reshape((N,3))
        q_coord = copy.deepcopy(p_coord)

        np.random.seed(6)
        np.random.shuffle(q_coord)

        review = self.reorder_brute(atoms, atoms, p_coord, q_coord)
        self.assertEqual(p_coord.tolist(), q_coord[review].tolist())


    def test_reorder_brute_ch(self):

        N = 6
        p_atoms = ["C"]*3 + ["H"]*3
        p_atoms = np.array(p_atoms)
        p_coord = np.arange(N*3)
        p_coord = p_coord.reshape((N,3))

        # random index
        np.random.seed(6)
        idx = np.arange(N, dtype=int)
        np.random.shuffle(idx)

        q_coord = copy.deepcopy(p_coord)
        q_atoms = copy.deepcopy(p_atoms)

        q_coord = q_coord[idx]
        q_atoms = q_atoms[idx]

        review = self.reorder_brute(p_atoms, q_atoms, p_coord, q_coord)

        self.assertEqual(p_coord.tolist(), q_coord[review].tolist())
        self.assertEqual(p_atoms.tolist(), q_atoms[review].tolist())


    def test_reorder_hungarian(self):
        N = 5
        atoms = np.array(["H"]*N)
        p_coord = np.arange(N*3)
        p_coord = p_coord.reshape((5,3))
        q_coord = copy.deepcopy(p_coord)

        np.random.seed(6)
        np.random.shuffle(q_coord)

        review = self.reorder_distance(atoms, atoms, p_coord, q_coord)
        self.assertEqual(p_coord.tolist(), q_coord[review].tolist())


    def test_reflections(self):

        atoms = np.array(["C", "H", "H", "H", "F"])

        p_coord = np.array(
            [[-0.000000,  -0.000000,  -0.000000],
            [  1.109398,  -0.000000,   0.000000],
            [ -0.3697920, -0.7362220, -0.7429600],
            [ -0.3698020,  1.011538,  -0.2661100],
            [ -0.3698020, -0.2753120,  1.009070]])
        q_coord = copy.deepcopy(p_coord)
        q_coord[:,[0, 2]] = q_coord[:,[2, 0]]

        min_rmsd, min_swap, min_reflection, min_review = self.check_reflections(atoms, atoms, p_coord, q_coord, reorder_method=None)

        assert np.isclose(min_rmsd, 0.0, atol=1e-6)


    def test_reflections_norotation(self):

        atoms = np.array(["C", "H", "H", "H", "F"])

        p_coord = np.array(
            [[-0.000000,  -0.000000,  -0.000000],
            [  1.109398,  -0.000000,   0.000000],
            [ -0.3697920, -0.7362220, -0.7429600],
            [ -0.3698020,  1.011538,  -0.2661100],
            [ -0.3698020, -0.2753120,  1.009070]])
        q_coord = copy.deepcopy(p_coord)
        q_coord[:,[0, 2]] = q_coord[:,[2, 0]]

        min_rmsd, min_swap, min_reflection, min_review = self.check_reflections(atoms, atoms, p_coord, q_coord, reorder_method=None, rotation_method=None)

        assert np.isclose(min_rmsd, 0.0, atol=1e-6)


    def test_reflections_reorder(self):

        p_atoms = np.array(["C", "H", "H", "H", "F"])

        p_coord = np.array(
            [[-0.000000,  -0.000000,  -0.000000],
            [  1.109398,  -0.000000,   0.000000],
            [ -0.3697920, -0.7362220, -0.7429600],
            [ -0.3698020,  1.011538,  -0.2661100],
            [ -0.3698020, -0.2753120,  1.009070]])

        # random reflection
        q_coord = copy.deepcopy(p_coord)
        q_coord[:,[0, 2]] = q_coord[:,[2, 0]]

        # random order
        # review = list(range(len(atoms)))
        review = [1, 4, 2, 0, 3]
        q_coord = q_coord[review]
        q_atoms = p_atoms[review]

        min_rmsd, min_swap, min_reflection, min_review = self.check_reflections(p_atoms, q_atoms, p_coord, q_coord)

        assert np.isclose(min_rmsd, 0.0, atol=1e-6)

    def test_reflections_keep_stereo(self):

        atoms = np.array(["C", "H", "H", "H", "F"])

        p_coord = np.array(
            [[-0.000000,  -0.000000,  -0.000000],
            [  1.109398,  -0.000000,   0.000000],
            [ -0.3697920, -0.7362220, -0.7429600],
            [ -0.3698020,  1.011538,  -0.2661100],
            [ -0.3698020, -0.2753120,  1.009070]])
        q_coord = copy.deepcopy(p_coord)
        q_coord[:,[0, 2]] = q_coord[:,[2, 0]] # swap [2,1,0]: prepare enantiomer coordinates (which is named as q_coord) of p_coord

        # If keep_stereo is off, enantiomer coordinates of q_coord are considered, resulting into identical coordinates of p_coord.
        min_rmsd, min_swap, min_reflection, min_review = self.check_reflections(atoms, atoms, p_coord, q_coord,
                                                                                reorder_method=None, keep_stereo=False)
        assert np.isclose(min_rmsd, 0.0, atol=1e-6) # the enantiomer of the enantiomer of the original molecule has zero RMSD with the original molecule.

        # No enantiomer coordinates, non-zero RMSD.
        min_rmsd, min_swap, min_reflection, min_review = self.check_reflections(atoms, atoms, p_coord, q_coord,
                                                                                reorder_method=None, keep_stereo=True)
        assert np.isclose(min_rmsd, 1.1457797, atol=1e-6)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRMSD)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    sys.exit(not result.wasSuccessful())
