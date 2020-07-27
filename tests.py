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
        self.reorder_inertia_hungarian = rmsd.reorder_inertia_hungarian

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

    def test_kabash_rmsd_xyz_lower_case(self):
        """Test capitalization of the first letter of the elements in the .xyz."""
        infile1 = self.xyzpath + 'lower_case_conformer_1.xyz'
        infile2 = self.xyzpath + 'lower_case_conformer_2.xyz'
        p_atoms, P = self.get_coordinates_xyz(infile1)
        q_atoms, Q = self.get_coordinates_xyz(infile2)
        Pc = self.centroid(P)
        Qc = self.centroid(Q)
        P -= Pc
        Q -= Qc
        rmsd = self.kabsch_rmsd(P, Q)
        self.assertAlmostEqual(2.1065, rmsd, places=3)

        U = self.kabsch_algo(P, Q)
        P -= Pc
        p_all = np.dot(P, U)

        with captured_output() as (out, err):
            self.print_coordinates(p_atoms, p_all, title="lower_case_conformer_1")
        output = out.getvalue().strip().split('\n')
        self.assertEqual(output[0:23],
                         ['21', 'lower_case_conformer_1',
                         'I      -1.38657986      2.39998120     -1.49939989',
                         'C      -0.98489618      1.53212019      0.38288255',
                         'C      -0.21428551      0.35944796      0.53203207',
                         'C      -1.53951234      2.18417939      1.49424119',
                         'C       0.37910374     -0.33668973     -0.61734070',
                         'C      -0.02298526     -0.12591089      1.84271696',
                         'C      -1.33498601      1.68115509      2.77593535',
                         'N       1.61296599      0.12226074     -0.96672215',
                         'C      -0.24832993     -1.38362432     -1.30101976',
                         'Cl      0.89562806     -1.56225120      2.16224578',
                         'C      -0.57784757      0.52718976      2.95019942',
                         'C       2.22824023     -0.47240505     -2.01154061',
                         'Br     -1.96729600     -2.03297392     -0.83513577',
                         'C       0.41504991     -1.97420771     -2.37599920',
                         'C       1.67604628     -1.51197881     -2.73899834',
                         'H      -2.13297589      3.08797399      1.37343959',
                         'H      -1.76400231      2.18817109      3.63678661',
                         'H      -0.41963291      0.13730923      3.95320530',
                         'H       3.21042774     -0.07689795     -2.25558754',
                         'H      -0.03959747     -2.78954454     -2.93191134',
                         'H       2.21568976     -1.95310058     -3.56993995'])

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


    def test_reorder_inertia_hungarian(self):
        # coordinates of scrambled and rotated butane
        atoms = np.array(["C","C","C","C","H","H","H","H","H","H","H","H","H","H"])
        
        p_coord = np.array(
            [[2.142e+00,  1.395e+00, -8.932e+00],
            [ 3.631e+00,  1.416e+00, -8.537e+00],
            [ 4.203e+00, -1.200e-02, -8.612e+00],
            [ 5.691e+00,  9.000e-03, -8.218e+00],
            [ 1.604e+00,  7.600e-01, -8.260e+00],
            [ 1.745e+00,  2.388e+00, -8.880e+00],
            [ 2.043e+00,  1.024e+00, -9.930e+00],
            [ 4.169e+00,  2.051e+00, -9.210e+00],
            [ 3.731e+00,  1.788e+00, -7.539e+00],
            [ 3.665e+00, -6.470e-01, -7.940e+00],
            [ 4.104e+00, -3.840e-01, -9.610e+00],
            [ 6.088e+00, -9.830e-01, -8.270e+00],
            [ 5.791e+00,  3.810e-01, -7.220e+00],
            [ 6.230e+00,  6.440e-01, -8.890e+00]])
        
        q_coord = np.array(
            [[6.71454, -5.53848, -3.50851],
            [ 6.95865, -6.22697, -2.15264],
            [ 8.16747, -5.57632, -1.45606],
            [ 5.50518, -6.19016, -4.20589],
            [ 5.33617, -5.71137, -5.14853],
            [ 7.58263, -5.64795, -4.12498],
            [ 6.51851, -4.49883, -3.35011],
            [ 6.09092, -6.11832, -1.53660],
            [ 5.70232, -7.22908, -4.36475],
            [ 7.15558, -7.26640, -2.31068],
            [ 8.33668, -6.05459, -0.51425],
            [ 7.97144, -4.53667, -1.29765],
            [ 4.63745, -6.08152, -3.58986],
            [ 9.03610, -5.68475, -2.07173]])

        p_coord -= self.centroid(p_coord)
        q_coord -= self.centroid(q_coord)

        review = self.reorder_inertia_hungarian(atoms, atoms, p_coord, q_coord)
        rmsd = self.kabsch_rmsd(p_coord, q_coord[review])
        self.assertAlmostEqual(0., rmsd, places=1)


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
