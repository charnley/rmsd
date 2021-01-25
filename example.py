#!/usr/bin/env python

from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np

import rmsd


def rotation_matrix(sigma):
    """

    https://en.wikipedia.org/wiki/Rotation_matrix

    """

    radians = sigma * np.pi / 180.0

    r11 = np.cos(radians)
    r12 = -np.sin(radians)
    r21 = np.sin(radians)
    r22 = np.cos(radians)

    R = np.array([[r11, r12], [r21, r22]])

    return R


def save_plot(A, B, filename):

    Ax = A[:, 0]
    Ay = A[:, 1]

    Bx = B[:, 0]
    By = B[:, 1]

    plt.plot(Ax, Ay, "o-", markersize=15, linewidth=3)
    plt.plot(Bx, By, "o-", markersize=15, linewidth=3)

    plt.ylim([-2.5, 2.5])
    plt.xlim([-2.5, 2.5])
    plt.grid(True)
    plt.tick_params(labelsize=15)
    plt.savefig(filename + ".png")

    plt.clf()

    return


A = np.array([[1.0, 1.0], [1.0, 2.0], [2.0, 1.5]])

# Same "molecule"
B = np.array([[1.0, 1.0], [1.0, 2.0], [2.0, 1.5]])

B *= 1.4

# Translate
B -= 3

# Rotate
B = np.dot(B, rotation_matrix(90))

print("Normal RMSD", rmsd.rmsd(A, B))
save_plot(A, B, "plot_beginning")

# Manipulate
A -= rmsd.centroid(A)
B -= rmsd.centroid(B)

print("Translated RMSD", rmsd.rmsd(A, B))
save_plot(A, B, "plot_translated")

U = rmsd.kabsch(A, B)
A = np.dot(A, U)

print("Rotated RMSD", rmsd.rmsd(A, B))
save_plot(A, B, "plot_rotated")
