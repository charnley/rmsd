#!/usr/bin/env python

import os

import setuptools  # type: ignore

__version__ = "1.5.1"

# Find the absolute path
here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, "README.rst")) as f:
    long_description = f.read()

short_description = (
    "Calculate root-mean-square deviation (RMSD) between two "
    "sets of cartesian coordinates (XYZ or PDB format), "
    "using rotation (fx. Kabsch algorithm), "
    "atom reordering (fx. Hungarian algorithm), "
    "and axis reflections, resulting in the minimal RMSD."
)


setuptools.setup(
    name="rmsd",
    version=__version__,
    url="https://github.com/charnley/rmsd",
    python_requires=">=3.8",
    install_requires=[],
    packages=["rmsd"],
    entry_points={"console_scripts": ["calculate_rmsd=rmsd.calculate_rmsd:main"]},
    classifiers=[],
)
