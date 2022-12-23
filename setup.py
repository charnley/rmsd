#!/usr/bin/env python

import os

import setuptools  # type: ignore

__version__ = "1.5.0"

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
    maintainer="Jimmy Kromann",
    maintainer_email="jimmy@charnley.dk",
    description=short_description,
    long_description=long_description,
    url="https://github.com/charnley/rmsd",
    license="BSD-2-Clause",
    python_requires=">=3.6",
    install_requires=[
        "numpy",
        "scipy",
    ],
    packages=["rmsd"],
    entry_points={"console_scripts": ["calculate_rmsd=rmsd.calculate_rmsd:main"]},
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: End Users/Desktop",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
