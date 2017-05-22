#!/usr/bin/env python

import setuptools
from setuptools import setup

# To use a consistent encoding
from os import path

# import version
from rmsd.calculate_rmsd import __version__

# Find the absolute path
here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

short_description = 'Calculate root-mean-square deviation (RMSD), using Kabsch or Quaternion algorithm for rotation, between two Cartesian coordinates in .xyz or .pdb format, resulting in the minimal RMSD.'

setup(name='rmsd',
      version=__version__,
      maintainer='Jimmy Kromann',
      maintainer_email='jimmy@charnley.dk',
      description=short_description,
      long_description=long_description,
      url='https://github.com/charnley/rmsd',
      license='BSD-2-Clause',
      install_requires=[
          'numpy',
      ],
      packages=['rmsd'],
      entry_points={
          'console_scripts': ['calculate_rmsd=rmsd.calculate_rmsd:main']
      },
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Environment :: Console",
          "Intended Audience :: End Users/Desktop",
          "License :: OSI Approved :: BSD License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 2.7"
      ])

