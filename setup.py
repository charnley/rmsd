#!/usr/bin/env python

import setuptools
import os

__version__ = '1.2.7'

# Find the absolute path
here = os.path.abspath(os.path.dirname(__file__))

# Get the long description from the README file
with open(os.path.join(here, 'README.rst')) as f:
    long_description = f.read()

short_description = 'Calculate root-mean-square deviation (RMSD), using Kabsch or Quaternion algorithm for rotation, between two Cartesian coordinates in .xyz or .pdb format, resulting in the minimal RMSD.'

setuptools.setup(name='rmsd',
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

