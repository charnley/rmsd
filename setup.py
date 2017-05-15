#!/usr/bin/env python

import setuptools
from setuptools import setup

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

short_description = 'Calculate RMSD using, translation and rotation, between molecules'

setup(name='rmsd',
      version='1.2.2',
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

