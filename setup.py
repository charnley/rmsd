#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.rst', 'rst')
except:
    long_description = ''

print long_description

quit()

short_description = 'Calculate RMSD using rotation algorithms between two molecules'

setup(name='rmsd',
      version='1.2.0',
      maintainer='Jimmy Kromann',
      maintainer_email='jimmy@charnley.dk',
      description=short_description,
      long_description=long_description,
      url='https://github.com/charnley/rmsd',
      license='BSD-2-Clause',
      packages=['rmsd'],
      install_requires=[
          'argparse',
          'numpy',
          're'
      ],
      extras_require={
      },
      classifiers=[
          "Development Status :: 5 - Production/Stable",
          "Environment :: Console",
          "Intended Audience :: End Users/Desktop",
          "License :: OSI Approved :: BSD-2-Clause",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 2.7"
      ])


