#!/usr/bin/env python
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

try:
    import pypandoc
    long_description = pypandoc.convert('README.md', 'rst',
                                        format='markdown_github',
                                        extra_args=("--no-wrap",))
except:
    long_description = ''

print long_description

setup(name='rmsd',
      version='1.2.0',
      maintainer='Jimmy Kromann',
      maintainer_email='jimmy@charnley.dk',
      description='Calculate RMSD between two XYZ/PDB molecules',
      long_description=long_description,
      url='https://github.com/charnley/rmsd',
      license='BSD-2-Clause',
      scripts=['calculate_rmsd.py'],
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


