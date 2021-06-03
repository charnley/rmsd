Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation
===========================================================================

The root-mean-square deviation (RMSD) is calculated, using Kabsch algorithm
(1976) or Quaternion algorithm (1991) for rotation, between two Cartesian
coordinates in either ``.xyz`` or ``.pdb`` format, resulting in the minimal
RMSD.

For more information please read RMSD_ and `Kabsch algorithm`_.

.. _RMSD: http://en.wikipedia.org/wiki/Root-mean-square_deviation
.. _Kabsch algorithm: http://en.wikipedia.org/wiki/Kabsch_algorithm

Motivation
----------

You have molecule A and B and want to calculate the structural difference
between those two. If you just calculate the RMSD_ straight-forward you might
get a too big of a value as seen below. You would need to first recenter the
two molecules and then rotate them unto each other to get the true minimal
RMSD. This is what this script does.

==========  ===========  ==========
No Changes  Re-centered  Rotated
----------  -----------  ----------
|begin|     |translate|  |rotate|
==========  ===========  ==========
RMSD 2.50   RMSD 1.07    RMSD 0.25
==========  ===========  ==========

.. |begin| image:: https://raw.githubusercontent.com/charnley/rmsd/master/img/plot_beginning.png
.. |translate| image:: https://raw.githubusercontent.com/charnley/rmsd/master/img/plot_translated.png
.. |rotate| image:: https://raw.githubusercontent.com/charnley/rmsd/master/img/plot_rotated.png


Citation
--------

- **Implementation**:
    Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation, GitHub,
    http://github.com/charnley/rmsd, <git commit hash or version number>

- **Kabsch algorithm**:
    Kabsch W., 1976,
    A solution for the best rotation to relate two sets of vectors,
    Acta Crystallographica, A32:922-923,
    doi: http://dx.doi.org/10.1107/S0567739476001873

- **Quaternion algorithm**:
    Michael W. Walker and Lejun Shao and Richard A. Volz, 1991,
    Estimating 3-D location parameters using dual number quaternions, CVGIP: Image Understanding, 54:358-367,
    doi: http://dx.doi.org/10.1016/1049-9660(91)90036-o

Please cite this project when using it for scientific publications.


Installation
------------

Easiest is to get the program vis PyPi under the package name ``rmsd``,

.. code-block:: bash

    pip install rmsd


or download the project from GitHub via

.. code-block:: bash

    git clone https://github.com/charnley/rmsd


There is only one Python file, so you can also download `calculate_rmsd.py` and
put it in your bin folder.

.. code-block:: bash

    wget -O calculate_rmsd https://raw.githubusercontent.com/charnley/rmsd/master/rmsd/calculate_rmsd.py
    chmod +x calculate_rmsd

Usage examples
--------------

Use ``calculate_rmsd --help`` to see all the features. Usage is pretty straight
forward, call ``calculate_rmsd`` with two structures in either ``.xyz`` or
``.pdb``. In this example Ethane has the exact same structure, but is
translated in space, so the RMSD should be zero.

.. code-block:: bash

    calculate_rmsd tests/ethane.xyz tests/ethane_translate.xyz

It is also possible to ignore all hydrogens (useful for larger molecules where
hydrogens move around indistinguishable) and print the rotated structure for
visual comparison. The output will be in XYZ format.

.. code-block:: bash

    calculate_rmsd --no-hydrogen --print tests/ethane.xyz tests/ethane_mini.xyz

If the atoms are scrambled and not aligned you can use the ``--reorder``
argument which will align the atoms from structure B unto A. Use
``--reorder-method`` to select what method for reordering. Choose between
Hungarian_ (default), distance (very approximate) and brute force (slow).

.. _Hungarian: https://en.wikipedia.org/wiki/Hungarian_algorithm

.. code-block:: bash

    calculate_rmsd --reorder tests/water_16.xyz tests/water_16_idx.xyz


It is also possible to use RMSD as a library in other scripts, see `example.py` for example usage.


Problems?
---------

Submit issues or pull requests on GitHub.


Contributions
-------------

Please note that we are using ``black`` with line length of 99. Easiest way to
abide to the code standard is to install the following package.

.. code-block:: bash

    pip install pre-commit

and run the following command in your repository

.. code-block:: bash

    pre-commit install

This will install a hook in your git and re-format your code to adhere to the standard.
As well as check for code quality.

Research notice
~~~~~~~~~~~~~~~

Please note that this repository is participating in a study into
sustainability of open source projects. Data will be gathered about this
repository for approximately the next 12 months, starting from June
2021.

Data collected will include number of contributors, number of PRs, time
taken to close/merge these PRs, and issues closed.

For more information, please visit `the informational
page <https://sustainable-open-science-and-software.github.io/>`__ or
download the `participant information
sheet <https://sustainable-open-science-and-software.github.io/assets/PIS_sustainable_software.pdf>`__.

