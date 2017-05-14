
Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation
===========================================================================

The root-mean-square deviation (RMSD) is calculated, using Kabsch algorithm
(1976) or Quaternion algorithm (1991) for rotation, between two Cartesian
coordinates in either ``.xyz`` or ``.pdb`` format, resulting in the prober
minimal RMSD.

For more information please read RMSD_ and `Kabsch algorithm`_.

.. _RMSD: http://en.wikipedia.org/wiki/Root-mean-square_deviation
.. _Kabsch algorithm: http://en.wikipedia.org/wiki/Kabsch_algorithm

Installation
------------

From PyPI
~~~~~~~~~


.. code-block:: bash

    $ pip install rmsd


Download Python file
~~~~~~~~~~~~~~~~~~~~

You can also just download the Python file, put it in your bin folder and make it executable.

.. code-block:: bash

    $ wget -O calculate_rmsd https://raw.githubusercontent.com/charnley/rmsd/master/calculate_rmsd.py
    $ chmod +x calculate_rmsd

Usage
-----

Type ``calculate_rmsd --help`` for all the arguments.

Pretty straight forward execution, clone and run as

.. code-block:: bash

    $ calculate_rmsd molecule1.xyz molecule2.xyz

or

.. code-block:: bash

    $ calculate_rmsd protein.pdb protein.pdb


Citation
--------

    - **Kabsch algorithm**:
      Kabsch W., 1976,
      A solution for the best rotation to relate two sets of vectors,
      Acta Crystallographica, A32:922-923,
      doi: http://dx.doi.org/10.1107/S0567739476001873

    - **Quaternion algorithm**:
      Michael W. Walker and Lejun Shao and Richard A. Volz, 1991,
      Estimating 3-D location parameters using dual number quaternions, CVGIP: Image Understanding, 54:358-367,
      doi: http://dx.doi.org/10.1016/1049-9660(91)90036-o

    - **Implementation**:
      Calculate RMSD for two XYZ structures, GitHub,
      http://github.com/charnley/rmsd, <commit hash or version number>

Please remmeber to cite this implementation when using it for scientific publications.

More usage examples
-------------------

Same structure, but translated in space, so the RMSD should be zero

.. code-block:: bash

    calculate_rmsd examples/ethane.xyz examples/ethane_translate.xyz

You can also output (stdout) ``molecule1``'s coordinates centered and rotated to
``molecule2``. Useful to visualize the difference. The output will be in XYZ
format.

.. code-block:: bash

    calculate_rmsd --output examples/ethane.xyz examples/ethane_translate.xyz

You can also use PDB format by using the argument `-f pdb` as seen:

.. code-block:: bash

    calculate_rmsd -f pdb examples/ci2_1.pdb examples/ci2_2.pdb

Problems?
---------

Submit issues on GitHub or submit pull requests.

Credit and Copyright
--------------------

Jimmy Charnley Kromann and Lars Andersen Bratholm

