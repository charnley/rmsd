Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation
===========================================================================

The root mean Square Deviation (RMSD) is the most common metric for measuring structural similarity between two structures. It is typically used in molecular biology, chemistry, and computational chemistry.

For more details, please see RMSD_ and `Kabsch algorithm`_.

.. _RMSD: http://en.wikipedia.org/wiki/Root-mean-square_deviation
.. _Kabsch algorithm: http://en.wikipedia.org/wiki/Kabsch_algorithm

.. contents:: Overview
    :depth: 1

Features
========

- Calculate the minimal RMSD between two molecules by applying translation and rotation.
- Supports both XYZ and PDB file formats.
- Offers atom reordering methods when input files have mismatched atom orders.
- Includes checks for conformer reflection and additional filtering options (Hydrogen, Alpha-carbon).

Motivation
==========

I want to know the minimal RMSD between two molecules
-----------------------------------------------------

To calculate the structural difference between two molecules, you might initially compute the RMSD directly (**Figure 1.a**). However, this straightforward approach could give you a misleadingly large value.
To get the true minimal RMSD, you must adjust for translations (**Figure 1.b**) and rotations (**Figure 1.c**). This process aligns the two molecules in the best possible way, ensuring the RMSD accurately reflects their structural similarity after optimal alignment.

.. list-table::
   :header-rows: 1

   * - 1.a
     - 1.b
     - 1.c

   * - |fig1.a|
     - |fig1.b|
     - |fig1.c|

   * - RMSD = 2.8
     - RMSD = 0.8
     - RMSD = 0.2

**Figure 1**: **a)** shows two molecules in space, unchanged. **b)** shows the molecules re-centered (translated) on top of each other. **c)** shows the molecules rotated to fit each other, finding the true RMSD.

I do not know the order of the atoms
------------------------------------

Atom reordering methods can be used in cases where the atoms in the two molecules are not in the same order (**Figure 2.a**). These algorithms find the optimal mapping of atoms between the two structures to minimize RMSD.

Each method has limitations because finding the best atom mapping depends on properly aligning structures. This is usually done by comparing atom-pair distances. Using the Hungarian_ cost minimization for atom distance works well if the molecules are aligned. If not, you use the molecular inertia eigenvectors (**Figure 2.b**), to rotate such the eigenvectors are aligned.
Or, use atomic descriptors (**Figure 2.c**), independent of the coordinate system, to reorder the atoms.
Note that all reordering methods have limitations and drawbacks, and the actual order might not be found.

.. _Hungarian: https://en.wikipedia.org/wiki/Hungarian_algorithm

.. list-table::
   :header-rows: 1

   * - 2.a
     - 2.b
     - 2.c

   * - |fig2.a|
     - |fig2.b|
     - |fig2.c|

**Figure 2**:
**a)** Two identical molecules, but not in the same atomic order, making it impossible to rotate correctly.
**b)** Illustrating the molecular inertia vectors used to align molecules to find the atom mapping.
**c)** Illustrating atomic representation for individual atoms. Structure-dependent but coordinate system-independent, which can be used for atom mapping.

Installation
============

Easiest is to get the program via ``pip``.

.. code-block:: bash

    pip install rmsd

There is only one Python file, so you can also download `calculate_rmsd.py` and
put it in your bin folder.

.. code-block:: bash

    wget -O calculate_rmsd https://raw.githubusercontent.com/charnley/rmsd/master/rmsd/calculate_rmsd.py
    chmod +x calculate_rmsd

Usage examples
==============

Use ``calculate_rmsd --help`` to see all the features. Usage is pretty straight
forward, call ``calculate_rmsd`` with two structures in either ``.xyz`` or
``.pdb``. In this example, Ethane has the same structure but is
translated in space, so the RMSD should be zero.

.. code-block:: bash

    calculate_rmsd tests/ethane.xyz tests/ethane_translate.xyz

It is also possible to ignore all hydrogens (useful for larger molecules where
hydrogens move around indistinguishable) and print the rotated structure for
visual comparison. The output will be in XYZ format.

.. code-block:: bash

    calculate_rmsd --no-hydrogen --print tests/ethane.xyz tests/ethane_mini.xyz

If the atoms are scrambled and not aligned, you can use the ``--reorder``
argument, which will align the atoms from structure B onto A.

Use ``--reorder-method`` to select the reordering method.
Choose between
Inertia_ aligned Hungarian_ distance ``inertia-hungarian`` (default),
Hungarian_ distance ``hungarian`` (if the structure is already aligned),
sorted distance ``distance``,
atomic representation ``qml``,
and brute force ``brute`` (for reference, don't use this).
More details on which to use in ``--help``.

.. _Hungarian: https://en.wikipedia.org/wiki/Hungarian_algorithm

.. _Inertia: https://en.wikipedia.org/wiki/Moment_of_inertia

.. code-block:: bash

    calculate_rmsd --reorder tests/water_16.xyz tests/water_16_idx.xyz

If you want to run multiple calculations simultaneously, it's best not to rely solely on the script. Instead, you can use GNU Parallel to handle this efficiently. For example, use two cores and compare all ``ethane_*`` molecules. Printing one file and the RMSD per line. Bash is good for stuff like that.

.. code-block:: bash

    find tests/resources -name "ethane_*xyz" | parallel -j2 "echo -n '{} ' && calculate_rmsd --reorder --no-hydrogen tests/resources/ethane.xyz {}"

It is also possible to use RMSD as a library in other scripts; see
``example.py`` and ``tests/*`` for example usage.


Problems?
=========

Submit issues or pull requests on GitHub.


Citation
========

Please cite this project when using it for scientific publications. And cite the relevant methods implemented.

**Implementation**:
Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation, GitHub,
http://github.com/charnley/rmsd, <git commit hash or version number>


.. list-table::
   :header-rows: 1

   * - Method
     - Argument
     - Citation

   * - **Kabsch**
     - ``--rotation-method kabsch`` (Default)
     - Wolfgang Kabsch (1976),
       Acta Crystallographica, A32:922-923

       http://dx.doi.org/10.1107/S0567739476001873

   * - **Quaternion**
     - ``--rotation-method quaternion``
     - Walker, Shao & Volz (1991),
       CVGIP: Image Understanding, 54:358-367,

       http://dx.doi.org/10.1016/1049-9660(91)90036-o

   * - **Distance Hungarian Assignment**
     - ``--reorder-method inertia-hungarian`` (Default)
     - Crouse (2016). Vol. 52, Issue 4, pp. 1679–1696, IEEE.

       http://dx.doi.org/10.1109/TAES.2016.140952

   * - **FCHL19**
     - ``--reorder-method qml``
     - Christensen et al (2020), J. Chem. Phys. 152, 044107

       https://doi.org/10.1063/1.5126701

References
==========

- https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html

A note on PDB
=============

Protein Data Bank format (PDB) is column-based; however, countless examples of non-standard ``.pdb`` files exist.
We try to read them, but if you have trouble reading the file, check if the file format is compliant with PDB.
For example, some hydrogens are noted as ``HG11``, which we assume is not mercury.

- https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM


.. |fig1.a| image:: https://raw.githubusercontent.com/charnley/rmsd/refs/heads/charnley/doc/docs/figures/fig_rmsd_nothing.png
.. |fig1.b| image:: https://raw.githubusercontent.com/charnley/rmsd/refs/heads/charnley/doc/docs/figures/fig_rmsd_recentered.png
.. |fig1.c| image:: https://raw.githubusercontent.com/charnley/rmsd/refs/heads/charnley/doc/docs/figures/fig_rmsd_rotated.png

.. |fig2.a| image:: https://raw.githubusercontent.com/charnley/rmsd/refs/heads/charnley/doc/docs/figures/fig_reorder_problem.png
.. |fig2.b| image:: https://raw.githubusercontent.com/charnley/rmsd/refs/heads/charnley/doc/docs/figures/fig_reorder_inertia.png
.. |fig2.c| image:: https://raw.githubusercontent.com/charnley/rmsd/refs/heads/charnley/doc/docs/figures/fig_reorder_qml.png
