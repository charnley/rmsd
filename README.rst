Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation
===========================================================================

The root mean Square Deviation (RMSD) is the most common metric for measuring structural similarity between two structures. It is typically used in molecular biology, chemistry, and computational chemistry.

However, the result can become misleadingly large unless the input data has pre-optimized translation and rotation of the molecules in question.
This solution will perform this optimization before calculating optimal (minimal) RMSD values.

Additionally, if the atoms in the molecules are not correctly ordered, optimal rotation is impossible to achieve.
This tool utilizes several ways to solve this problem.

For more details, see below and read RMSD_ and `Kabsch algorithm`_.

.. contents:: Overview
    :depth: 1

Installation
============

The easiest is to get the program via ``pip``.

.. code-block:: bash

    pip install rmsd

There is only one Python file, so you can also download `calculate_rmsd.py` and put it in your bin folder.

.. code-block:: bash

    wget -O calculate_rmsd https://raw.githubusercontent.com/charnley/rmsd/master/rmsd/calculate_rmsd.py
    chmod +x calculate_rmsd

Details
=======

To calculate the structural difference between two molecules, you might initially compute the RMSD directly (**Figure 1.a**).
However, this straightforward approach could give you a misleadingly large value.
To get the true minimal RMSD, you must adjust for translation (**Figure 1.b**) and rotation (**Figure 1.c**). This process aligns the two molecules best, ensuring the RMSD accurately reflects their structural similarity after optimal alignment.

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

Atom reordering methods are used when the atoms in two molecules are not in the same order (**Figure 2.a**).
While brute-force through all possible atom combinations and calculating the optimal rotation for each is possible, this approach is computationally infeasible for large structures, as it scales $O(N!)$.
Instead, the implemented algorithms efficiently find the optimal mapping of atoms between the two structures using smarter techniques.

Each method has limitations because finding the best atom mapping depends on properly aligning structures.
This is usually done by comparing atom-pair distances. If the molecules are aligned, using the Hungarian_ cost minimization for atom distance works well.
If not, you can align the Inertia_ eigenvectors (**Figure 2.b**) as an approximation to align the molecules.
Or, use atomic descriptors (**Figure 2.c**), independent of the coordinate system, to reorder the atoms. Note that all reordering methods have limitations and drawbacks, and the actual order might not be found.

.. list-table::
   :header-rows: 1

   * - 2.a
     - 2.b
     - 2.c

   * - |fig2.a|
     - |fig2.b|
     - |fig2.c|

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

.. code-block:: bash

    calculate_rmsd --reorder tests/water_16.xyz tests/water_16_idx.xyz

If you want to run multiple calculations simultaneously, it's best not to rely solely on the script.
Instead, you can use GNU Parallel to handle this efficiently. For example, compare all ``ethane_*`` molecules using two cores and print one file and the RMSD per line.
Bash is good for stuff like that.

.. code-block:: bash

    find tests/resources -name "ethane_*xyz" | parallel -j2 "echo -n '{} ' && calculate_rmsd --reorder --no-hydrogen tests/resources/ethane.xyz {}"

It is also possible to use RMSD as a library in other scripts; see ``tests/*`` for example usage.

Known problems
==============

Found a bug? Submit issues or pull requests on GitHub.

**Note on PDB format.** Protein Data Bank format (PDB) is column-based; however, countless examples of non-standard ``.pdb`` files exist.
We try to read them, but if you have trouble reading the file, check if the file format is compliant with PDB.
For example, some hydrogens are noted as ``HG11``, which we assume is not mercury.

- https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

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

- https://en.wikipedia.org/wiki/Root-mean-square_deviation
- https://en.wikipedia.org/wiki/Kabsch_algorithm
- https://en.wikipedia.org/wiki/Hungarian_algorithm
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html

.. _RMSD: https://en.wikipedia.org/wiki/Root-mean-square_deviation
.. _Kabsch algorithm: https://en.wikipedia.org/wiki/Kabsch_algorithm
.. _Hungarian: https://en.wikipedia.org/wiki/Hungarian_algorithm
.. _Inertia: https://en.wikipedia.org/wiki/Moment_of_inertia


.. |fig1.a| image:: https://raw.githubusercontent.com/charnley/rmsd/master/docs/figures/fig_rmsd_nothing.png
.. |fig1.b| image:: https://raw.githubusercontent.com/charnley/rmsd/master/docs/figures/fig_rmsd_recentered.png
.. |fig1.c| image:: https://raw.githubusercontent.com/charnley/rmsd/master/docs/figures/fig_rmsd_rotated.png

.. |fig2.a| image:: https://raw.githubusercontent.com/charnley/rmsd/master/docs/figures/fig_reorder_problem.png
.. |fig2.b| image:: https://raw.githubusercontent.com/charnley/rmsd/master/docs/figures/fig_reorder_inertia.png
.. |fig2.c| image:: https://raw.githubusercontent.com/charnley/rmsd/master/docs/figures/fig_reorder_qml.png
