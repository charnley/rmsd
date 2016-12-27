
# Root-mean-square deviation (RMSD) of two XYZ structures.

The root-mean-square deviation (RMSD) is calculated, using Kabsch algorithm (1976) for
rotation, between two Cartesian coordinates (.xyz) or (.pdb) files.

For more information please read
[RMSD](http://en.wikipedia.org/wiki/Root-mean-square_deviation) and
[Kabsch algoritm](http://en.wikipedia.org/wiki/Kabsch_algorithm).

## Usage

Type `calculate_rmsd --help` for all the arguments.

Pretty straight forward execution, clone and run as

    calculate_rmsd molecule1.xyz molecule2.xyz

If it isn't then run it with Python 2.7

## Citation

1. Kabsch algorithm: Kabsch W., 1976, A solution for the best rotation to relate two sets of vectors, Acta Crystallographica, A32:922-923, doi:[10.1107/S0567739476001873](http://dx.doi.org/10.1107/S0567739476001873)

2. Quaternion algorithm: Michael W. Walker and Lejun Shao and Richard A. Volz, 1991, Estimating 3-D location parameters using dual number quaternions, CVGIP: Image Understanding, 54:358-367, doi:[10.1016/1049-9660(91)90036-o](http://dx.doi.org/10.1016/1049-9660\(91\)90036-o)

2. Calculate RMSD for two XYZ structures, GitHub, http://github.com/charnley/rmsd, doi: [10.5281/zenodo.46697](http://dx.doi.org/10.5281/zenodo.46697)

Note: some journals may require GitHub commit id.

## Examples

Same structure, but translated in space, so the RMSD should be zero

    calculate_rmsd examples/ethane.xyz examples/ethane_translate.xyz

You can also output (stdout) `molecule1`'s coordinates centered and rotated to
`molecule2`. Useful to visualize the difference. The output will be in XYZ
format.

    calculate_rmsd --output examples/ethane.xyz examples/ethane_translate.xyz

You can also use PDB format by using the argument `-f pdb` as seen:

    calculate_rmsd -f pdb examples/ci2_1.pdb examples/ci2_2.pdb


## Problems?

Make a issue or fork and fix it.


## Credit and Copyright

Jimmy Charnley Kromann and Lars Andersen Bratholm

