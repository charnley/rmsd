Root-mean-square deviation (RMSD) of two XYZ structures.
====

The root-mean-square deviation (RMSD) is calculated, using Kabsch algorithm (1976) for
rotation, between two Cartesian coordinates (.xyz) files.

For more information please read
[RMSD](http://en.wikipedia.org/wiki/Root-mean-square_deviation) and
[Kabsch algoritm](http://en.wikipedia.org/wiki/Kabsch_algorithm).

## Usage

The code should be executable, so run it via terminal:

    ./calculate_rmsd molecule1.xyz molecule2.xyz

If it isn't then run it with Python 2.7

    python calculate_rmsd molecule1.xyz molecule2.xyz

## Citation

1. Kabsch W., 1976, A solution for the best rotation to relate two sets of vectors, Acta Crystallographica, A32:922-923, doi:[10.1107/S0567739476001873](http://dx.doi.org/10.1107/S0567739476001873)

2. GitHub: Calculate RMSD for two XYZ structures, http://github.com/charnley/rmsd

Note: some journals may require GitHub commit id.

## Examples

Same structure, but translated in space

    ./calculate_rmsd examples/ethane.xyz examples/ethane_trans.xyz

Same structure, two forcefield minimizations

    ./calculate_rmsd examples/ethane.xyz examples/ethane_mini.xyz

Same structure, different bondlengths

    ./calculate_rmsd examples/ethane.xyz examples/ethane_bond.xyz

You can also output (stdout) `molecule1`'s coordinates centered and rotated to
`molecule2`. Useful to visualize the difference.

    ./calculate_rmsd --output examples/ethane.xyz examples/ethane_trans.xyz

## Help

Type `./calculate_rmsd --help` or look in the source file.

## Credit and Copyright

Jimmy Charnley Kromann and Lars Bratholm

