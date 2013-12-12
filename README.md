Root-mean-square deviation (RMSD) of two XYZ structures.
====

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/de4affc944cd23d22da90d39df19af82 "githalytics.com")](http://githalytics.com/charnley/rmsd)

The root-mean-square deviation (RMSD) is calculated, using Kabsch algorithm (1976) for
rotation, between two cartesian coordinates (.xyz) files.

Based on
[RMSD](http://en.wikipedia.org/wiki/Root-mean-square_deviation) and
[Kabsch algoritm](http://en.wikipedia.org/wiki/Kabsch_algorithm).

## Usage

    python calculate_rmsd.py molecule1.xyz molecule2.xyz

## Examples

Same molecule, but translated in space

    python calculate_rmsd.py examples/ethane.xyz examples/ethane_trans.xyz

Same molecule, two forcefield minimizations

    python calculate_rmsd.py examples/ethane.xyz examples/ethane_mini.xyz

Same structure, different bondlengths

    python calculate_rmsd.py examples/ethane.xyz examples/ethane_bond.xyz

## Credit and Copyright

Jimmy Charnley Kromann and Lars Bratholm

