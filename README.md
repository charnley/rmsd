Root-mean-square deviation of two XYZ files based on Kabsch algorithm.
====

The root-mean-square deviation (RMSD) is calculated, using Kabsch algorithm for rotation, between two XYZ files.

Based on:[Root-mean-square deviation](http://en.wikipedia.org/wiki/Root-mean-square_deviation_(bioinformatics)) and [Kabsch algoritm](http://en.wikipedia.org/wiki/Kabsch_algorithm).

## Usage

    python calculate_rmsd.py molecule1.xyz molecule2.xyz

## Example

Same molecule translated in space

    python calculate_rmsd.py ethane.xyz ethane_translate.xyz

Same molecule two minimization

    python calculate_rmsd.py ethane.xyz ethane_move.xyz

Same structure, different bondlength

    python calculate_rmsd.py ethane.xyz ethane_far.xyz

## Credit and Copyright

Jimmy Charnley Kromann and Lars Bratholm


