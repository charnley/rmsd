
from rmsd.calculate_rmsd import *
from rmsd.calculate_rmsd import __version__
from rmsd.calculate_rmsd import __doc__

__all__ = [\
    "kabsch_rmsd",
    "kabsch_rotate",
    "kabsch",
    "quaternion_rmsd",
    "quaternion_rotate",
    "centroid",
    "rmsd",
    "write_coordinates",
    "get_coordinates",
    "get_coordinates_pdb",
    "get_coordinates_xyz"]

if __name__ == "__main__":
    main()

