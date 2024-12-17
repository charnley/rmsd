import matplotlib
import numpy as np
from matplotlib import patheffects
from matplotlib import pyplot as plt
from matplotlib.patches import FancyArrowPatch

import rmsd as rmsdlib

outline = patheffects.withStroke(linewidth=2, foreground="w")


def do_dot(ax, coord1, coord2, size=24, **kwargs):
    ax.plot([coord1], [coord2], "o", markersize=size, linewidth=3, color="k", **kwargs)
    return


def do_arrow(ax, pos_start, pos_end, rad=0.0, color="k"):
    edge_width = 2.0
    arrowstyle = "fancy,head_length={},head_width={},tail_width={}".format(
        2 * edge_width, 3 * edge_width, edge_width
    )
    arrow = FancyArrowPatch(
        posA=pos_start,
        posB=pos_end,
        arrowstyle=arrowstyle,
        connectionstyle=f"arc3,rad=-{rad:f}",
        color=color,
    )
    ax.add_artist(arrow)


def plot_coord(ax, atoms, coords, color="k", show_hydrogens=False):

    offset = 0.25

    for idx, (atom, coord) in enumerate(zip(atoms, coords)):

        _atom = rmsdlib.str_atom(atom)

        if not show_hydrogens and _atom == "H":
            continue

        ax.plot([coord[0]], [coord[1]], "o", markersize=20, linewidth=0, color=color)

        if _atom != "C":
            ax.text(
                coord[0],
                coord[1],
                _atom,
                verticalalignment="center",
                horizontalalignment="center",
                color="w",
                fontsize=12,
                fontweight="bold",
            )

        ax.text(
            coord[0] + offset,
            coord[1] - offset,
            idx,
            verticalalignment="center",
            horizontalalignment="center",
            color="k",
            fontsize=8,
            fontweight="bold",
            path_effects=[outline],
        )

    return


def set_global_style():
    font = {"weight": "bold", "size": 18}
    matplotlib.rc("font", **font)
    matplotlib.rc("axes", labelweight="bold")

    # Thicker spines
    matplotlib.rcParams["axes.linewidth"] = 2

    matplotlib.rcParams["xtick.major.width"] = 2
    matplotlib.rcParams["ytick.major.width"] = 2


def get_plot(n_ax=1, size=8):
    """Get a jupyter-sized plot"""
    fig, axs = plt.subplots(1, n_ax, sharey=True, sharex=True, figsize=(size * n_ax, size))

    fig.tight_layout()
    fig.subplots_adjust(hspace=0, wspace=0)
    return fig, axs


def fix_borders(ax, visibles=[False, False, True, True], fix_bounds=True):
    """Make border pretty"""

    directions = ["top", "right", "bottom", "left"]

    spines = ax.spines.items()
    spines = dict(spines)

    xticks = ax.get_xticks()
    yticks = ax.get_yticks()
    min_x, max_x = ax.get_xlim()
    min_y, max_y = ax.get_ylim()

    # Correct to the actual ticks

    (x_idxs,) = np.where((xticks > min_x) & (xticks < max_x))
    (y_idxs,) = np.where((yticks > min_y) & (yticks < max_y))
    xticks = xticks[x_idxs]
    yticks = yticks[y_idxs]

    min_x = np.min(xticks)
    max_x = np.max(xticks)

    min_y = np.min(yticks)
    max_y = np.max(yticks)

    # TODO Better ax.set_xlim()

    for direction, visible in zip(directions, visibles):

        spine = spines[direction]
        spine.set_visible(visible)

        if not visible:
            continue

        if not fix_bounds:
            continue

        if direction == "left" or direction == "right":
            spine.set_bounds(min_y, max_y)

        else:
            spine.set_bounds(min_x, max_x)
