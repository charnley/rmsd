import matplotlib
import numpy as np
from matplotlib import patheffects
from matplotlib import pyplot as plt
from qmllib.representations import generate_fchl19  # type: ignore

import rmsd as rmsdlib

outline = patheffects.withStroke(linewidth=5, foreground="w")

FIGURE_SIZE = 4
MARKER_SIZE = 400
REP_SIZE = 0.1
REP_HEIGHT = 0.15


def plot_molecule(ax, atoms, coords, hatch="/////"):

    X = coords[:, 0]
    Y = coords[:, 1]

    ax.plot(X, Y, "-", color="k", linewidth=1)

    ax.scatter(
        X,
        Y,
        s=MARKER_SIZE,
        hatch=hatch,
        facecolor="#fff",
        edgecolor="#000",
        zorder=10,
        path_effects=[outline],
    )

    for atom, coord in zip(atoms, coords):

        if isinstance(atom, int):
            atom = rmsdlib.str_atom(atom)

        ax.text(
            *coord[:2],
            f"{atom}",
            verticalalignment="center",
            horizontalalignment="center",
            color="k",
            fontsize=12,
            fontweight="bold",
            path_effects=[outline],
            zorder=20,
        )

    return


def plot_representation(ax, atoms, coord):

    # qml vectors
    parameters = {
        "elements": np.unique(atoms),
        "pad": len(atoms),
        "rcut": 20,
        "acut": 20,
    }

    vecs = generate_fchl19(atoms, coord, **parameters)

    s = np.sum(vecs, axis=0)
    (non_zero_indicies,) = np.where(s > 10**-1)
    vecs = vecs[:, non_zero_indicies]

    offset = 0.15

    for idx, _ in enumerate(atoms):
        c = coord[idx][:2]
        c += offset
        x, y = ax.transLimits.transform(c)

        ins3 = ax.inset_axes(
            [x, y, REP_SIZE, REP_HEIGHT],
        )
        ins3.imshow(np.expand_dims(vecs[idx], axis=1), cmap="gray")
        ins3.xaxis.set_ticks([])
        ins3.yaxis.set_ticks([])

    return


# def do_dot(ax, coord1, coord2, size=24, **kwargs):
#     ax.plot([coord1], [coord2], "o", markersize=size, linewidth=3, color="k", **kwargs)
#     return


# def do_arrow(ax, pos_start, pos_end, rad=0.0, color="k"):
#     edge_width = 2.0
#     arrowstyle = "fancy,head_length={},head_width={},tail_width={}".format(
#         2 * edge_width, 3 * edge_width, edge_width
#     )
#     arrow = FancyArrowPatch(
#         posA=pos_start,
#         posB=pos_end,
#         arrowstyle=arrowstyle,
#         connectionstyle=f"arc3,rad=-{rad:f}",
#         color=color,
#     )
#     ax.add_artist(arrow)


def set_axis_default(ax, lim=2.0, use_grid=True):

    ax.set_box_aspect(1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(use_grid)

    ax.tick_params(axis="both", which="major", labelsize=0)

    for tick in ax.xaxis.get_major_ticks():

        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False)

        tick.label1.set_visible(False)
        tick.label2.set_visible(False)

    for tick in ax.yaxis.get_major_ticks():

        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False)

        tick.label1.set_visible(False)
        tick.label2.set_visible(False)


def set_global_style():
    font = {"weight": "bold", "size": 18, "family": "serif"}
    matplotlib.rc("font", **font)
    matplotlib.rc("axes", labelweight="bold")

    # Thicker spines
    # matplotlib.rcParams["axes.linewidth"] = 1

    # matplotlib.rcParams["xtick.major.width"] = 2
    # matplotlib.rcParams["ytick.major.width"] = 2


def get_plot(n_ax=1, size=FIGURE_SIZE):
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
