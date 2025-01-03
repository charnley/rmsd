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


def plot_inertia(ax, pos, atoms, coord):

    center = rmsdlib.get_cm(atoms, coord)
    coord = coord - center

    pos = center[:2]

    inertia = rmsdlib.get_inertia_tensor(atoms, coord)
    eigval, eigvec = np.linalg.eig(inertia)

    eigvec = eigvec.T
    eigvec = eigvec[np.argsort(eigval)]
    # eigvec = eigvec.T

    eigvec *= -1

    print(pos)
    print(eigvec)

    arrow_options = dict(
        zorder=10,
        # path_effects=[outline],
        head_width=0.1,
        head_length=0.1,
        fc="k",
    )

    arrow1 = ax.arrow(*pos, *eigvec[0, :2] * 0.8, **arrow_options)
    arrow2 = ax.arrow(*pos, *eigvec[1, :2] * 0.8, **arrow_options)

    arrow1.set_path_effects([outline])
    arrow2.set_path_effects([outline])

    arrow1 = ax.arrow(*pos, *eigvec[0, :2] * 0.8, **arrow_options)
    arrow2 = ax.arrow(*pos, *eigvec[1, :2] * 0.8, **arrow_options)

    ax.text(
        pos[0] + eigvec[0, 0] / 2.5,
        pos[1] + eigvec[0, 1] / 2.5,
        "",
        verticalalignment="center",
        horizontalalignment="center",
        color="k",
        fontsize=12,
        fontweight="bold",
        path_effects=[outline],
        zorder=20,
    )

    ax.text(
        pos[0] + eigvec[1, 0] / 2.5,
        pos[1] + eigvec[1, 1] / 2.5,
        "",
        verticalalignment="center",
        horizontalalignment="center",
        color="k",
        fontsize=12,
        fontweight="bold",
        path_effects=[outline],
        zorder=20,
    )


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
