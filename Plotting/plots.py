import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.cm import get_cmap
from cycler import cycler
import h5py
import scipy.io
from celluloid import Camera
import torch

import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, TwoSlopeNorm


# print(plt.style.available)


# import plotly

# mpl.style.use("fast")

# mpl.rcParams["axes.prop_cycle"] = cycler(
# color=plt.style.library["ggplot"]["axes.prop_cycle"].by_key()["color"])

plt.style.use(["tableau-colorblind10", "seaborn-paper"])
mpl.rcParams["axes.prop_cycle"] = cycler(
    color=["#F05039", "#E57A77", "#EEBAB4", "#1F449C", "#3D65A5", "#7CA1CC", "#A8B6CC"]
) + cycler(
    linestyle=["-", "--", "-.", ":", "-", "--", "-."]
)  # From: https://www.datylon.com/blog/data-visualization-for-colorblind-readers  and https://ranocha.de/blog/colors/ , respectively
# cycler(color=plt.style.library["tab10"]["axes.prop_cycle"].by_key()["color"])
# ggplot seaborn-colorblind
w = 1
mpl.rcParams["axes.linewidth"] = w
mpl.rcParams["xtick.major.width"] = w
mpl.rcParams["xtick.minor.width"] = w
mpl.rcParams["ytick.major.width"] = w
mpl.rcParams["ytick.minor.width"] = w

mpl.rcParams["lines.markersize"] = 6
mpl.rcParams["lines.linewidth"] = 3
mpl.rcParams["font.size"] = 12
mpl.rcParams["legend.fontsize"] = 14
mpl.rcParams["figure.titlesize"] = 20
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["axes.labelsize"] = 12
mpl.rcParams["figure.figsize"] = (8, 6)
mpl.rcParams["figure.constrained_layout.use"] = True
mpl.rcParams["axes.formatter.use_mathtext"] = True
# mpl.rcParams['text.usetex'] = True # Use Latex
"""
mpl.rcParams['legend.loc'] = "upper_right" # Suggestion
mpl.rcParams['']
mpl.rcParams['']
"""


# Borrow some code to create colormap from color palette https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72 and https://www.delftstack.com/howto/matplotlib/custom-colormap-using-python-matplotlib/#use-rgba-values-to-create-custom-listed-colormap-in-python


def hex_to_rgb(value):
    """
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values"""
    value = value.strip("#")  # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i : i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    """
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values"""
    return [v / 256 for v in value]


def get_continuous_cmap(hex_list, float_list=None):
    """creates and returns a color map that can be used in heat map figures.
    If float_list is not provided, colour map graduates linearly between each color in hex_list.
    If float_list is provided, each color in hex_list is mapped to the respective location in float_list.

    Parameters
    ----------
    hex_list: list of hex code strings
    float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.

    Returns
    ----------
    colour map"""
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0, 1, len(rgb_list)))

    cdict = dict()
    for num, col in enumerate(["red", "green", "blue"]):
        col_list = [
            [float_list[i], rgb_list[i][num], rgb_list[i][num]]
            for i in range(len(float_list))
        ]
        cdict[col] = col_list
    cmp = LinearSegmentedColormap("XRDCT_palette_cmp", segmentdata=cdict, N=256)
    return cmp


# Diverging by using divnorm, Use TwoSlopeNorm to define min, center and max of data. Add this together with the colormap to the contourf plot.

XRDCT_palette_cmp = get_continuous_cmap(
    ["#F05039", "#E57A77", "#EEBAB4", "#1F449C", "#3D65A5", "#7CA1CC", "#A8B6CC"]
    # Might be a good choice. Need consultants # ["#F05039", "#E57A77", "#EEBAB4", "#1F449C", "#3D65A5", "#7CA1CC", "#A8B6CC"]
    # # [     "#E57A77","#EEBAB4","#7CA1CC",]
)


def plot_visualisation_grads(
    mat_paths: list, keys: list, titles: list, bins=200, range=(-4e-5, -1e-9)
) -> None:
    """
    Loads a .mat-file as a python dict.
    Extracts data, reshapes, and executes visualisation in the form of a 1D histogram

    """

    if len(mat_paths) == 4:  # Some piss about figure format
        rows = 2
        cols = 2
    else:
        rows = len(mat_paths)
        cols = 1

    fig, axs = plt.subplots(
        rows, cols, sharex=True, sharey=True, figsize=(4 * cols, 4 * rows)
    )
    cycler = plt.rcParams["axes.prop_cycle"].by_key()[
        "color"
    ]  # Use this to retrieve cycle
    for i, (ax, path) in enumerate(
        zip(np.reshape(axs, (rows * cols,)), mat_paths)
    ):  # np.reshape(axs, (len(axs[0])*len(axs[1]), )

        mat_dict = scipy.io.loadmat(path)
        data = np.reshape(
            mat_dict[keys[i]], -1
        )  # * 1e5 # Scaled, though a bit hard-coded.

        ax.hist(np.abs(data), bins=bins, range=range, density=False, alpha=0.69)
        ax.ticklabel_format(style="sci", scilimits=(0, 0))
        ax.set_ylabel("Count")
        ax.set_title(titles[i])

        t = ax.yaxis.get_offset_text()
        t.set_x(-0.2)

    ax.set_xlabel("Gradient value")
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    plt.show()
    return


def plot_comparison_grads(
    mat_paths: list, keys: list, labels: list, bins=200, range=(-4e-5, -1e-9)
) -> None:
    """
    Plots gradient fingerprints in same figure for comparison purpouses
    """
    colors = ["#F05039", "#F05039", "#1F449C", "#A8B6CC"]
    fig, ax = plt.subplots()
    for i, path in enumerate(mat_paths[1:]):
        i += 1

        mat_dict = scipy.io.loadmat(path)
        data = np.reshape(mat_dict[keys[i]], -1)

        ax.hist(
            data,
            bins=bins,
            range=range,
            density=False,
            alpha=0.69,
            color=colors[i],
            label=labels[i],
        )

    ax.ticklabel_format(style="sci", scilimits=(0, 0))
    ax.set_ylabel("Count")
    t = ax.yaxis.get_offset_text()
    t.set_x(-0.2)
    ax.set_xlabel("Gradient value")
    ax.legend()
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))
    plt.show()
    return


def GD_cost_func(x, y):
    return (
        torch.pi / 2 * (torch.sin(x**2 + y**2) + torch.pi) / (x**2 + y**2 + 1)
        + 0.69 * (x**2 + y**2)
        - torch.pi / 2 * torch.exp(-((x - 1) ** 2 + (y - 1) ** 2) / (4))
    )


def plot_3D_GD(func, x, y, first=True, fig=None, ax=None):

    if first:
        fig = plt.figure()
        ax = plt.axes(projection="3d")

    else:
        fig = fig
        ax = ax

    X, Y = np.meshgrid(x, y)
    X, Y = torch.tensor(X), torch.tensor(Y)
    Z = func(X, Y)

    ax.plot_surface(X, Y, Z, alpha=0.8, cmap="Reds")  # [u"#F05039"] )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=23, azim=0)

    return fig, ax


def animate_3D_GD(func, x, y, init_pos, frames=100, step=0.1):

    fig, ax = plot_3D_GD(func, x, y)
    camera = Camera(fig)
    traj = np.zeros((frames, 2))
    traj[0] = init_pos
    ax.plot(
        traj[0, 0],
        traj[0, -1],
        func(torch.tensor(traj[0, 0]), torch.tensor(traj[0, -1])),
        "o",
        c="#1F449C",
        markersize=14,
    )

    camera.snap()

    for epoch in range(frames - 1):

        x1, x2 = traj[epoch]
        input = torch.tensor(np.array([x1, x2]), requires_grad=True)
        print(input)
        err = func(input[0], input[-1])
        err.backward()
        grad = input.grad

        traj[epoch + 1] = traj[epoch] - step * grad.numpy()
        fig, ax = plot_3D_GD(func, x, y, first=False, fig=fig, ax=ax)
        ax.plot(
            traj[epoch + 1, 0],
            traj[epoch + 1, -1],
            func(torch.tensor(traj[epoch + 1, 0]), torch.tensor(traj[epoch + 1, -1])),
            "o",
            c="#1F449C",
            markersize=14,
        )
        ax.plot(traj[: epoch + 2, 0], traj[: epoch + 2, -1], 0, c="#1F449C")

        camera.snap()

    animation = camera.animate(interval=50, repeat=True)
    animation.save(f"GD_test.mp4")
    return


# plt.rc("text", usetex=True)


def plot_SH(ax, l=0, m=0):
    """
    Much code copied from https://scipython.com/blog/visualizing-the-real-forms-of-the-spherical-harmonics/
    """

    # Grids of polar and azimuthal angles
    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    # Create a 2-D meshgrid of (theta, phi) angles.
    theta, phi = np.meshgrid(theta, phi)
    # Calculate the Cartesian coordinates of each point in the mesh.
    xyz = np.array(
        [np.sin(theta) * np.sin(phi), np.sin(theta) * np.cos(phi), np.cos(theta)]
    )

    Y = sph_harm(0, l, phi, theta)
    print(Y)
    Yx, Yy, Yz = np.abs(Y) * xyz

    cmap = plt.cm.ScalarMappable(
        cmap=XRDCT_palette_cmp  # cmap=plt.get_cmap("PRGn")
    )  # Find out how to use own palette
    cmap.set_clim(-0.5, 0.5)
    divnorm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5)

    ax.plot_surface(
        Yx, Yy, Yz, facecolors=cmap.to_rgba(Y.real), norm=divnorm, rstride=2, cstride=2
    )

    # Draw a set of x, y, z axes for reference.
    ax_lim = 0.5
    ax.plot([-ax_lim, ax_lim], [0, 0], [0, 0], c="0.5", lw=1, zorder=10)
    ax.plot([0, 0], [-ax_lim, ax_lim], [0, 0], c="0.5", lw=1, zorder=10)
    ax.plot([0, 0], [0, 0], [-ax_lim, ax_lim], c="0.5", lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.
    ax.set_title(r"$Y_{{{},{}}}$".format(l, m))
    ax_lim = 0.5
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis("off")


def plot_debugging_coeffs(mat_paths: list, titles: list, normalise=False):
    rows = 2
    cols = 2
    fig, axs = plt.subplots(rows, cols, sharex=True, sharey=False)

    for i, (ax, path) in enumerate(zip(np.reshape(axs, (rows * cols,)), mat_paths)):

        mat = scipy.io.loadmat(path)
        coeffs = np.ndarray.flatten(mat["s"]["a"][0, 0][0, 0][0])

        if normalise:
            coeffs /= np.max(coeffs)

        x = np.arange(0, len(coeffs))
        ax.bar(x, coeffs)
        ax.set_title(titles[i])
        ax.set_xlabel("Coefficient index")
        ax.set_ylabel("Coefficient value")
        if normalise:
            ax.set_ylim(0, 1)
        else:
            ax.set_ylim(-5, 69)

    correct = np.ndarray.flatten(
        scipy.io.loadmat(mat_paths[0])["s"]["a"][0, 0][0, 0][0]
    )
    incorrect_AD = np.ndarray.flatten(
        scipy.io.loadmat(mat_paths[1])["s"]["a"][0, 0][0, 0][0]
    )
    incorrect_python = np.ndarray.flatten(
        scipy.io.loadmat(mat_paths[2])["s"]["a"][0, 0][0, 0][0]
    )
    ax = axs[-1, -1]
    if normalise:
        error_AD = np.abs(correct - incorrect_AD) / np.abs(correct)
        error_python = np.abs(correct - incorrect_python) / np.abs(correct)
        ax.set_ylabel("Relative Error")
        ax.set_ylim(0, 2)
    else:
        error_AD = np.abs(correct - incorrect_AD)
        error_python = np.abs(correct - incorrect_python)
        ax.set_ylabel("Absolute Error")

    ax.bar(x, error_AD, width=0.5, label="AD")
    ax.bar(x + 0.5, error_python, width=0.5, label="Python")
    ax.legend()

    ax.set_xlabel("Coefficient index")
    ax.set_title("Compared to symbolic")

    plt.show()
    return


def plot_debugging_orientations(mat_paths: list, titles: list, shape):
    rows = shape[0]
    cols = shape[1]

    fig, axs = plt.subplots(rows, cols, sharex=True, sharey=True)

    for i, (ax, path) in enumerate(zip(np.reshape(axs, (rows * cols,)), mat_paths)):

        mat = scipy.io.loadmat(path)
        theta = np.ndarray.flatten(mat["s"]["theta"][0, 0][0, 0][0])
        phi = np.ndarray.flatten(mat["s"]["phi"][0, 0][0, 0][0])
        x = np.arange(len(phi))
        ax.bar(x, theta, width=0.5, label="Theta")
        ax.bar(x + 0.5, phi, width=0.5, label="Phi")
        ax.set_title(titles[i])
        ax.set_xlabel("Voxel index")
        ax.set_ylabel("Orientation (radians)")
        ax.legend()

    plt.show()
    return


def plot_result_heatmap(tt_result_AD, tt_result_symbolic, title):

    ylabels = ["a0", "a2", "a4", "a6", "\u03B8", "\u03C6"]
    x_length = len(tt_result_AD)
    num_coeffs = 4

    Ad_a = np.vstack(
        (
            np.ndarray.flatten(tt_result_AD.a0),
            np.ndarray.flatten(tt_result_AD.a2),
            np.ndarray.flatten(tt_result_AD.a4),
            np.ndarray.flatten(tt_result_AD.a6),
        )
    )  # np.ndarray.flatten(tt_result_AD.a)
    Ad_theta = np.ndarray.flatten(tt_result_AD.theta)
    Ad_phi = np.ndarray.flatten(tt_result_AD.phi)

    # Sym_a = np.ndarray.flatten(tt_result_symbolic.a)
    Sym_a = np.vstack(
        (
            np.ndarray.flatten(tt_result_symbolic.a0),
            np.ndarray.flatten(tt_result_symbolic.a2),
            np.ndarray.flatten(tt_result_symbolic.a4),
            np.ndarray.flatten(tt_result_symbolic.a6),
        )
    )
    Sym_theta = np.ndarray.flatten(tt_result_symbolic.theta)
    Sym_phi = np.ndarray.flatten(tt_result_symbolic.phi)

    err_a = np.abs(Ad_a - Sym_a) / np.abs(Sym_a + 1e-10)
    err_theta = np.abs(Ad_theta - Sym_theta) / np.abs(Sym_theta + 1e-10)
    err_phi = np.abs(Ad_phi - Sym_phi) / np.abs(Sym_phi + 1e-10)
    print(Ad_a.shape)
    img = np.vstack(
        (err_a[0, :], err_a[1, :], err_a[2, :], err_a[3, :], err_theta, err_phi)
    )  # Bit manual, but but

    fig, ax = plt.subplots(figsize=(20, 8))
    ax.imshow(img, cmap=XRDCT_palette_cmp, vmin=0, vmax=2)
    ax.set_yticks(np.arange(len(ylabels)), labels=ylabels)

    for i in range(x_length):
        for j in range(num_coeffs):
            text = ax.text(
                i,
                j,
                f"{err_a[j, i]:.2f}",
                ha="center",
                va="center",
                color="w",
            )
        text = ax.text(
            i, j + 1, f"{err_theta[i]:.2f}", ha="center", va="center", color="w"
        )
        text = ax.text(
            i, j + 2, f"{err_phi[i]:.2f}", ha="center", va="center", color="w"
        )
    # cbar = ax.figure.colorbar(img, ax=ax, location="bottom")
    # cbar.ax.set_ylabel("Relative error", rotation=-90, va="bottom")
    ax.figure.colorbar(img, ax=ax, orientation="horizontal")
    ax.set_title(title)
    plt.show()


def plot_tt_3D_quiver(tt_result):

    return
