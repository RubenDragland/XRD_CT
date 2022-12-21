import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.cm import get_cmap
from cycler import cycler
import h5py
import scipy.io
from celluloid import Camera
import torch
from scipy.optimize import curve_fit

import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from scipy.special import sph_harm

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.colors import LogNorm


# print(plt.style.available)


# import plotly

# mpl.style.use("fast")

# mpl.rcParams["axes.prop_cycle"] = cycler(
# color=plt.style.library["ggplot"]["axes.prop_cycle"].by_key()["color"])

plt.style.use(["tableau-colorblind10", "seaborn-paper"])
mpl.rcParams["axes.prop_cycle"] = cycler(
    color=[
        "#3D65A5",
        "#E57A77",
        "#7CA1CC",
        "#F05039",
        "#1F449C",
        "#A8B6CC",
        "#EEBAB4",
    ]  # ["#F05039", "#E57A77", "#EEBAB4", "#1F449C", "#3D65A5", "#7CA1CC", "#A8B6CC"]
) + cycler(
    linestyle=["-", "--", "-.", ":", "-", "--", "-."]
)  # From: https://www.datylon.com/blog/data-visualization-for-colorblind-readers  and https://ranocha.de/blog/colors/ , respectively
# cycler(color=plt.style.library["tab10"]["axes.prop_cycle"].by_key()["color"])
# ggplot seaborn-colorblind
DEFAULT_FIGSIZE = (5.69, 3.9)
w = 1
mpl.rcParams["axes.linewidth"] = w
mpl.rcParams["xtick.major.width"] = w
mpl.rcParams["xtick.minor.width"] = w / 2
mpl.rcParams["ytick.major.width"] = w
mpl.rcParams["ytick.minor.width"] = w / 2

mpl.rcParams["lines.markersize"] = 6
mpl.rcParams["lines.linewidth"] = 3
mpl.rcParams["font.size"] = 12
mpl.rcParams["legend.fontsize"] = 14
mpl.rcParams["figure.titlesize"] = 20
mpl.rcParams["axes.titlesize"] = 16
mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["axes.labelsize"] = 12
mpl.rcParams["figure.figsize"] = DEFAULT_FIGSIZE  # (8, 6)
mpl.rcParams["figure.constrained_layout.use"] = True
# mpl.rcParams["axes.formatter.use_mathtext"] = True
# mpl.rcParams["text.usetex"] = True  # Use Latex
# mpl.rcParams.update(
#     {
#         "text.usetex": True,
#         "font.family": "serif",
#         "font.serif": ["Palatino"],
#     }
# )

# Try this to get font similar to latex
"""
mpl.rcParams['legend.loc'] = "upper_right" # Suggestion
mpl.rcParams['']
mpl.rcParams['']
"""


# Borrow some code to create colormap from color palette https://towardsdatascience.com/beautiful-custom-colormaps-with-matplotlib-5bab3d1f0e72 and https://www.delftstack.com/howto/matplotlib/custom-colormap-using-python-matplotlib/#use-rgba-values-to-create-custom-listed-colormap-in-python


def choose_formatter(incscape=True):
    if incscape:
        mpl.rcParams["svg.fonttype"] = "none"
        plt.rcParams["svg.fonttype"] = "none"
        plt.rcParams["axes.unicode_minus"] = False
        return
    else:
        mpl.rcParams["axes.formatter.use_mathtext"] = True
        mpl.rcParams["text.usetex"] = True  # Use Latex
        mpl.rcParams.update(
            {
                "text.usetex": True,
                "font.family": "serif",
                "font.serif": ["Palatino"],
            }
        )
    return


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
    ["#1F449C", "#3D65A5", "#7CA1CC", "#A8B6CC", "#EEBAB4", "#E57A77", "#F05039"]
    # Might be a good choice. Need consultants # ["#F05039", "#E57A77", "#EEBAB4", "#1F449C", "#3D65A5", "#7CA1CC", "#A8B6CC"]
    # # [     "#E57A77","#EEBAB4","#7CA1CC",]
)

# XRDCT_cyclic_cmp = get_continuous_cmap(
#     [
#         "#1F449C",
#         "#3D65A5",
#         "#7CA1CC",
#         "#A8B6CC",
#         "#EEBAB4",
#         "#E57A77",
#         "#F05039",
#         "#F05039",
#         "#E57A77",
#         "#EEBAB4",
#         "#A8B6CC",
#         "#7CA1CC",
#         "#3D65A5",
#         "#1F449C",
#     ]
# )

XRDCT_cyclic_cmp = get_continuous_cmap(
    [
        "#1F449C",
        "#A8B6CC",
        "#7CA1CC",
        "#E57A77",
        "#F05039",
        "#EEBAB4",
        "A8B6CC",
        "#1F449C",
    ]
)


############# DATA ANAL SECTION ####################


def gaussian(x, A, mu, sig):
    return (
        A
        / (np.sqrt(2 * np.pi) * sig)
        * np.exp(-np.power(x - mu, 2.0) / (2 * np.power(sig, 2.0)))
    )


def lorentzian(x, A, mu, gamma):
    return A / (np.pi * gamma) * gamma**2 / ((x - mu) ** 2 + gamma**2)


def exp_decay(x, A, b):
    return A * b**x


def find_vmin_vmax(AD, symbolic, std_c):
    """
    Determines lower and upper bounds for plots given input data.
    Gives much space on left for legend
    """
    # Additional space to legend.
    val1_lower = np.mean(AD) - std_c * np.std(AD)
    val2_lower = np.mean(symbolic) - std_c * np.std(symbolic)

    val1_upper = np.mean(AD) + 1.2 * std_c * np.std(AD)
    val2_upper = np.mean(symbolic) + 1.2 * std_c * np.std(symbolic)

    vmin = min(val1_lower, val2_lower)
    vmax = max(val1_upper, val2_upper)
    return vmin, vmax


def plot_SH_aligned_distribution(
    AD,
    symbolic,
    dummy_values,
    slice,
    attribute,
    title="SH ",
    bins=100,
    std_c=2,
    shape="square",
    save=False,
    save_name="SH_aligned_distribution_",
    size_fraction=0,
    shareaxis=True,
    morefigs=False,
    incscape=True,
    sharey=True,
    ticks=5,
):

    choose_formatter(incscape=incscape)

    # Set up figure
    if attribute == "coeffs":
        keys = ["a0", "a2", "a4", "a6"]
        titles = keys
        title += "Coefficients"
        save_name += "coeffs"
        if shape == "square":
            rows, cols = 2, 2
        else:
            rows, cols = 4, 1
    elif attribute == "angles":
        keys = ["theta", "phi"]
        titles = [r"$\theta$", r"$\varphi$"]
        rows, cols = 1, 2
        title += "Orientation"
        save_name += "angles"
        if shape == "horizontal":
            rows, cols = 1, 2
        else:
            rows, cols = 2, 1

    slice1, slice2 = slice
    ind = np.arange(slice1, slice2)
    X, Y, Z = np.meshgrid(ind, ind, ind)

    return_dict = {}

    if not size_fraction:
        cols_fraction = 1
        rows_fraction = 1
    else:
        cols_fraction = cols / size_fraction
        rows_fraction = rows / size_fraction

    size1, size2 = (
        cols_fraction * DEFAULT_FIGSIZE[0],
        rows_fraction * DEFAULT_FIGSIZE[1],
    )

    fig, axs = plt.subplots(rows, cols, figsize=(size1, size2), sharey=sharey)
    for i, (ax, key) in enumerate(zip(np.reshape(axs, -1), keys)):

        if attribute == "angles" and shareaxis:
            vmin, vmax = (
                dummy_values[i] - std_c * np.pi / 4,
                dummy_values[i] + std_c * np.pi / 4,
            )
        elif attribute == "coeffs" and morefigs:
            vmin, vmax = (
                dummy_values[i] - 1 / std_c * dummy_values[i],
                dummy_values[i] + 1.2 / std_c * dummy_values[i],
            )
        else:
            vmin, vmax = find_vmin_vmax(AD[key][X, Y, Z], symbolic[key][X, Y, Z], std_c)

        AD_vals = AD[key][X, Y, Z].flatten()
        symbolic_vals = symbolic[key][X, Y, Z].flatten()
        ax.set_title(titles[i])

        AD_ax = ax.hist(
            AD_vals,
            label="AD",
            alpha=0.69,
            bins=bins,
            range=(vmin, vmax),
        )
        SYM_ax = ax.hist(
            symbolic_vals,
            label="SYM",
            alpha=0.69,
            bins=bins,
            range=(vmin, vmax),
        )

        ax.axvline(
            x=dummy_values[i],
            color="black",
            linestyle="--",
            alpha=0.69,
            label=f"{dummy_values[i]:.4f}",
        )

        AD_popt, cov1 = curve_fit(lorentzian, AD_ax[1][:-1], AD_ax[0])
        SYM_popt, cov2 = curve_fit(lorentzian, SYM_ax[1][:-1], SYM_ax[0])

        AD_fit = lorentzian(AD_ax[1], *AD_popt)
        SYM_fit = lorentzian(SYM_ax[1], *SYM_popt)

        return_dict[key] = {
            "AD": AD_popt,
            "SYM": SYM_popt,
        }

        AD_color = AD_ax[-1].patches[0].get_facecolor()
        ax.plot(
            AD_ax[1],
            AD_fit,
            color=AD_color,
            alpha=1,
            label="AD Fit",  # $\mu$={:.2f}, $\sigma$={:.3f}".format(*AD_popt[1:]),
        )

        SYM_color = SYM_ax[-1].patches[0].get_facecolor()
        ax.plot(
            SYM_ax[1],
            SYM_fit,
            color=SYM_color,
            alpha=1,
            label="SYM Fit",  # $\mu$={:.2f}, $\sigma$={:.3f}".format(*SYM_popt[1:]),
        )
        handles, labels = ax.get_legend_handles_labels()
        sol_handle = handles[2]
        sol_label = labels[2]
        ax.legend([sol_handle], [sol_label], loc="upper right", frameon=False)
        ax.xaxis.set_major_locator(plt.MaxNLocator(ticks))
        ax.yaxis.set_major_locator(plt.MaxNLocator(ticks))

    fig.suptitle(title)
    fig.text(-0.02, 0.5, "Count", rotation="vertical", va="center")
    fig.text(0.5, -0.03, "Value", ha="center")

    bbox_handles = [handles[i] for i in [0, 1]]  # , 3, 4]]
    bbox_labels = [labels[i] for i in [0, 1]]  # , 3, 4]]

    fig.legend(
        bbox_handles,
        bbox_labels,
        loc="lower center",
        ncol=len(bbox_labels),
        bbox_to_anchor=(0.5, -0.6 / size2),
        frameon=False,
    )

    if save:
        plt.savefig(r"../Plotting/thesis_plots/" + save_name + ".svg")
    plt.show()
    return return_dict


def plot_SH_deviation_grad_scatter(
    AD_filename,
    SYM_filename,
    slice=None,
    title="Gradient Comparison",
    save=False,
    save_name="SH_diff_grads",
    dir=r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Results\thesis_res",
    incscape=True,
):
    choose_formatter(incscape=incscape)

    AD_grad = np.squeeze(scipy.io.loadmat(dir + r"/" + AD_filename)["grad_a"])
    SYM_grad = np.squeeze(scipy.io.loadmat(dir + r"/" + SYM_filename)["grad_a"])

    if slice is not None:
        slice1, slice2 = slice
        AD_grad = AD_grad[slice1:slice2, slice1:slice2, slice1:slice2]
        SYM_grad = SYM_grad[slice1:slice2, slice1:slice2, slice1:slice2]

    radius = np.zeros(SYM_grad.shape)
    center = np.array(radius.shape) // 2
    for index, v in np.ndenumerate(radius):
        radius[index] = np.sqrt(
            (index[0] - center[0]) ** 2
            + (index[1] - center[1]) ** 2
            + (index[2] - center[2]) ** 2
        )
    radius = np.ndarray.flatten(radius)

    difference = (AD_grad - SYM_grad) / np.abs(SYM_grad)

    fig, axs = plt.subplots(
        2,
        1,
        figsize=(DEFAULT_FIGSIZE[0], 2 * DEFAULT_FIGSIZE[1]),
        dpi=600,
        sharex=True,
    )

    axs[0].scatter(
        SYM_grad.flatten(),
        difference.flatten(),
        alpha=0.69,
        c=radius,
        cmap=XRDCT_palette_cmp,
    )
    axs[-1].scatter(
        SYM_grad.flatten(),
        np.abs(difference.flatten()),
        alpha=0.69,
        c=radius,
        cmap=XRDCT_palette_cmp,
    )
    axs[-1].set_yscale("log")

    axs[-1].xaxis.set_major_locator(plt.MaxNLocator(5))

    # axs[-1].set_xlabel("Gradient Value Sorted")
    # ax.set_ylabel("Relative Error")

    fig.text(-0.02, 0.5, "Relative Deviation", rotation="vertical", va="center")
    fig.text(0.5, -0.01, "Gradient Value Sorted", ha="center")

    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

    cmap = XRDCT_palette_cmp
    norm = mpl.colors.Normalize(vmin=0, vmax=radius.max())
    # cbar = fig.add_subplot(133, frameon=False)
    plt.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, label="Radius"
    )

    plt.suptitle(title)
    if save:
        fig.savefig(r"../Plotting/thesis_plots/" + save_name + ".png")
    plt.show()
    return


def plot_SH_diff_grads(
    AD_filename,
    SYM_filename,
    bins=100,
    title="Gradient Comparison",
    save=False,
    save_name="SH_diff_grads",
    dir=r"C:\Users\Bruker\OneDrive\Dokumenter\NTNU\XRD_CT\Results\thesis_res",
):

    AD_grad = np.squeeze(scipy.io.loadmat(dir + r"/" + AD_filename)["grad_a"])
    SYM_grad = np.squeeze(scipy.io.loadmat(dir + r"/" + SYM_filename)["grad_a"])

    difference = (AD_grad - SYM_grad) / np.abs(SYM_grad)

    fig, axs = plt.subplots(
        4, 1, figsize=(DEFAULT_FIGSIZE[0], 4 * DEFAULT_FIGSIZE[1]), dpi=100
    )
    axs[0].set_xlabel("Relative Error")
    axs[0].set_ylabel("Count")
    axs[0].hist(difference.flatten(), bins=bins, alpha=0.69)
    axs[0].set_yscale("log")

    axs[1].set_xlabel("Gradient Value")
    axs[1].set_ylabel("Count")
    axs[1].hist(SYM_grad.flatten(), bins=bins, alpha=0.69)
    axs[1].hist(AD_grad.flatten(), bins=bins, alpha=0.69)
    axs[1].set_yscale("log")

    axs[2].set_xlabel("Error Sorted")
    axs[2].set_ylabel("Relative Error")
    axs[2].scatter(SYM_grad.flatten(), np.abs(difference.flatten()), alpha=0.19)
    axs[2].set_yscale("log")
    # axs[2].set_xscale("lin")

    # col = axs[3].matshow((AD_grad[:,:,16]-SYM_grad[:,:,16]) /  np.abs(SYM_grad[:,:,16]), cmap=XRDCT_palette_cmp, LogNorm=() )
    # fig.colorbar(col)

    mean = np.mean(difference)
    std = np.std(difference)

    # ax.set_xlim(mean - 3 * std, mean + 3 * std)

    plt.suptitle(title)

    if save:
        fig.savefig(r"thesis_plots/" + save_name + ".svg")
    plt.show()

    return


def plot_performance_curves(
    sizes,
    AD_times,
    SYM_times,
    SYM_dict=None,
    help_coeff=1,
    help_line=None,
    type="convergence",
    AD_initially=None,
    title="Performance Curves",
    save=False,
    save_name="performance_curves",
    incscape=True,
):

    choose_formatter(incscape=incscape)

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE, dpi=100)

    if type == "convergence":
        label1 = "AD"
        label2 = "SYM"
        sizes = np.array(sizes) ** 3
    else:
        label1 = "AD laptop"
        label2 = "AD CPU/GPU"
        sizes = np.array(sizes)

    x = np.linspace(sizes[0], sizes[-1], 1000)
    ax.plot(sizes, AD_times, label=label1, marker="o")
    ax.plot(sizes, SYM_times, label=label2, marker="d")
    if help_line is not None:
        help_line_str = str(help_line)
        help_coeff_str = str(help_coeff)
        ax.plot(
            x,
            help_coeff * x ** (help_line),
            ":",
            label=f"${help_coeff_str}N^{{{help_line_str}}}$",
        )

    if AD_initially is not None:
        ax.scatter(
            AD_initially[0] ** 3,
            AD_initially[1],
            100,
            label=f"AD {AD_initially[1]}s ",
        )

    ax.set_xlabel("Voxels")
    ax.set_xticks(
        [1e2, 1e3, 1e4, 1e5], labels=["\$10^2$", "\$10^3$", "\$10^4$", "\$10^5$"]
    )

    if SYM_dict is not None:
        ax.plot(
            sizes,
            SYM_dict["SYM_ParFor"],
            label="SYM Parallel",
            marker="P",
            color="#EEBAB4",
        )
        ax.plot(
            sizes,
            SYM_dict["SYM"],
            label="SYM",
            marker="X",
            color="#A8B6CC",
            linestyle="dashed",
        )

    ax.set_ylabel("Time [s]")
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.set_title(title)
    # ax.legend(loc = "upper center")
    handles = ax.get_legend_handles_labels()[0]
    labels = ax.get_legend_handles_labels()[1]

    help_leg_handle = handles[2]
    help_leg_label = labels[2]
    help_leg = ax.legend(
        loc="upper left",
        handles=[help_leg_handle],
        labels=[help_leg_label],
        frameon=False,
    )
    rest_handles = handles[0:2] + handles[3:]
    rest_labels = labels[0:2] + labels[3:]
    plt.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, -0.45),
        ncol=2,
        frameon=False,
        handles=rest_handles,
        labels=rest_labels,
        # fancybox=True,
        # shadow=True,
    )
    plt.gca().add_artist(help_leg)

    if save:
        fig.savefig(r"../Plotting/thesis_plots/" + save_name + ".svg")
    plt.show()

    return


def plot_loss_curves(
    data_dict,
    title="Loss Curve",
    save=False,
    save_name="loss_curves",
    DPI=100,
    size_fraction=1.2,
    logcurve=False,
    yloglim=(1e-1, 1e1),
    incscape=True,
):
    """
    Plots convergence curves for the given data_dict
    """

    choose_formatter(incscape=incscape)  # RSD: Working?

    if logcurve:
        cols = 1
        rows = 2
        size1, size2 = (
            cols / size_fraction * DEFAULT_FIGSIZE[0],
            rows / size_fraction * DEFAULT_FIGSIZE[1],
        )
        fig, axs = plt.subplots(rows, cols, figsize=(size1, size2), dpi=DPI)
        axs[1].set_yscale("log")
        axs[1].set_ylim(yloglim)
        # axs[1].set_xscale("log")

    else:
        size1, size2 = DEFAULT_FIGSIZE
        fig, axs = plt.subplots(figsize=(size1, size2), dpi=DPI)
        axs = np.expand_dims(axs, axis=0)

    for ax in np.reshape(axs, -1):
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Loss")

        for i, (key, value) in enumerate(data_dict.items()):

            cycler = ax._get_lines.prop_cycler
            this_color = next(cycler)["color"]
            x = np.arange(0, len(value[0]), dtype=np.float64)
            params = curve_fit(exp_decay, x, value[0])[0]

            plot_value = (
                (value[0] / value[0][0]) if (len(data_dict.keys()) > 3) else value[0]
            )

            ax.plot(
                plot_value,
                color=this_color,
                label=f"{key} {value[-1]:.0f}s",
            )

        ax.legend(loc="upper right", frameon=False)

    plt.suptitle(title)
    if save:
        fig.savefig(r"../Plotting/thesis_plots/" + save_name + ".svg")
    plt.show()
    return


def plot_slices(
    symbolic,
    automatic,
    expsin,
    key_SH,
    key_EXPSIN,
    slice,
    shape="vertical",
    title="Reconstruction Slices",
    save=False,
    save_name="reconstruction_slices",
    incscape=False,
    size_fraction=1,
    DoA_tol=-1,
):

    choose_formatter(incscape=incscape)
    if shape == "vertical":
        cols = 1
        rows = 3
    elif shape == "horizontal":
        cols = 3
        rows = 1
    else:
        raise ValueError("shape must be either vertical or horizontal")

    slice_y, slice_x, slice_z = slice

    sym_data = np.squeeze(
        symbolic[key_SH][
            slice_y[0] : slice_y[1], slice_x[0] : slice_x[1], slice_z[0] : slice_z[1]
        ]
    )

    aut_data = np.squeeze(
        automatic[key_SH][
            slice_y[0] : slice_y[1], slice_x[0] : slice_x[1], slice_z[0] : slice_z[1]
        ]
    )

    exp_data = np.squeeze(
        expsin[key_EXPSIN][
            slice_y[0] : slice_y[1], slice_x[0] : slice_x[1], slice_z[0] : slice_z[1]
        ]
    ) * (expsin.sign if key_EXPSIN == "B" else 1)

    min_value = np.min([sym_data, aut_data, exp_data])
    max_value = np.max([sym_data, aut_data, exp_data])

    cmap = (
        XRDCT_cyclic_cmp
        if (key_SH == "theta" or key_SH == "phi")
        else XRDCT_palette_cmp
    )
    norm = (
        mpl.colors.Normalize(vmin=0, vmax=np.pi)
        if (key_SH == "theta" or key_SH == "phi")
        else mpl.colors.Normalize(vmin=min_value, vmax=max_value)
    )
    if key_SH == "DoA":
        norm = mpl.colors.Normalize(vmin=-0.5, vmax=1)

        # sym_data = np.where(
        #     np.squeeze(
        #         symbolic["a0"][
        #             slice_y[0] : slice_y[1],
        #             slice_x[0] : slice_x[1],
        #             slice_z[0] : slice_z[1],
        #         ]
        #     )
        #     > DoA_tol,
        #     sym_data,
        #     0,
        # )
        # aut_data = np.where(
        #     np.squeeze(
        #         automatic["a0"][
        #             slice_y[0] : slice_y[1],
        #             slice_x[0] : slice_x[1],
        #             slice_z[0] : slice_z[1],
        #         ]
        #     )
        #     > DoA_tol,
        #     aut_data,
        #     0,
        # )
        # exp_data = np.where(
        #     np.squeeze(
        #         expsin["A"][
        #             slice_y[0] : slice_y[1],
        #             slice_x[0] : slice_x[1],
        #             slice_z[0] : slice_z[1],
        #         ]
        #     )
        #     > DoA_tol,
        #     exp_data,
        #     0,
        # )

    if not size_fraction:
        cols_fraction = 1
        rows_fraction = 1
    else:
        cols_fraction = cols / size_fraction
        rows_fraction = rows / size_fraction

    size1, size2 = (
        cols_fraction * DEFAULT_FIGSIZE[0],
        rows_fraction * DEFAULT_FIGSIZE[1],
    )
    fig, axs = plt.subplots(rows, cols, figsize=(size1, size2), dpi=100)
    # fig.suptitle(title)
    axs[0].matshow(sym_data, cmap=cmap, norm=norm)
    axs[1].matshow(aut_data, cmap=cmap, norm=norm)
    axs[2].matshow(exp_data, cmap=cmap, norm=norm)

    axs[0].set_title("SYM")
    axs[1].set_title("AD")
    axs[2].set_title("EXPSIN")

    all_axes = axs.reshape(-1)
    for a in all_axes:
        a.set_xticks([])
        a.set_yticks([])
        a.set_aspect("equal")

    fig.subplots_adjust(bottom=0.1)
    cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.05])

    if key_SH == "phi":
        label = r"$\varphi$"
        plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            label=f"{label} [rad]",
            orientation="horizontal",
        )
    elif key_SH == "theta":
        label = r"$\theta$"
        plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            label=f"{label} [rad]",
            orientation="horizontal",
        )
    elif key_SH == "DoA":
        label = r"DoA [a.u.]"
        plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            label=f"{label}",
            orientation="horizontal",
        )
    else:
        plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            label=f"{key_SH}/{key_EXPSIN} [a.u.]",
            orientation="horizontal",
        )

    if save:
        fig.savefig(r"../Plotting/thesis_plots/" + save_name + ".svg")

    plt.show()

    return


def plot_individial_reconstructions(
    symbolic,
    automatic,
    expsin,
    slice,
    save=False,
    save_name="reconstruction_slices",
    incscape=True,
    size_fraction=1,
    arrow_length=0.1,
    intensity_threshold=2,
    view=None,
):

    choose_formatter(incscape=incscape)
    res = {"SYM": symbolic, "AD": automatic, "EXPSIN": expsin}
    slice_y, slice_x, slice_z = slice

    if not size_fraction:
        cols_fraction = 1
        rows_fraction = 1
    else:
        cols_fraction = 1 / size_fraction
        rows_fraction = 1 / size_fraction

    size1, size2 = (
        cols_fraction * DEFAULT_FIGSIZE[0],
        rows_fraction * DEFAULT_FIGSIZE[1],
    )

    cmap = XRDCT_palette_cmp
    norm = mpl.colors.Normalize(vmin=-0.5, vmax=1)
    colors = [
        elem.DoA[
            slice_y[0] : slice_y[-1],
            slice_x[0] : slice_x[-1],
            slice_z[0] : slice_z[-1],
        ]
        for elem in res.values()
    ]
    for i, key in enumerate(res.keys()):

        fig = plt.figure(figsize=(size1, size2))
        ax = fig.add_subplot(111, projection="3d")
        # ax.set_title(key)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.grid(False)
        plt.axis("off")

        if view == None:
            ax.set_xlim(0, len(symbolic.a0))
            ax.set_ylim(0, len(symbolic.a0[0]))
            ax.set_zlim(0, len(symbolic.a0[0, 0]))
            ax.text(
                len(symbolic.a0) / 2, len(symbolic.a0[0]), len(symbolic.a0[0, 0]), key
            )
        else:
            ax.set_xlim(view[0][0], view[0][1])
            ax.set_ylim(view[1][0], view[1][1])
            ax.set_zlim(view[2][0], view[2][1])
            ax.text2D(0.45, 0.7, s=key, fontsize=16, transform=ax.transAxes)

        x, y, z = np.mgrid[
            slice_y[0] : slice_y[-1], slice_x[0] : slice_x[-1], slice_z[0] : slice_z[-1]
        ]

        if key == "EXPSIN":
            check_value = res[key].A[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ]
        else:
            check_value = np.abs(
                res[key].a0[
                    slice_y[0] : slice_y[-1],
                    slice_x[0] : slice_x[-1],
                    slice_z[0] : slice_z[-1],
                ]
            )

        above_threshold = np.argwhere(check_value >= intensity_threshold)
        print(f"{key}: {above_threshold.shape[0]} voxels above threshold")
        x, y, z = (
            x[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            y[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            z[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        # RSD: Nb! Convention of theta and phi.
        dirx, diry, dirz = (
            (np.sin(res[key].theta) * np.cos(res[key].phi))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            (np.sin(res[key].phi) * np.sin(res[key].theta))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            np.cos(res[key].theta)[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
        )

        dirx, diry, dirz = (
            dirx[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            diry[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            dirz[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        q = ax.quiver(
            x,
            y,
            z,
            dirx,
            diry,
            dirz,
            cmap=cmap,
            norm=norm,
            length=arrow_length,
            pivot="middle",
            arrow_length_ratio=0,
            normalize=True,
        )
        color_floats = colors[i][
            above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]
        ].reshape(
            -1
        )  # RSD: Still not sure about this one.
        q.set_array(color_floats)

        # S = np.array([x, y, z, [dirx, diry, dirz, color_floats]])
        # scipy.io.savemat(f"{save_name}_{key}.mat", mdict={"S": S})

        if save:
            fig.savefig(r"../Plotting/thesis_plots/" + save_name + key + ".svg")
        plt.show()

    fig, cbarax = plt.subplots(figsize=(size1, size2 / 5))
    plt.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cbarax,
        label=f"Degree of Anisotropy",
        orientation="horizontal",
    )
    if save:
        fig.savefig(r"../Plotting/thesis_plots/" + save_name + "colorbar.svg")
    plt.show()


def plot_reconstruction(
    symbolic,
    automatic,
    expsin,
    slice,
    shape="horizontal",
    degree_of_anisotropy=1,
    title="Reconstruction Slices",
    save=False,
    save_name="reconstruction_slices",
    incscape=True,
    size_fraction=1,
    arrow_length=0.1,
    intensity_threshold=2,
    view=None,
):

    choose_formatter(incscape=incscape)

    res = {"SYM": symbolic, "AD": automatic, "EXPSIN": expsin}

    slice_y, slice_x, slice_z = slice

    if shape == "vertical":
        cols = 1
        rows = 3
    elif shape == "horizontal":
        cols = 3
        rows = 1
    else:
        raise ValueError("shape must be either vertical or horizontal")

    if not size_fraction:
        cols_fraction = 1
        rows_fraction = 1
    else:
        cols_fraction = cols / size_fraction
        rows_fraction = rows / size_fraction

    size1, size2 = (
        cols_fraction * DEFAULT_FIGSIZE[0],
        rows_fraction * DEFAULT_FIGSIZE[1],
    )

    if degree_of_anisotropy:
        cmap = XRDCT_palette_cmp
        norm = mpl.colors.Normalize(vmin=-0.5, vmax=1)
        colors = [
            elem.DoA[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ]
            for elem in res.values()
        ]
        # RSD: Something weird here.

    else:
        cmap = XRDCT_palette_cmp
        theta_re, phi_img = np.mgrid[0 : np.pi : 100j, 0 : np.pi : 100j]
        angle = np.angle(theta_re + 1j * phi_img)
        norm = mpl.colors.Normalize(vmin=np.min(angle), vmax=np.max(angle))
        colors = [
            np.angle(elem.theta + 1j * elem.phi)[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ]
            for elem in res.values()
        ]

    # fig, axs = plt.subplots(rows, cols, figsize=(size1, size2))
    fig = plt.figure(figsize=(size1, size2), constrained_layout=True)

    # for i, (key, ax) in enumerate(zip(["SYM", "AD", "EXPSIN"], axs.reshape(-1))):
    for i, key in enumerate(["SYM", "AD", "EXPSIN"]):
        ax = fig.add_subplot(rows, cols, i + 1, projection="3d")
        # ax.set_title(key)

        x, y, z = np.mgrid[
            slice_y[0] : slice_y[-1], slice_x[0] : slice_x[-1], slice_z[0] : slice_z[-1]
        ]

        if key == "EXPSIN":
            check_value = res[key].A[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ]
        else:
            check_value = np.abs(
                res[key].a0[
                    slice_y[0] : slice_y[-1],
                    slice_x[0] : slice_x[-1],
                    slice_z[0] : slice_z[-1],
                ]
            )

        above_threshold = np.argwhere(check_value >= intensity_threshold)
        print(f"{key}: {above_threshold.shape[0]} voxels above threshold")
        x, y, z = (
            x[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            y[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            z[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        # RSD: Nb! Convention of theta and phi.
        dirx, diry, dirz = (
            (np.sin(res[key].theta) * np.cos(res[key].phi))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            (np.sin(res[key].phi) * np.sin(res[key].theta))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            np.cos(res[key].theta)[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
        )

        dirx, diry, dirz = (
            dirx[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            diry[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            dirz[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        q = ax.quiver(
            x,
            y,
            z,
            dirx,
            diry,
            dirz,
            cmap=cmap,
            norm=norm,
            length=arrow_length,
            pivot="middle",
            arrow_length_ratio=0,
            normalize=True,
        )
        color_floats = colors[i][
            above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]
        ].reshape(
            -1
        )  # RSD: Still not sure about this one.
        q.set_array(color_floats)

        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.grid(False)
        plt.axis("off")
        # ax.set_xlabel("x")
        # ax.set_ylabel("y")
        # ax.set_zlabel("z")

        if view == None:
            ax.set_xlim(0, len(symbolic.a0))
            ax.set_ylim(0, len(symbolic.a0[0]))
            ax.set_zlim(0, len(symbolic.a0[0, 0]))
            ax.text(
                len(symbolic.a0) / 2, len(symbolic.a0[0]), len(symbolic.a0[0, 0]), key
            )
        else:
            ax.set_xlim(view[0][0], view[0][1])
            ax.set_ylim(view[1][0], view[1][1])
            ax.set_zlim(view[2][0], view[2][1])
            ax.text2D(
                0.45, 0.7, s=key, fontsize=16, transform=ax.transAxes
            )  # 0*(view[0][0] + view[0][1])/3, view[1][1] , view[2][1]

    # RSD: Path set vertices to draw a marker manually

    if degree_of_anisotropy:
        cbar_ax = fig.add_axes([0.2 / cols, -0.1 / rows, 0.6 / cols, 0.05 / rows])
        plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            label=f"Degree of Anisotropy",
            orientation="horizontal",
        )
    else:
        cbar_ax = fig.add_axes([0.33 / cols, -0.4 / rows, 0.4 / cols, 0.4 / rows])
        cbar_ax.imshow(angle, cmap=cmap, norm=norm)
        cbar_ax.set_xlabel(r"$\theta$")
        cbar_ax.set_ylabel(r"$\varphi$")
        cbar_ax.set_xticks(np.linspace(0, 100, 3))
        cbar_ax.set_xticklabels([r"$0$", r"$\pi/2$", r"$\pi$"])
        cbar_ax.set_yticks(np.linspace(0, 100, 3))
        cbar_ax.set_yticklabels([r"$0$", r"$\pi/2$", r"$\pi$"])

    fig.subplots_adjust(left=0, bottom=-0.2, right=1, top=1, wspace=-0.4, hspace=-0.4)
    # width = fig.get_figwidth()
    # height = fig.get_figheight()
    # fig.set_figwidth(width)
    # fig.set_figheight(height)

    if save:
        fig.savefig(r"../Plotting/thesis_plots/" + save_name + ".svg")

    plt.show()
    return


def plot_2D_vector_field(
    symbolic,
    automatic,
    expsin,
    slice,
    save=False,
    save_name="2D_vector_field",
    incscape=True,
    size_fraction=1,
    arrow_length=0.1,
    intensity_threshold=2,
    view=None,
    shape="horizontal",
):
    import matplotlib.path as mpath
    import matplotlib.lines as mlines
    import matplotlib.patches as mpatches
    from matplotlib.collections import PatchCollection
    from matplotlib.collections import EllipseCollection

    choose_formatter(incscape=incscape)

    if shape == "vertical":
        cols = 1
        rows = 3
    elif shape == "horizontal":
        cols = 3
        rows = 1
    else:
        raise ValueError("shape must be either vertical or horizontal")

    res = {"SYM": symbolic, "AD": automatic, "EXPSIN": expsin}
    slice_y, slice_x, slice_z = slice

    if not size_fraction:
        cols_fraction = 1
        rows_fraction = 1
    else:
        cols_fraction = cols / size_fraction
        rows_fraction = rows / size_fraction

    size1, size2 = (
        cols_fraction * DEFAULT_FIGSIZE[0],
        rows_fraction * DEFAULT_FIGSIZE[1],
    )

    cmap = XRDCT_palette_cmp
    norm = mpl.colors.Normalize(vmin=-0.5, vmax=1)
    colors = [
        elem.DoA[
            slice_y[0] : slice_y[-1],
            slice_x[0] : slice_x[-1],
            slice_z[0] : slice_z[-1],
        ]
        for elem in res.values()
    ]

    fig, axs = plt.subplots(rows, cols, figsize=(size1, size2))

    for i, (key, ax) in enumerate(zip(res.keys(), axs.flatten())):

        dirx, diry, dirz = (
            (np.sin(res[key].theta) * np.cos(res[key].phi))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            (np.sin(res[key].phi) * np.sin(res[key].theta))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            np.cos(res[key].theta)[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
        )

        y = np.arange(slice_y[0], slice_y[-1])
        x = np.arange(slice_x[0], slice_x[-1])
        angle = res[key].phi[
            slice_y[0] : slice_y[-1],
            slice_x[0] : slice_x[-1],
            slice_z[0] : slice_z[-1],
        ]  # np.arctan2(diry, dirx)

        X, Y = np.meshgrid(x, y)
        XY = np.column_stack((X.ravel(), Y.ravel()))

        ww = arrow_length - dirz * arrow_length
        # hh = arrow_length #- dirz * arrow_length
        hh = ww * (np.exp(1 + np.abs(colors[i])))
        # ww = hh * colors[i]
        aa = angle

        ec = EllipseCollection(
            ww,
            hh,
            aa,
            units="xy",
            offsets=XY,
            transOffset=ax.transData,
            cmap=cmap,
            norm=norm,
        )  # RSD: Cmap?
        ec.set_array(colors[i].ravel())
        ax.add_collection(ec)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(slice_x[0], slice_x[-1])
        ax.set_ylim(slice_y[0], slice_y[-1])
        ax.set_title(key)
        ax.set_xticklabels([])
        ax.set_yticklabels([])

        fig.subplots_adjust(bottom=0.4)
        cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.05])
        label = r"DoA [a.u.]"
        plt.colorbar(
            mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
            cax=cbar_ax,
            label=f"{label}",
            orientation="horizontal",
        )

    return


def plot_orthogonal_projections(
    symbolic,
    automatic,
    expsin,
    slice,
    degree_of_anisotropy=1,
    title="Orthogonal Projections",
    save=False,
    save_name="reconstruction_slices",
    incscape=False,
    size_fraction=1,
    arrow_length=0.1,
    intensity_threshold=2,
):
    """Plot orthogonal projections of the reconstruction."""

    choose_formatter(incscape=incscape)

    res = {"SYM": symbolic, "AD": automatic, "EXPSIN": expsin}

    slice_y, slice_x, slice_z = slice

    cols, rows = 3, 3

    if not size_fraction:
        cols_fraction = 1
        rows_fraction = 1
    else:
        cols_fraction = cols / size_fraction
        rows_fraction = rows / size_fraction

    size1, size2 = (
        cols_fraction * DEFAULT_FIGSIZE[0],
        rows_fraction * DEFAULT_FIGSIZE[1],
    )

    cmap = XRDCT_palette_cmp
    norm = mpl.colors.Normalize(vmin=-0.5, vmax=0)
    colors_list = [
        elem.DoA[
            slice_y[0] : slice_y[-1],
            slice_x[0] : slice_x[-1],
            slice_z[0] : slice_z[-1],
        ]
        for elem in res.values()
    ]
    colors = {"SYM": colors_list[0], "AD": colors_list[1], "EXPSIN": colors_list[2]}

    views = [
        ("SXY", (90, -90, 0)),
        ("AXY", (90, -90, 0)),
        ("EXY", (90, -90, 0)),
        ("SXZ", (0, -90, 0)),
        ("AXZ", (0, -90, 0)),
        ("EXZ", (0, -90, 0)),
        ("SYZ", (0, 0, 0)),
        ("AYZ", (0, 0, 0)),
        ("EYZ", (0, 0, 0)),
    ]
    layout = [
        ["SXY", "AXY", "EXY"],
        ["SXZ", "AXZ", "EXZ"],
        ["SYZ", "AYZ", "EYZ"],
    ]
    fig, axs = plt.subplot_mosaic(
        layout,
        figsize=(size1, size2),
        subplot_kw=dict(projection="3d"),
        constrained_layout=True,
    )

    for i, (plane, angles) in enumerate(views):
        ax = axs[plane]
        ax.set_proj_type("ortho")
        ax.view_init(elev=angles[0], azim=angles[1])  # -, roll=angles[2])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        lim_scalar = 8
        ax.set_xlim(slice_x[0] + lim_scalar, slice_x[-1] - lim_scalar)
        ax.set_ylim(slice_y[0] + lim_scalar, slice_y[-1] - lim_scalar)
        ax.set_zlim(slice_z[0] + lim_scalar, slice_z[-1] - lim_scalar)
        # ax.set_title(plane)
        ax.text2D(
            0.45, 0.65, plane, transform=ax.transAxes, fontsize=14 / size_fraction
        )
        # ax.set_aspect("equal")
        ax.axis("off")
        # ax.axis("tight")
        # ax.axis("image")

        x, y, z = np.mgrid[
            slice_y[0] : slice_y[-1], slice_x[0] : slice_x[-1], slice_z[0] : slice_z[-1]
        ]
        if plane[0] == "S":
            key = "SYM"
        elif plane[0] == "A":
            key = "AD"
        elif plane[0] == "E":
            key = "EXPSIN"

        if key == "EXPSIN":
            check_value = res[key].A[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ]
        else:
            check_value = np.abs(
                res[key].a0[
                    slice_y[0] : slice_y[-1],
                    slice_x[0] : slice_x[-1],
                    slice_z[0] : slice_z[-1],
                ]
            )

        above_threshold = np.argwhere(check_value >= intensity_threshold)
        print(f"{key}: {above_threshold.shape[0]} voxels above threshold")
        x, y, z = (
            x[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            y[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            z[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        # RSD: Nb! Convention of theta and phi.
        dirx, diry, dirz = (
            (np.sin(res[key].theta) * np.cos(res[key].phi))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            (np.sin(res[key].phi) * np.sin(res[key].theta))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            np.cos(res[key].theta)[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
        )

        dirx, diry, dirz = (
            dirx[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            diry[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            dirz[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        q = ax.quiver(
            x,
            y,
            z,
            dirx,
            diry,
            dirz,
            cmap=cmap,
            norm=norm,
            length=arrow_length,
            pivot="middle",
            arrow_length_ratio=0,
            normalize=True,
        )
        color_floats = colors[key][
            above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]
        ].reshape(
            -1
        )  # RSD: Still not sure about this one.
        q.set_array(color_floats)

    # fig.subplots_adjust(bottom=0.4)
    # cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.05])
    # label = r"DoA [a.u.]"
    # plt.colorbar(
    #     mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    #     cax=cbar_ax,
    #     label=f"{label}",
    #     orientation="horizontal",
    # )
    fig.subplots_adjust(wspace=0, hspace=0)
    if save:
        fig.savefig(
            r"../Plotting/thesis_plots/" + save_name + ".pdf",
            bbox_inches="tight",
            pad_inches=0,
        )

    plt.show()

    choose_formatter()

    fig, cbarax = plt.subplots(figsize=(DEFAULT_FIGSIZE[0], DEFAULT_FIGSIZE[1] / 5))
    plt.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cbarax,
        label=f"Degree of Anisotropy",
        orientation="horizontal",
    )
    if save:
        fig.savefig(
            r"../Plotting/thesis_plots/" + save_name + "colorbar_equatorial.svg"
        )
    plt.show()

    return


def plot_3D_ellipsoids(
    symbolic,
    automatic,
    expsin,
    slice,
    degree_of_anisotropy=1,
    title="Orthogonal Projections",
    save=False,
    save_name="reconstruction_slices",
    incscape=False,
    size_fraction=1,
    arrow_length=0.1,
    intensity_threshold=2,
):
    # from mpl_toolkits.mplot3d.art3d import Ellipsoid
    from mpl_toolkits.mplot3d.art3d import EllipsoidCollection

    choose_formatter(incscape=incscape)

    res = {"SYM": symbolic, "AD": automatic, "EXPSIN": expsin}

    slice_y, slice_x, slice_z = slice

    cols, rows = 3, 1

    if not size_fraction:
        cols_fraction = 1
        rows_fraction = 1
    else:
        cols_fraction = cols / size_fraction
        rows_fraction = rows / size_fraction

    size1, size2 = (
        cols_fraction * DEFAULT_FIGSIZE[0],
        rows_fraction * DEFAULT_FIGSIZE[1],
    )

    cmap = XRDCT_palette_cmp
    norm = mpl.colors.Normalize(vmin=-0.5, vmax=0)
    colors_list = [
        elem.DoA[
            slice_y[0] : slice_y[-1],
            slice_x[0] : slice_x[-1],
            slice_z[0] : slice_z[-1],
        ]
        for elem in res.values()
    ]
    colors = {"SYM": colors_list[0], "AD": colors_list[1], "EXPSIN": colors_list[2]}

    # fig, axs = plt.subplots(1, 3, figsize=(size1, size2), subplot_kw={"projection": "3d"})
    fig = plt.figure(figsize=(size1, size2), constrained_layout=True)
    for i, key in enumerate(res.keys()):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        ax.set_title(key, fontsize=20)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        lim_scalar = 8
        ax.set_xlim(slice_x[0] + lim_scalar, slice_x[-1] - lim_scalar)
        ax.set_ylim(slice_y[0] + lim_scalar, slice_y[-1] - lim_scalar)
        ax.set_zlim(slice_z[0] + lim_scalar, slice_z[-1] - lim_scalar)
        ax.set_view(azim=90, elev=-90)

        if key == "EXPSIN":
            check_value = res[key].A[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ]
        else:
            check_value = np.abs(
                res[key].a0[
                    slice_y[0] : slice_y[-1],
                    slice_x[0] : slice_x[-1],
                    slice_z[0] : slice_z[-1],
                ]
            )

        above_threshold = np.argwhere(colors[key] > check_value)

        x, y, z = np.mgrid[
            slice_y[0] : slice_y[-1], slice_x[0] : slice_x[-1], slice_z[0] : slice_z[-1]
        ]
        x, y, z = (
            x[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            y[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            z[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        # RSD: Nb! Convention of theta and phi.
        dirx, diry, dirz = (
            (np.sin(res[key].theta) * np.cos(res[key].phi))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            (np.sin(res[key].phi) * np.sin(res[key].theta))[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
            np.cos(res[key].theta)[
                slice_y[0] : slice_y[-1],
                slice_x[0] : slice_x[-1],
                slice_z[0] : slice_z[-1],
            ],
        )

        dirx, diry, dirz = (
            dirx[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            diry[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
            dirz[above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]],
        )

        color_floats = colors[key][
            above_threshold[:, 0], above_threshold[:, 1], above_threshold[:, 2]
        ]  # .reshape(
        #     -1
        # )

        ellipsoids = EllipsoidCollection(
            (dirx, diry, dirz), (x, y, z), color=color_floats
        )

        ax.add_collection3d(ellipsoids)

    return


###########################################

# Additional plots

# def plot_DoA_EXPSIN_example(): RSD: Maybe baby.


def plot_ck_expected_scattering():
    N = 16
    theta = np.linspace(0, 2 * np.pi, 360)
    radii = np.ones(360)
    width = 2 * np.pi / 16
    ax = plt.subplot(111, projection="polar")
    cmap = XRDCT_palette_cmp
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    sectors = np.arange(0, 2 * np.pi, 2 * np.pi / N)  # * 180 / np.pi
    radii = np.ones(len(sectors))

    # sector_highlighted = np.zeros(500)
    # sector_highlighted[sectors[2] : sectors[3]] = 1
    # sector_highlighted[sectors[6] : sectors[7]] = 1

    sector_highlighted = np.ones(N) * 0.5
    sector_highlighted[3:4] = 1
    sector_highlighted[11:12] = 1
    colors = cmap(sector_highlighted)
    ax.bar(sectors, radii, width=width, alpha=0.69, color=colors)
    ax.set_xticks([])
    ax.set_yticks([])
    # ax.set_yticks(sectors)
    # ax.set_xticklabels(["", "", "", "", "", "", "", "", "", "", "", "", "", "", "", ""])
    plt.savefig(r"../Figures/ck_expected_scattering.svg")

    return


def plot_sinogram(image, save=True, save_name="sinogram.svg"):
    from skimage.transform import radon, rescale

    choose_formatter(True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=False)

    ax1.set_title("Image")
    ax1.imshow(image, cmap=plt.cm.Greys_r)

    theta = np.linspace(0.0, 180.0, max(image.shape), endpoint=False)
    sinogram = radon(image, theta=theta)
    dx, dy = 0.5 * 180.0 / max(image.shape), 0.5 / sinogram.shape[0]
    ax2.set_title("Sinogram")
    ax2.set_xlabel("Projection angle [deg]")
    ax2.set_ylabel("Projection position [pixels]")
    ax2.imshow(
        sinogram,
        cmap=plt.cm.Greys_r,
        extent=(-dx, 180.0 + dx, -dy, sinogram.shape[0] + dy),
        aspect="auto",
    )

    fig.tight_layout()
    if save:
        plt.savefig(r"./thesis_plots/" + save_name)

    plt.show()


###########################################


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
    alt_cmap = XRDCT_palette_cmp

    ax.plot_surface(X, Y, Z, alpha=0.8, cmap=alt_cmap)  # cmap="Reds")  # [u"#F05039"] )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=23, azim=0)

    return fig, ax


def animate_3D_GD(
    func, x, y, init_pos, frames=100, step=0.1, animate=False, filename="GD_3D"
):

    fig, ax = plot_3D_GD(func, x, y)
    if animate:
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
    if animate:
        camera.snap()

    for epoch in range(frames - 1):

        x1, x2 = traj[epoch]
        input = torch.tensor(np.array([x1, x2]), requires_grad=True)
        err = func(input[0], input[-1])
        err.backward()
        grad = input.grad

        traj[epoch + 1] = traj[epoch] - step * grad.numpy()

        if animate:
            fig, ax = plot_3D_GD(func, x, y, first=False, fig=fig, ax=ax)
            ax.plot(
                traj[epoch + 1, 0],
                traj[epoch + 1, -1],
                func(
                    torch.tensor(traj[epoch + 1, 0]), torch.tensor(traj[epoch + 1, -1])
                ),
                "o",
                c="#1F449C",
                markersize=14,
            )
            ax.plot(traj[: epoch + 2, 0], traj[: epoch + 2, -1], 0, c="#1F449C")

            camera.snap()

    if animate:
        animation = camera.animate(interval=50, repeat=True)
        animation.save(f"{filename}.mp4")
    else:

        fig, ax = plot_3D_GD(func, x, y, first=False, fig=fig, ax=ax)
        ax.plot(
            traj[-1, 0],
            traj[-1, -1],
            func(torch.tensor(traj[-1, 0]), torch.tensor(traj[-1, -1])),
            "o",
            c="#A8B6CC",
            markersize=14,
        )
        ax.plot(traj[:, 0], traj[:, -1], 0, c="#A8B6CC")
        fig.savefig(f"{filename}.svg")

    return


def animate_3D_CGD(
    func,
    x,
    y,
    init_pos,
    frames=100,
    animate=False,
    interval=50,
    max_step=2,
    linesearches=1000,
    filename="CGD_3D",
):

    fig, ax = plot_3D_GD(func, x, y)
    if animate:
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
    if animate:
        camera.snap()

    alphas = np.linspace(1e-6, max_step, linesearches)
    x1, x2 = traj[0]
    for epoch in range(frames - 1):

        input = torch.tensor(np.array([x1, x2]), requires_grad=True)
        err = func(input[0], input[-1])
        err.backward()
        grad = input.grad.numpy()

        if epoch == 0:
            d = -grad
        else:
            beta = np.max(
                np.dot(grad.T, grad - grad_old) / np.dot(grad_old.T, grad_old), 0
            )
            # if np.any(beta) < 0:
            #     beta = 0
            d = -grad + beta * d_old

        grad_old = grad
        func_eval_old = 1000
        next_alpha = alphas[0]
        for i, a in enumerate(alphas):
            input = torch.tensor(
                np.array([x1 + a * d[0], x2 + a * d[-1]]), requires_grad=False
            )
            func_eval = func(input[0], input[-1])

            if func_eval <= func_eval_old:
                func_eval_old = func_eval
                next_alpha = a

        if next_alpha == alphas[0]:
            print(f"Converged in {epoch +1} steps")
            break

        traj[epoch + 1] = traj[epoch] + next_alpha * d
        d_old = d
        x1, x2 = traj[epoch + 1]

        if animate:
            fig, ax = plot_3D_GD(func, x, y, first=False, fig=fig, ax=ax)
            ax.plot(
                traj[epoch + 1, 0],
                traj[epoch + 1, -1],
                func(
                    torch.tensor(traj[epoch + 1, 0]), torch.tensor(traj[epoch + 1, -1])
                ),
                "o",
                c="#1F449C",
                markersize=14,
            )
            ax.plot(traj[: epoch + 2, 0], traj[: epoch + 2, -1], 0, c="#1F449C")

            camera.snap()
    if animate:
        animation = camera.animate(interval=interval, repeat=True)
        animation.save(f"{filename}.mp4")
    else:

        fig, ax = plot_3D_GD(func, x, y, first=False, fig=fig, ax=ax)
        ax.plot(
            traj[-1, 0],
            traj[-1, -1],
            func(torch.tensor(traj[-1, 0]), torch.tensor(traj[-1, -1])),
            "o",
            c="#A8B6CC",
            markersize=14,
        )
        ax.plot(traj[:, 0], traj[:, -1], 0, "X", c="#A8B6CC", markersize=10)
        fig.savefig(f"{filename}.svg")


# plt.rc("text", usetex=True)


def plot_SH(ax, l=0, m=0, sphere=False, zenith=(0, 0), intensity=1):
    """
    Much code copied from https://scipython.com/blog/visualizing-the-real-forms-of-the-spherical-harmonics/
    """

    # Grids of polar and azimuthal angles
    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    t0 = 0  # zenith[0]
    p0 = 0  # zenith[1]
    # Create a 2-D meshgrid of (theta, phi) angles.
    theta, phi = np.meshgrid(theta, phi)
    # Calculate the Cartesian coordinates of each point in the mesh.
    xyz = np.array(
        [
            np.sin(theta - t0) * np.sin(phi - p0),
            np.sin(theta - t0) * np.cos(phi - p0),
            np.cos(theta - t0),
        ]
    )

    Y = sph_harm(0, l, phi, theta) * intensity
    # print(Y)
    Yx, Yy, Yz = np.abs(Y) * xyz

    cmap = plt.cm.ScalarMappable(
        cmap=XRDCT_palette_cmp  # cmap=plt.get_cmap("PRGn")
    )  # Find out how to use own palette
    cmap.set_clim(-0.5, 0.5)
    divnorm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5)

    if sphere:
        ax.plot_surface(
            xyz[0],
            xyz[1],
            xyz[2],
            facecolors=cmap.to_rgba(Y.real),
            norm=divnorm,
            rstride=2,
            cstride=2,
        )
        ax_lim = 1.25
    else:
        ax.plot_surface(
            Yx,
            Yy,
            Yz,
            facecolors=cmap.to_rgba(Y.real),
            norm=divnorm,
            rstride=2,
            cstride=2,
        )
        ax.set_title(r"$Y_{{{},{}}}$".format(l, m))
        ax_lim = 0.5

    # Draw a set of x, y, z axes for reference.
    # ax_lim = 0.5
    ax.plot([-ax_lim, ax_lim], [0, 0], [0, 0], c="0.5", lw=1, zorder=10)
    ax.plot([0, 0], [-ax_lim, ax_lim], [0, 0], c="0.5", lw=1, zorder=10)
    ax.plot([0, 0], [0, 0], [-ax_lim, ax_lim], c="0.5", lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.

    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis("off")


def plot_exp_sin(ax, A, B, sphere=False):

    theta = np.linspace(0, np.pi, 100)
    phi = np.linspace(0, 2 * np.pi, 100)
    # Create a 2-D meshgrid of (theta, phi) angles.
    theta, phi = np.meshgrid(theta, phi)
    # Calculate the Cartesian coordinates of each point in the mesh.
    xyz = np.array(
        [np.sin(theta) * np.sin(phi), np.sin(theta) * np.cos(phi), np.cos(theta)]
    )
    # COS_2_THETA = np.cos(np.sin(theta) * np.sin(phi)) ** 2
    # SIN_2_THETA = (
    #     1 - np.cos(np.sin(theta) * np.cos(phi) ) ** 2
    # )
    F = A**2 * np.exp(-B * np.sin(theta) ** 2)

    Fx, Fy, Fz = F * xyz

    cmap = plt.cm.ScalarMappable(cmap=XRDCT_palette_cmp)
    cmap.set_clim(-0.5, 0.5)
    divnorm = TwoSlopeNorm(vmin=-0.5, vcenter=0, vmax=0.5)

    if sphere:
        ax.plot_surface(
            xyz[0],
            xyz[1],
            xyz[2],
            facecolors=cmap.to_rgba(F),
            norm=divnorm,
            rstride=2,
            cstride=2,
        )
        ax_lim = 1.25
    else:
        ax.plot_surface(
            Fx, Fy, Fz, facecolors=cmap.to_rgba(F), norm=divnorm, rstride=2, cstride=2
        )
        ax.set_title(r"$F_{{{},{}}}$".format(A, B))
        ax_lim = 0.5

    # Draw a set of x, y, z axes for reference.

    ax.plot([-ax_lim, ax_lim], [0, 0], [0, 0], c="0.5", lw=1, zorder=10)
    ax.plot([0, 0], [-ax_lim, ax_lim], [0, 0], c="0.5", lw=1, zorder=10)
    ax.plot([0, 0], [0, 0], [-ax_lim, ax_lim], c="0.5", lw=1, zorder=10)
    # Set the Axes limits and title, turn off the Axes frame.

    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_zlim(-ax_lim, ax_lim)
    ax.axis("off")


def plot_err_hist(
    fasit_rec,
    AD_rec,
    sym_rec,
    title="",
    bins=100,
    limits=(0, 6.5),
    stacked=True,
    xscale="linear",
):

    metric = [AD_rec, sym_rec]
    subtitles = ["Automatic", "Symbolic"]
    legends = ["a0", "a2", "a4", "a6", r"$\theta$ ($\pi$)", r"$\varphi$ ($\pi$)"]

    fasit_params = [
        fasit_rec.get_1D_array(fasit_rec.params[p])
        for p in range(len(fasit_rec.params))
    ]

    fig, ax = plt.subplots(2, 1, sharey=True, sharex=True)
    for i in range(len(ax)):

        compare_params = [
            metric[i].get_1D_array(metric[i].params[p])
            for p in range(len(metric[i].params))
        ]

        ax[i].set_xlabel(r"Abs error (wrt. $\pi$)")
        ax[i].set_ylabel("Counts")
        ax[i].set_xscale(xscale)
        ax[i].set_title(subtitles[i])

        err = np.zeros((len(fasit_params), len(fasit_params[0]))).T

        for j in range(len(fasit_params)):

            if j > 3:
                err1 = (
                    np.abs(fasit_params[j] % np.pi - compare_params[j] % np.pi) / np.pi,
                )
                err2 = (
                    np.abs(
                        np.pi
                        - np.abs(fasit_params[j] % np.pi - (compare_params[j] % np.pi))
                    )
                    / np.pi
                )

                err[:, j] = np.minimum(err1, err2)
                assert np.all(err[:, j] <= 0.5)

            else:
                err[:, j] = np.abs(fasit_params[j] - compare_params[j])

        #     ax[i].hist(
        #         err,
        #         label=legends[j],
        #         alpha=0.69,
        #         bins=bins,
        #         range=limits,
        #         stacked=stacked,
        #         histtype="bar",
        #     )

        ax[i].hist(
            err, label=legends, alpha=0.69, bins=bins, range=limits, stacked=stacked
        )

        ax[i].legend()

    fig.suptitle(title)
    fig.savefig(title + ".svg")
    plt.show()
    return


def plot_result_heatmap(tt_result_AD, tt_result_symbolic, title, slice=None):

    ylabels = ["a0", "a2", "a4", "a6", "\u03B8", "\u03C6"]
    if slice is None:
        x_length = len(tt_result_AD)
        slice = [0, x_length]
    else:
        x_length = len(
            tt_result_AD.get_1D_array(tt_result_AD.a0)[slice[0] : slice[-1]]
        )  # len(tt_result_AD.a0[slice])

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

    err_a = (np.abs(Ad_a - Sym_a) / np.abs(Sym_a))[:, slice[0] : slice[1]]
    err_theta = (np.abs(Ad_theta % np.pi - Sym_theta % np.pi) / np.abs(np.pi))[
        slice[0] : slice[1]
    ]
    err_theta_2 = np.abs(1 - err_theta)
    err_theta = np.minimum(err_theta, err_theta_2)

    err_phi = (np.abs(Ad_phi % np.pi - Sym_phi % np.pi) / np.abs(np.pi))[
        slice[0] : slice[1]
    ]
    err_phi_2 = np.abs(1 - err_phi)
    err_phi = np.minimum(err_phi, err_phi_2)

    img = np.vstack(
        (err_a[0, :], err_a[1, :], err_a[2, :], err_a[3, :], err_theta, err_phi)
    )  # Bit manual, but but

    fig, ax = plt.subplots(figsize=(20, 8))
    ax.imshow(img, cmap=XRDCT_palette_cmp, vmin=0, vmax=1)
    ax.set_yticks(np.arange(len(ylabels)), labels=ylabels)
    ax.set_xlabel("Voxel number")

    for i in range(x_length):
        for j in range(num_coeffs):
            text = ax.text(
                i,
                j,
                f"{err_a[j, i]:.2e}",
                ha="center",
                va="center",
                color="w",
            )
        text = ax.text(
            i, j + 1, f"{err_theta[i]:.2e}", ha="center", va="center", color="w"
        )
        text = ax.text(
            i, j + 2, f"{err_phi[i]:.2e}", ha="center", va="center", color="w"
        )
    # cbar = ax.figure.colorbar(img, ax=ax, location="bottom")
    # cbar.ax.set_ylabel("Relative error", rotation=-90, va="bottom")
    # ax.figure.colorbar(img, ax=ax, orientation="horizontal")
    ax.set_title(title)
    fig.savefig(title + ".svg")
    plt.show()


def plot_angles_distribution(
    fasit_rec, AD_rec, sym_rec, title="", bins=30, limits=(0, np.pi)
):

    fasit_params = [
        fasit_rec.get_1D_array(fasit_rec.params[p])
        for p in range(len(fasit_rec.params))
    ]
    AD_params = [
        AD_rec.get_1D_array(AD_rec.params[p]) for p in range(len(AD_rec.params))
    ]
    sym_params = [
        sym_rec.get_1D_array(sym_rec.params[p]) for p in range(len(sym_rec.params))
    ]

    titles = ["a0", "a2", "a4", "a6", r"$\theta$", r"$\varphi$"]

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for i, ax in enumerate(np.reshape(axs, -1)):
        ax.set_xticks(np.arange(0, 5 * np.pi / 4, np.pi / 4))
        ax.set_xticklabels(["0", "1/4", "1/2", "3/4", "1"])

        ax.set_xlabel(r"Value [$\pi$]")
        ax.set_ylabel("Counts")
        ax.set_title(titles[i + 4])
        ax.hist(
            fasit_params[i + 4] % np.pi,
            label="Solution",
            alpha=0.69,
            bins=bins,
            range=limits,
            # stacked=True,
        )
        ax.hist(
            AD_params[i + 4] % np.pi,
            label="Automatic",
            alpha=0.69,
            bins=bins,
            range=limits,
            # stacked=True,
        )
        ax.hist(
            sym_params[i + 4] % np.pi,
            label="Symbolic",
            alpha=0.69,
            bins=bins,
            range=limits,
            # stacked=True,
        )
        ax.legend(loc="upper left")

    fig.suptitle(title)
    fig.savefig(title + ".svg")
    plt.show()
    return


def plot_convergence_curve(AD_rec, sym_rec, title="", scale="linear"):

    x1 = np.arange(1, len(AD_rec.error_data) + 1)
    x2 = np.arange(1, len(sym_rec.error_data) + 1)
    fig, ax = plt.subplots()
    ax.plot(
        x1, AD_rec.error_data, label=f"Automatic {np.float64(AD_rec.timing_data):.1f} s"
    )
    ax.plot(
        x2,
        sym_rec.error_data,
        label=f"Symbolic {np.float64(sym_rec.timing_data):.1f} s",
    )
    ax.legend()
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Projection Residual")
    ax.set_title(title)
    ax.set_xscale(scale)
    ax.set_yscale(scale)
    fig.savefig(title + ".svg")
    plt.show()
    return


def plot_tt_3D_quiver(tt_result):

    return


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
