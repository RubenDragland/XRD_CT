

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.cm import get_cmap
from cycler import cycler
import h5py
import scipy.io


#print(plt.style.available)


#import plotly

# mpl.style.use("fast")

# mpl.rcParams["axes.prop_cycle"] = cycler(
    # color=plt.style.library["ggplot"]["axes.prop_cycle"].by_key()["color"])

plt.style.use([ "tableau-colorblind10", "seaborn-paper"])
mpl.rcParams["axes.prop_cycle"] = ( cycler(color = [u"#F05039", u"#E57A77", u"#EEBAB4", u"#1F449C", u"#3D65A5", u"#7CA1CC", u"#A8B6CC"]) + cycler(linestyle=["-", "--", "-.", ":", "-", "--", "-."]) )  #From: https://www.datylon.com/blog/data-visualization-for-colorblind-readers  and https://ranocha.de/blog/colors/ , respectively
#cycler(color=plt.style.library["tab10"]["axes.prop_cycle"].by_key()["color"])
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
mpl.rcParams['axes.formatter.use_mathtext'] = True
#mpl.rcParams['text.usetex'] = True # Use Latex
"""
mpl.rcParams['legend.loc'] = "upper_right" # Suggestion
mpl.rcParams['']
mpl.rcParams['']
"""


def plot_visualisation_grads(mat_paths: list, keys: list, titles: list, bins=200, range=(-4e-5,-1e-9)  )-> None:
    """
    Loads a .mat-file as a python dict. 
    Extracts data, reshapes, and executes visualisation in the form of a 1D histogram

    """

    fig, axs = plt.subplots(3,1, sharex=True, sharey=True, figsize= (4,12))
    cycler = plt.rcParams["axes.prop_cycle"].by_key()["color"] # Use this to retrieve cycle
    for i, (ax, path) in enumerate( zip( axs, mat_paths ) ):
        mat_dict = scipy.io.loadmat(path)
        data = np.reshape(mat_dict[ keys[i] ], -1 ) #* 1e5 # Scaled, though a bit hard-coded. 

        ax.hist(data, bins=bins, range = range, density = False, alpha = 0.69 )   
        ax.ticklabel_format( style = "sci", scilimits = (0, 0))     
        ax.set_ylabel("Count")
        ax.set_title(titles[i])

        t = ax.yaxis.get_offset_text()
        t.set_x(-0.2)

    ax.set_xlabel("Gradient value")
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.show()
    return

def plot_comparison_grads(mat_paths: list, keys: list, labels: list, bins = 200, range = (-4e-5, -1e-9))-> None:
    """
    Plots gradient fingerprints in same figure for comparison purpouses
    """
    colors = [u"#F05039", u"#F05039",  u"#1F449C", u"#A8B6CC"]
    fig, ax = plt.subplots()
    for i , path in enumerate(mat_paths[1:]):
        i+=1

        mat_dict = scipy.io.loadmat(path)
        data = np.reshape(mat_dict[ keys[i] ], -1 ) 

        ax.hist(data, bins=bins, range = range, density = False, alpha = 0.69, color = colors[i], label = labels[i] )   


    ax.ticklabel_format( style = "sci", scilimits = (0, 0))   
    ax.set_ylabel("Count")
    t = ax.yaxis.get_offset_text()
    t.set_x(-0.2)
    ax.set_xlabel("Gradient value")
    ax.legend()
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    plt.show()
    return





#print(plt.rcParams)
#print(plt.style.library.keys())