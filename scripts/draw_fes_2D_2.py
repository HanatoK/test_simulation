#!/usr/bin/env python3
import matplotlib
import argparse
import pandas as pd
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import matplotlib.pyplot as plt
matplotlib.use("Agg")
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "FreeSans",  # use serif/main font for text elements
    "mathtext.fontset": "stix",
    "text.usetex": False,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
    'axes.unicode_minus': False
})


# parser = argparse.ArgumentParser()
# parser.add_argument("pmf", help = "specify the PMF file")
# parser.add_argument("-o", "--output", default = "output.png", help = "specify the PNG output image file")
# parser.add_argument("--xtitle", default = "X", help = "title along X axis")
# parser.add_argument("--ytitle", default = "Y", help = "title along Y axis")
# parser.add_argument("--levels", default = 25, type = int, help = "number of levels")
# args = parser.parse_args()

def plotfes(pmffilename, pngfilename, xtitle, ytitle, title=''):
    x, y, z = np.genfromtxt(pmffilename, unpack=True)
    # z = z
    z = z - np.min(z)
    # w, h = figaspect(1/1.75)
    # plt.figure(figsize=(w, h))
    # z = np.clip(z, 0, 16.0)
    path = pd.read_csv('../cmake-build-release/PCV_bias_10_10_pcv.traj', delimiter=r'\s+', comment='#', header=None)
    path = path[::10]
    binx = len(set(x))
    biny = len(set(y))
    xi = x.reshape(binx, biny)
    yi = y.reshape(binx, biny)
    zi = z.reshape(binx, biny)
    cf = plt.contourf(xi, yi, zi, np.linspace(0, 12, 49), cmap='nipy_spectral')
    plt.scatter(path[3], path[4], s=0.5, alpha=0.6, color='black')
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    #plt.title('Free Energy Surface')
    ax = plt.gca()
    ax.set_title(title)
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
    ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    ax.yaxis.set_major_locator(plt.MaxNLocator(6))
    clb = plt.colorbar(cf)
    clb.ax.set_title("kcal/mol", fontsize=20, pad=10.0)
    clb.ax.xaxis.get_major_formatter()._usetex = False
    clb.ax.yaxis.get_major_formatter()._usetex = False
    # clbticks = [float(i.get_position()[1]) for i in clb.ax.get_yticklabels()]
    # clbticksstr = ['{:.1f}'.format(i) for i in clbticks]
    # print(clbticksstr)
    # clb.ax.set_yticklabels(clbticksstr, fontsize = 20)
    clb.ax.yaxis.set_major_locator(plt.MultipleLocator(2))
    plt.savefig(pngfilename, dpi=400, bbox_inches='tight', transparent=False)
    plt.close()
    return


plotfes('bias_100_10_300000000.czar.grad.pmf', 'pmf_100_10.png', 'X', 'Y', r'$\gamma_x/\gamma_y=10.0$')
plotfes('bias_10_100_300000000.czar.grad.pmf', 'pmf_10_100.png', 'X', 'Y', r'$\gamma_x/\gamma_y=0.1$')
plotfes('bias_10_10_300000000.czar.grad.pmf', 'pmf_10_10.png', 'X', 'Y', r'$\gamma_x/\gamma_y=1.0$')
