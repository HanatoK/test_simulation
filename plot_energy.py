#!/usr/bin/env python3
import matplotlib
import os
import argparse
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "Arimo",  # use serif/main font for text elements
    "text.usetex": False,     # use inline math for ticks
    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
    "axes.labelsize": 28,
    "axes.linewidth": 2.0,
    "font.size": 24,
    "pgf.preamble": '\n'.join([
         "\\usepackage{units}",          # load additional packages
         "\\usepackage{metalogo}",
         "\\usepackage{unicode-math}",   # unicode math setup
         r"\setmathfont{MathJax_Math}",
         r"\setmainfont{Arimo}",  # serif font via preamble
         ])
})
import numpy as np
from matplotlib.figure import figaspect
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator

kineticEnergy = []
potentialEnergy = []
with open('XY.traj', 'r') as fInput:
    for line in fInput:
        fields = line.split()
        if fields:
            if fields[0] != '#':
                kineticEnergy.append(float(fields[9]))
                potentialEnergy.append(float(fields[10]))
totalEnergy_all = np.array(kineticEnergy) + np.array(potentialEnergy)
totalEnergy = totalEnergy_all[::1000]
w, h = figaspect(1/2)
plt.figure(figsize = (w,h))
X = np.linspace(0, len(totalEnergy), len(totalEnergy))
plt.plot(X, totalEnergy)
plt.xlabel(r'X')
plt.ylabel(r'Energy')
plt.savefig('energy.png', bbox_inches='tight', dpi=300, transparent=False)
plt.close()

plt.figure(figsize = (w,h))
plt.hist(totalEnergy_all, bins=40, density=True)
plt.xlabel(r'Total energy')
plt.title('Histogram')
plt.savefig('hist.png', bbox_inches='tight', dpi=300, transparent=False)
plt.close()
