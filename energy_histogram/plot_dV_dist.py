#!/usr/bin/env python3
import matplotlib
import os
import argparse
matplotlib.use("pgf")
import matplotlib.pyplot as plt
plt.rcParams.update({
    "pgf.texsystem": "lualatex",
    "font.family": "serif",  # use serif/main font for text elements
    "text.usetex": True,     # use inline math for ticks
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

data = np.genfromtxt('dV_dist.dat', unpack=True)
w, h = figaspect(1/1.1)
plt.figure(figsize = (w,h))
x = data[0]
y = data[1]
plt.bar(x, y, width=0.2, edgecolor='black', linewidth=1.5)
plt.xlabel('$V$ (kcal/mol)')
plt.ylabel('$p(V)$')
ax = plt.gca()
ax.set_xlim(0, 10)
ax.tick_params(direction = 'in', which = 'major', length=6.0, width = 1.0, top = True, right = True)
ax.tick_params(direction = 'in', which = 'minor', length=3.0, width = 1.0, top = True, right = True)
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False
ax.xaxis.set_major_locator(plt.FixedLocator(np.linspace(0, 10, 6)))
ax.xaxis.set_minor_locator(plt.FixedLocator(np.linspace(0, 10, 21)))
ax.yaxis.set_minor_locator(AutoMinorLocator())
#plt.text(5, 0.25, r'$\gamma=6.47\times 10^{-3}$')
plt.savefig('dV_dist.png', dpi=300, bbox_inches='tight', transparent=False)
