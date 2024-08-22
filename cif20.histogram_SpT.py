# %%

# PLOT: HISTOGRAM > SPECTRAL TYPES
# Cifuentes et al. 2020

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
import random

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data

Mother_version = '01'
filename = 'cif20_hist_SpT'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    SpTnum = []
    for row in Mother:
        try:
            SpTnum.append(float(row['SpTnum']))
        except ValueError:
            next(Mother)

SpTypes = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
           'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', ]

# PLOTTING

# Labels

xlabel = r'Spectral type'
ylabel = r'Number of stars'

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18
bins_size = np.arange(-2, len(SpTypes)-1, 1)

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cm = plt.cm.get_cmap('magma_r')

# Plots: histogram

n, bins, patches = plt.hist(SpTnum, bins=bins_size, color='green')
bin_centers = 1 * (bins_size[:-1] + bins_size[1:])

col = bin_centers - min(bin_centers)  # Scale values to interval [0,1]
col = col/max(col)

for c, p in zip(col, patches):
    plt.setp(p, 'facecolor', cm(c))


# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.set_xticklabels(SpTypes)

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize-2, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=False, top=False)
plt.xticks(np.arange(-2, 21, 1))
plt.setp(ax.xaxis.get_majorticklabels(), rotation=-45, ha="left")

# Axes: range & scale

plt.xlim(-2, 19)
plt.ylim(1e0, 1e3)
ax.set_yscale('log')

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
