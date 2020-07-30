# %%

# PLOT: HISTOGRAM > DISTANCES
# Cifuentes et al. 2020

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_hist_distance'

Bool = ['Bool_dphot']

# Mother
# Appending data for each spectral type (KML).

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    d_pc_K = []
    d_pc_M = []
    d_pc_L = []
    for row in Mother:
        if row[Bool[0]] == 'false':
            if float(row['SpTnum']) < 0:
                try:
                    d_pc_K.append(float(row['d_pc']))
                except ValueError:
                    next(Mother)
            if float(row['SpTnum']) >= 0 and float(row['SpTnum']) < 10:
                try:
                    d_pc_M.append(float(row['d_pc']))
                except ValueError:
                    next(Mother)
            if float(row['SpTnum']) >= 10:
                try:
                    d_pc_L.append(float(row['d_pc']))
                except ValueError:
                    next(Mother)


# PLOTTING

# Labels

xlabel = r'$d$ (pc)'
ylabel = r'Number of stars'
cbarlabel = r'Spectral type'

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18
linewidth = 3

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma')

# Plots: histogram

plt.hist(d_pc_K, bins=np.logspace(np.log10(min(d_pc_K)), np.log10(max(d_pc_K)), int(
    len(d_pc_K)*0.15)), density=False, histtype='step', color=cmap(0.9), label='K', linewidth=linewidth)
plt.hist(d_pc_M, bins=np.logspace(np.log10(min(d_pc_M)), np.log10(max(d_pc_M)), int(len(
    d_pc_M)*0.02)), density=False, histtype='step', color=cmap(0.6), label='M', linewidth=linewidth, zorder=0)
plt.hist(d_pc_L, bins=np.logspace(np.log10(min(d_pc_L)), np.log10(max(d_pc_L)), int(len(
    d_pc_L)*0.15)), density=False, histtype='step', color=cmap(0.3), label='L', linewidth=linewidth)

# Aesthetics

# Axes: range & scale

ax.set_xscale('log')
ax.set_yscale('log')

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)

ax.set_xticks((1, 2, 3, 5, 10, 20, 30, 50, 100))

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
