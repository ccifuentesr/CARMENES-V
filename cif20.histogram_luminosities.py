# %%

# PLOT: HISTOGRAM > LUMINOSITIES
# Cifuentes et al. 2020

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_hist_luminosities'

# Mother
# Appending data for each spectral type (KML).

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    Lbol_K = []
    Lbol_M = []
    Lbol_L = []
    for row in Mother:
        if row['Sample_source'] == 'REC100' or row['Sample_source'] == 'Lep13' or row['Sample_source'] == 'AF15':
            try:
                Lbol_K.append(float(row['Lbol']))
            except ValueError:
                next(Mother)
        if row['Sample_source'] == 'Karmn':
            try:
                Lbol_M.append(float(row['Lbol']))
            except ValueError:
                next(Mother)
        if row['Sample_source'] == 'Sma17':
            try:
                Lbol_L.append(float(row['Lbol']))
            except ValueError:
                next(Mother)

# PLOTTING

# Labels

xlabel = r'$L [L_\odot]$'
ylabel = r'Number of stars'

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

# Plot: histogram

plt.hist(Lbol_K, bins=np.logspace(np.log10(min(Lbol_K)), np.log10(max(Lbol_K)), int(
    len(Lbol_K)*0.15)), density=False, histtype='step', color=cmap(0.9), label='K', linewidth=linewidth)
plt.hist(Lbol_M, bins=np.logspace(np.log10(min(Lbol_M)), np.log10(max(Lbol_M)), int(len(
    Lbol_M)*0.02)), density=False, histtype='step', color=cmap(0.6), label='M', linewidth=linewidth, zorder=0)
plt.hist(Lbol_L, bins=np.logspace(np.log10(min(Lbol_L)), np.log10(max(Lbol_L)), int(len(
    Lbol_L)*0.15)), density=False, histtype='step', color=cmap(0.3), label='L', linewidth=linewidth)

# Aesthetics

# Axes: range & scale

ax.set_xscale('log')

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)

# Axes: ticks

ax.set_xticks((0.0001, 0.001, 0.01, 0.1, 1))

ax.tick_params(axis='x', labelsize=tickssize-2, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
