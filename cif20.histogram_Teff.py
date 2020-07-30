# %%

# PLOT: HISTOGRAM > EFFECTIVETEMPERATURES
# Cifuentes et al. 2020

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib import ticker
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data

Mother_version = '01'
filename = 'cif20_hist_Teff'

# Mother
# Appending data for each spectral type (KML).

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    Teff_K = []
    Teff_M = []
    Teff_L = []
    for row in Mother:
        if row['Sample_source'] == 'REC100' or row['Sample_source'] == 'Lep13' or row['Sample_source'] == 'AF15':
            try:
                Teff_K.append(float(row['Teff']))
            except ValueError:
                next(Mother)
        if row['Sample_source'] == 'Karmn':
            try:
                Teff_M.append(float(row['Teff']))
            except ValueError:
                next(Mother)
        if row['Sample_source'] == 'Sma17':
            try:
                Teff_L.append(float(row['Teff']))
            except ValueError:
                next(Mother)

# PLOTTING

# Labels

xlabel = r'$T_{\rm eff}$ [K]'
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

# Plots: histogram

plt.hist(Teff_K, bins=list(np.arange(min(Teff_K), max(Teff_K), 100)),
         density=False, histtype='step', color=cmap(0.9), label='K', linewidth=linewidth)
plt.hist(Teff_M, bins=list(np.arange(min(Teff_M), max(Teff_M), 100)), density=False,
         histtype='step', color=cmap(0.6), label='M', linewidth=linewidth, zorder=0)
plt.hist(Teff_L, bins=list(np.arange(min(Teff_L), max(Teff_L), 100)),
         density=False, histtype='step', color=cmap(0.3), label='L', linewidth=linewidth)

"""
Aesthetics
----------
"""

# Axes: range & scale

# ...

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)

# Axes: ticks

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
