# %%

# PLOT: HISTOGRAM > LOGG
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
filename = 'cif20_hist_logg'

# Mother
# Appending data for each spectral type (KML).

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    logg_K = []
    logg_M = []
    logg_L = []
    for row in Mother:
        if row['Sample_source'] == 'REC100' or row['Sample_source'] == 'Lep13' or row['Sample_source'] == 'AF15':
            try:
                logg_K.append(float(row['logg']))
            except ValueError:
                next(Mother)
        if row['Sample_source'] == 'Karmn':
            try:
                logg_M.append(float(row['logg']))
            except ValueError:
                next(Mother)
        if row['Sample_source'] == 'Sma17':
            try:
                logg_L.append(float(row['logg']))
            except ValueError:
                next(Mother)

# Variables

x = [logg_K, logg_M, logg_L]

n45 = []
n50 = []
n55 = []

for i in range(len(x)):
    n45.append(x[i].count(4.5))
    n50.append(x[i].count(5.0))
    n55.append(x[i].count(5.5))

# PLOTTING

# Labels


xlabel = r'$\log{g}$ [dex]'
ylabel = r'Number of stars'

index = [4.5, 5.0, 5.5]
width = 0.1       # the width of the bars: can also be len(x) sequence
sep = 0.12        # the separation of the bars

index0 = [index[i]-sep for i in range(len(index))]
index1 = [index[i]+sep for i in range(len(index))]

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
colours = [0.9, 0.6, 0.3]

# Plots: stacked bars
# 4.5, 5.0 and 5.5

plt.bar(index1, n45, width, edgecolor=cmap(
    0.9), fill=False, linewidth=linewidth)
plt.bar(index, n50, width, edgecolor=cmap(
    0.6), fill=False, linewidth=linewidth)
plt.bar(index0, n55, width, edgecolor=cmap(
    0.3), fill=False, linewidth=linewidth)

# Aesthetics

# Axes: range & scale

plt.xlim(4, 6)

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)

# Axes: ticks

ax.set_xticks((4.5, 5.0, 5.5))

ax.tick_params(axis='x', labelsize=tickssize-2, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.xaxis.set_tick_params(which='minor', bottom=False, top=False)
ax.minorticks_on()

ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
