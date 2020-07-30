# %%

# PLOT: METALLICITY
# Cifuentes et al. 2020
# Literature values against BT-Settl output.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import csv
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams.update({'errorbar.capsize': 4})

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif01_plot_box_meta'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Meta_lit = []
    Meta_VOSA = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
            if row['FeH_lit'] != '' and row['Meta_meta'] != '':
                Meta_lit.append(float(row['FeH_lit']))
                Meta_VOSA.append(float(row['Meta_meta']))

# Mother (Bis)

with open('cif01.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Meta_lit15 = []
    Meta_lit10 = []
    Meta_lit05 = []
    Meta_lit00 = []
    Meta_lit003 = []
    Meta_lit005 = []
    Meta_VOSA15 = []
    Meta_VOSA10 = []
    Meta_VOSA05 = []
    Meta_VOSA00 = []
    Meta_VOSA003 = []
    Meta_VOSA005 = []
    for row in reader:
        if row['FeH_lit'] != '':
            if row['Meta_meta'] == '-1.5':
                Meta_lit15.append(float(row['FeH_lit']))
                Meta_VOSA15.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '-1.0':
                Meta_lit10.append(float(row['FeH_lit']))
                Meta_VOSA10.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '-0.5':
                Meta_lit05.append(float(row['FeH_lit']))
                Meta_VOSA05.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '0.0':
                Meta_lit00.append(float(row['FeH_lit']))
                Meta_VOSA00.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '0.3':
                Meta_lit003.append(float(row['FeH_lit']))
                Meta_VOSA003.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '0.5':
                Meta_lit005.append(float(row['FeH_lit']))
                Meta_VOSA005.append(float(row['Meta_meta']))

# Variables

x_lit = [Meta_lit15, Meta_lit10, Meta_lit05,
         Meta_lit00, Meta_lit003, Meta_lit005]
x_VOSA = [Meta_VOSA15, Meta_VOSA10, Meta_VOSA05,
          Meta_VOSA00, Meta_VOSA003, Meta_VOSA005]

# PLOTTING

# Labels

xlabel = r'[Fe/H]$_{\rm BT-Settl}$'
ylabel = r'[Fe/H]$_{\rm Literature}$'

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('copper')

# Boxplot variables

medianprops = dict(linestyle='-', linewidth=3.0)
c_box = 'white'
c_whisker = 'black'
c_flier = 'grey'
c_median = 'blue'
# c_median = cmap(3/len(labels))

# Plot: box

# (len(x_lit[0])/len(x_lit[5]))
widths = (len(x_lit[0])/len(x_lit[5]), len(x_lit[1])/len(x_lit[5]), len(x_lit[2])/len(x_lit[5]),
          len(x_lit[3])/len(x_lit[5]), len(x_lit[4])/len(x_lit[5]), len(x_lit[5])/len(x_lit[5]))

plt.boxplot(x_lit, showfliers=True, notch=True, widths=widths, patch_artist=True,
            boxprops=dict(facecolor=c_box, color=c_whisker),
            capprops=dict(color=c_whisker),
            whiskerprops=dict(color=c_whisker),
            flierprops=dict(color=c_flier, markeredgecolor=c_flier),
            medianprops=dict(color=c_median, linewidth=3)
            )

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
labels = [item.get_text() for item in ax.get_xticklabels()]
# labels[0] = r'BT-Settl'
# labels[1] = r'Literature'
labels = ['$-1.5$', '$-1.0$', '$-0.5$', '$0.0$', '$+0.3$', '$+0.5$']
# labels = [-1.5, -1.0, -0.5, 0.0, 0.3, 0.5, 1.0]
ax.set_xticklabels(labels)

# Axes: ticks

# ticks = [1, 2]
# ax.set_xticks(ticks)
ax.tick_params(axis='x', labelsize=tickssize,
               top=False, bottom=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=False, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
# ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)

# Axes: range & scale

# ax.set_xlim(-1.5, 1.0)
# ax.set_ylim([1E-3, 2E-1])
# plt.ylim(9E-4, 1)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.invert_xaxis()
# ax.invert_yaxis()

# Show & Save

plt.savefig(filename+'.png', bbox_inches='tight')
plt.show()
