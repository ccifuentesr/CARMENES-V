# %%

# PLOT: HISTOGRAM > RUWE
# Cifuentes et al. 2020

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import matplotlib.ticker as mticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data

Mother_version = '01'
filename = 'cif20_hist_RUWE'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    RUWE = []
    for row in Mother:
        if row['RUWE'] != '':
            RUWE.append(float(row['RUWE']))

# PLOTTING

# Labels

xlabel = r'RUWE'
ylabel = r'Number of stars'
cbarlabel = r'Spectral type'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

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

# Plots: histogram

plt.hist(RUWE, bins=np.logspace(np.log10(min(RUWE)), np.log10(max(RUWE)), 80),
         density=False, alpha=1, histtype='step', color='blue', linewidth=linewidth)
plt.axvline(1.41, color='red', linestyle='dashed', linewidth=linewidth)

# Aesthetics

# Axes: range & scale

ax.set_xscale('log')

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
ax.xaxis.get_major_formatter().set_scientific(False)

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.minorticks_on()

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
