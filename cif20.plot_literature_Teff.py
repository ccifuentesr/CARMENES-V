# %%

# PLOT: LITERATURE > EFFECTIVE TEMPERATURE
# Cifuentes et al. 2020
# Passegger et al. 2019

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams.update({'errorbar.capsize': 4})


# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_literature_Teff'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    Teff = []
    TeffVN = []
    e_TeffVN = []
    TeffV = []
    e_TeffV = []
    TeffN = []
    e_TeffN = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
            if row['Teff'] != '':
                if row['Teff_Pas19'] != '':
                    SpT.append(float(row['SpTnum']))
                    Teff.append(float(row['Teff']))
                    TeffVN.append(float(row['Teff_Pas19']))
                    e_TeffVN.append(float(row['eTeff_Pas19']))

# Variables

x = Teff
y = TeffVN
z = SpT

xerr = [50 for i in range(len(x))]  # VOSA
yerr = e_TeffVN

xp = np.linspace(1000, 5000, 10000)
yp = xp

# PLOTTING

# Labels

xlabel = r'T$_{\rm eff}$ [K] $_{\rm This~work}$'
ylabel = r'T$_{\rm eff}$ [K] $_{\rm Pas19}$'
cbarlabel = r'Spectral type'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma_r')

# Plots: distribution

# Mapping colors for errorbars
norm = matplotlib.colors.Normalize(vmin=min(y), vmax=max(y), clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap='magma')
x_color = np.array([(mapper.to_rgba(v)) for v in y])

sc = ax.scatter(x, y, c=z, cmap='magma_r', s=pointsize,
                marker='o', zorder=1)  # INVISIBLE

# Loop over each data point to plot
for x, y, ex, ey, color in zip(x, y, xerr, yerr, x_color):
    plt.errorbar(x, y,  xerr=ex, yerr=ey, lw=1,
                 capsize=3, color='gainsboro', zorder=0)

# Plots: equality

plt.plot(xp, yp, 'k--', lw=2, zorder=0)

# Aesthetics

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

# Axes: range & scale
# Limited in x by Sch19

ax.set_xlim(2700, 4300)
ax.set_ylim(2700, 4300)


# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax)  # Colorbar
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)
# cbar.set_ticks(np.arange(-2, 19, 2))
# cbar.ax.set_yticklabels(SpT_half)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
