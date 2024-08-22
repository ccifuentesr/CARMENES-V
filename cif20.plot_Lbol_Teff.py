# %%

# PLOT: LUMINOSITY - EFFECTIVE TEMPERATURE
# Cifuentes et al. 2020

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_plot_Lbol_Teff'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

x_axis = 'Teff'
y_axis = 'Lbol'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    SpT_young = []
    d_pc = []
    d_pc_ = []
    d_pc_young = []
    Teff = []
    Teff_ = []
    Teff_young = []
    Lbol = []
    Lbol_ = []
    Lbol_young = []
    eLbol = []
    eTeff = []
    for row in reader:
        if row['d_pc'] != '' and row['Lbol'] != '' and row['Teff'] != '':
            SpT_.append(float(row['SpTnum']))
            d_pc_.append(float(row['d_pc']))
            Lbol_.append(float(row['Lbol']))
            Teff_.append(float(row['Teff']))
            if row[booleans[1]] == 'true':  # Young stars
                SpT_young.append(float(row['SpTnum']))
                d_pc_young.append(float(row['d_pc']))
                Lbol_young.append(float(row['Lbol']))
                Teff_young.append(float(row['Teff']))
            if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                SpT.append(float(row['SpTnum']))
                d_pc.append(float(row['d_pc']))
                Lbol.append(float(row['Lbol']))
                eLbol.append(float(row['Lberr']))
                Teff.append(float(row['Teff']))

# Variables

x_cif_ = Teff_
y_cif_ = Lbol_

x_cif_young = Teff_young
y_cif_young = Lbol_young

x_cif = Teff
y_cif = Lbol
z_cif = SpT

error_bars = [4000, 1E-4]
x_error = 50  # K
y_error = np.mean(eLbol)

# PLOTTING

# Labels

xlabel = r'T$_{\rm eff}$ [K]'
ylabel = r'$L$ [$L_{\rm sol}$]'
cbarlabel = r'Spectral type'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (12, 10)
pointsize = 30
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma_r')

# Plots: distribution

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap, zorder=1)

# Plots: Young stars sample

ax.scatter(x_cif_young, y_cif_young, marker='o', s=pointsize,
           facecolors='', edgecolors='k', zorder=2)

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

plt.ylim(1e-5, 1e0)
ax.set_yscale('log')
ax.invert_xaxis()

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)  # Colorbar
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
