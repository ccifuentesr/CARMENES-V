# %%

# PLOT: LITERATURE > ABSOLUTE MAGNITUDE - EFFECTIVE TEMPERATURE
# MJ vs. Teff
# Cifuentes et al. 2020
# Dahn et al. 2002
# Lepine et al. 2014
# Gaidos et al. 2014

from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter
import csv
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
sys.path

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_literature_Mabs_Teff'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

x_axis = 'Teff'
y_axis = 'J'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    d_pc = []
    d_pc_ = []
    Teff = []
    Teff_ = []
    Lbol = []
    Lbol_ = []
    eLbol = []
    eTeff = []
    mag = []
    mag_ = []
    for row in reader:
        if row[y_axis+'_mag'] != '':
            if row['d_pc'] != '' and row['Lbol'] != '' and row['Teff'] != '':
                SpT_.append(float(row['SpTnum']))
                d_pc_.append(float(row['d_pc']))
                Lbol_.append(float(row['Lbol']))
                Teff_.append(float(row['Teff']))
                mag_.append(float(row[y_axis+'_mag']))
                if row['Qf_'+y_axis] == 'false':
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        SpT.append(float(row['SpTnum']))
                        d_pc.append(float(row['d_pc']))
                        Lbol.append(float(row['Lbol']))
                        eLbol.append(float(row['Lberr']))
                        Teff.append(float(row['Teff']))
                        mag.append(float(row[y_axis+'_mag']))

# Lepine et al. 2013
# plx_1: Trigonometric; _2: Spectroscopic; _3: Photometric.
# Teff from PHOENIX models.

with open('Literature/lepine13.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    plx_lep13 = []
    eplx_lep13 = []
    Teff_lep13 = []
    mag_lep13 = []
    for row in reader:
        if row[y_axis+'mag'] != '' and row['plx'] != '' and row['Teff'] != '':
            plx_lep13.append(float(row['plx']))
            eplx_lep13.append(float(row['e_plx']))
            Teff_lep13.append(float(row['Teff']))
            mag_lep13.append(float(row[y_axis+'mag']))

# Gaidos et al. 2014
# Chosing spectroscopic determinations of Teff only ('PHX').
# Teff from PHOENIX models.

with open('Literature/gaidos14.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    plx_gai14 = []
    eplx_gai14 = []
    Teff_gai14 = []
    eTeff_gai14 = []
    emag_gai14 = []
    mag_gai14 = []
    for row in reader:
        if row[y_axis+'mag'] != '' and row['e_'+y_axis+'mag'] != '' and row['plx'] != '' and row['Teff'] != '':
            if row['n_Teff'] == "PHX":
                plx_gai14.append(float(row['plx']))
                eplx_gai14.append(float(row['e_plx']))
                Teff_gai14.append(float(row['Teff']))
                eTeff_gai14.append(float(row['e_Teff']))
                mag_gai14.append(float(row[y_axis+'mag']))
                emag_gai14.append(float(row['e_'+y_axis+'mag']))

# Dahn et al. 2002
# Table 5: Effective temperatures

Mabs_dah02 = [10.05, 10.43, 10.46, 10.47, 10.72, 10.77, 10.83, 10.91, 10.98, 10.88, 11.13, 11.34, 11.44, 11.44, 11.63,
              11.49, 11.55, 11.42, 11.66, 11.77, 12.03, 11.72, 11.97, 11.81, 11.78, 11.87, 12.24, 12.88, 12.58, 12.90,
              12.71, 12.95, 13.20, 13.23, 13.38, 13.77, 13.56, 13.49, 13.85, 14.68, 14.79, 14.97, 14.89]

Teff_dah02 = [2957, 2813, 2811, 2800, 2707, 2687, 2676, 2645, 2638, 2613, 2583, 2496, 2495, 2471, 2467, 2451, 2441, 2439,
              2367, 2364, 2354, 2319, 2319, 2302, 2302, 2273, 2138, 2108, 2092, 2031, 1993, 1945, 1923, 1885, 1853, 1792,
              1734, 1703, 1601, 1429, 1376, 1372, 1335]

# Variables

Mabs_ = [mag_[i] - 5*np.log10(d_pc_[i]) + 5 for i in range(len(mag_))]
Mabs = [mag[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag))]

d_pc_gai14 = [1/plx_gai14[i] for i in range(len(plx_gai14))]
Mabs_gai14 = [mag_gai14[i] - 5 *
              np.log10(d_pc_gai14[i]) + 5 for i in range(len(mag_gai14))]

d_pc_lep13 = [1/plx_lep13[i] for i in range(len(plx_lep13))]
Mabs_lep13 = [mag_lep13[i] - 5 *
              np.log10(d_pc_lep13[i]) + 5 for i in range(len(mag_lep13))]

x_cif_ = Teff_
y_cif_ = Mabs_

x_cif = Teff
y_cif = Mabs
z_cif = SpT

x_dah02 = Teff_dah02
y_dah02 = Mabs_dah02

x_gai14 = Teff_gai14
y_gai14 = Mabs_gai14

x_lep13 = Teff_lep13
y_lep13 = Mabs_lep13

# PLOTTING

# Labels

xlabel = r'T$_{\rm eff}$ [K]'
ylabel = r'M$_J$ [mag]'
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
ax.scatter(x_dah02, y_dah02, s=pointsize,
           facecolors='', edgecolors='blue', zorder=2)
ax.scatter(x_lep13, y_lep13, s=pointsize,
           facecolors='', edgecolors='green', alpha=0.5, zorder=4)
ax.scatter(x_gai14, y_gai14, s=pointsize,
           facecolors='', edgecolors='magenta', alpha=0.5, zorder=3)

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

plt.ylim(3, 16)
ax.invert_xaxis()
ax.invert_yaxis()

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT
# Only keep (#) to hide colorbar with consistent colouring.

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
