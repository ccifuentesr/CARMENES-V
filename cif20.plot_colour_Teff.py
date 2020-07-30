# %%

# PLOT: COLOUR - EFFECTIVE TEMPERATURE
# G-RP vs. Teff
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
plt.rcParams.update({'errorbar.capsize': 4})

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif01_literature_colour_Teff'
booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

x_axis = 'Teff'
y_axis = ['GG', 'RP']


# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_ = []
    Teff_ = []
    SpT = []
    Teff = []
    mag0 = []
    mag1 = []
    mag0_ = []
    mag1_ = []
    for row in reader:
        if row[x_axis] != '' and row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '':
            SpT_.append(float(row['SpTnum']))
            Teff_.append(float(row[x_axis]))
            mag0_.append(float(row[y_axis[0]+'_mag']))
            mag1_.append(float(row[y_axis[1]+'_mag']))
            if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == 'false':
                if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                    SpT.append(float(row['SpTnum']))
                    Teff.append(float(row[x_axis]))
                    mag0.append(float(row[y_axis[0]+'_mag']))
                    mag1.append(float(row[y_axis[1]+'_mag']))

colour = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]

# Average colours based on Teff

T = np.arange(min(Teff), max(Teff)+100, 100)
avg_cols = [[] for _ in range(len(T))]

for i in range(len(Teff)):
    for j in range(len(T)):
        if Teff[i] == T[j]:
            avg_cols[j].append(mag0[i] - mag1[i])

avg_col = []
n_col = []
for n in range(len(avg_cols)):
    avg_col.append(np.mean(avg_cols[n]))
    n_col.append(len(avg_cols[n]))

# Variables

x_cif = Teff
y_cif = colour
z_cif = SpT

x_cif_ = Teff_
y_cif_ = colour_

# PLOTTING

# Labels

xlabel = r'$T_{\rm eff}$ [K]'
ylabel = r'$G-G_{RP}$ [mag]'
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

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma_r')

# Plots: distribution

ax.scatter(x_cif_, y_cif_, s=pointsize, marker='o', facecolors='',
           edgecolors='grey', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap, zorder=1)

# Plots: average colours

for n in range(len(avg_col)):
    ax.scatter(T[n], avg_col[n], s=n_col[n],
               facecolors='none', edgecolors='black', lw=1.5)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

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

ax.set_ylim(0.6, 1.9)

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

# Colorbar

plt.savefig(filename+'.png', bbox_inches='tight')
plt.show()
