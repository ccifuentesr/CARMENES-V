# %%

# PLOT: COLOUR - COLOUR
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
# Uncomment to see Gaia/2MASS diagrams.

Mother_version = '01'
filename = 'cif20_plot_BPRP_GJ'  # Gaia
# filename = 'cif20_plot_JH_HKs' # 2MASS

booleans = ['Bool_delta', 'Bool_young']

x_axis = ['GG', 'J']  # Gaia
y_axis = ['BP', 'RP']  # Gaia
# x_axis = ['H', 'Ks']  # 2MASS
# y_axis = ['J', 'H']  # 2MASS

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    mag0 = []
    mag0_ = []
    mag1 = []
    mag1_ = []
    mag2 = []
    mag2_ = []
    mag3 = []
    mag3_ = []
    emag0 = []
    emag1 = []
    emag2 = []
    emag3 = []
    for row in reader:
        if row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
            SpT_.append(float(row['SpTnum']))
            mag0_.append(float(row[x_axis[0]+'_mag']))
            mag1_.append(float(row[x_axis[1]+'_mag']))
            mag2_.append(float(row[y_axis[0]+'_mag']))
            mag3_.append(float(row[y_axis[1]+'_mag']))
            if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                if row[booleans[0]] == row[booleans[1]] == 'false':
                    SpT.append(float(row['SpTnum']))
                    mag0.append(float(row[x_axis[0]+'_mag']))
                    mag1.append(float(row[x_axis[1]+'_mag']))
                    mag2.append(float(row[y_axis[0]+'_mag']))
                    mag3.append(float(row[y_axis[1]+'_mag']))
                    emag0.append(float(row['e'+x_axis[0]+'_mag']))
                    emag1.append(float(row['e'+x_axis[1]+'_mag']))
                    emag2.append(float(row['e'+y_axis[0]+'_mag']))
                    emag3.append(float(row['e'+y_axis[1]+'_mag']))

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_y = [mag2[i] - mag3[i] for i in range(len(mag2))]
ecolour_x = [np.sqrt(emag0[i]**2 + emag1[i]**2) for i in range(len(emag0))]
ecolour_y = [np.sqrt(emag2[i]**2 + emag3[i]**2) for i in range(len(emag2))]
colour_x_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]
colour_y_ = [mag2_[i] - mag3_[i] for i in range(len(mag2_))]

# Variables

x_cif_ = colour_x_
y_cif_ = colour_y_

x_cif = colour_x
y_cif = colour_y
z_cif = SpT

# PLOTTING

# Labels
# Uncomment to see Gaia/2MASS diagrams.

xlabel = r'$G-J$ [mag]'  # Gaia
ylabel = r'$G_{BP}-G_{RP}$ [mag]'  # Gaia
# xlabel = r'$H-K_s$ [mag]' # 2MASS
# ylabel = r'$J-H$ [mag]' # 2MASS
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

sc_ = plt.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                  edgecolors='gainsboro', zorder=0)
sc = plt.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))  # Gaia
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f')) # 2MASS
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
# Uncomment to see Gaia/2MASS diagrams.

plt.xlim(1.3, 5.6)  # Gaia
plt.ylim(1.0, 5.8)  # Gaia
# plt.xlim(0.0, 0.9)  # 2MASS
# plt.ylim(0.3, 1.3)  # 2MASS

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
