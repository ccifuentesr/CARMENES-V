# %%

# PLOT: PLOT > COLOUR - COLOUR
# For Appendix 6-figure plots in
# Cifuentes et al. 2020

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import scipy.stats as stats
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
import pyperclip
import math

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.
# Uncomment to select colour against G - J.

Mother_version = '01'
booleans = ['Bool_delta', 'Bool_young']
x_axis = ['GG', 'J']

# y_axis = ['NUV', 'RP']
# filename = 'cif20_plot_NUVRP_GJ'

# y_axis = ['B', 'V']
# filename = 'cif20_plot_BV_GJ'

# y_axis = ['BP', 'RP']
# filename = 'cif20_plot_BPRP_GJ'

# y_axis = ['g', 'i']
# filename = 'cif20_plot_gi_GJ'

# y_axis = ['r', 'i']
# filename = 'cif20_plot_ri_GJ'

# y_axis = ['J', 'H']
# filename = 'cif20_plot_JH_GJ'

y_axis = ['J', 'W2']
filename = 'cif20_plot_JW2_GJ'

# y_axis = ['W1', 'W3']
# filename = 'cif20_plot_W1W3_GJ'


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

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_y = [mag2[i] - mag3[i] for i in range(len(mag2))]
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
# Uncomment when plotting different colour_y

# ylabel = r'$NUV-G_{RP}$ [mag]'
# ylabel = r'$B-V$ [mag]'
# ylabel = r'$G_{BP}-G_{RP}$ [mag]'
# ylabel = r'$g-i$ [mag]'
# ylabel = r'$r-i$ [mag]'
# ylabel = r'$J-H$ [mag]'
ylabel = r'$J-W2$ [mag]'
# ylabel = r'$W1-W3$ [mag]'

xlabel = r'$G-J$ [mag]'
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
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
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
# Uncomment for the corresponding colour

# plt.ylim(6.0, 13.5)  # NUV-BP
# plt.ylim(0.5,2.8) # B-V
# Automatic for BP-RP
# plt.ylim(0.5, 6.2)  # g-i
# plt.ylim(-0.1, 3.2) # r-i
# plt.ylim(0.3, 1.3) # J-H
plt.ylim(0.4, 3.6)  # J-W2
# plt.ylim(-0.4, 1.9)  # W1-W3

plt.xlim(1.2, 5.7)  # G-J (ALL)

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax) # Colorbar
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)  # Keep uncommented to hide the colorbar
# cbar.set_ticks(np.arange(-2, 19, 2)) #
# cbar.ax.set_yticklabels(SpT_half) #
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
