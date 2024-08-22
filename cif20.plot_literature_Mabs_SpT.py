# %%

# PLOT: LITERATURE > ABSOLUTE MAGNITUDE - SPECTRAL TYPE
# FUV - B vs. B - V
# Cifuentes et al. 2020
# Hawley et al. 2002
# Kiman et al. 2019

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
filename = 'cif01_literature_MJ_SpT'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

y_axis = 'J'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    d_pc = []
    d_pc_ = []
    mag0 = []
    mag0_ = []
    for row in reader:
        if row[y_axis+'_mag'] != '':
            if row['d_pc'] != '':
                SpT_.append(float(row['SpTnum']))
                d_pc_.append(float(row['d_pc']))
                mag0_.append(float(row[y_axis+'_mag']))
                if row['Qf_'+y_axis] == 'false':
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        SpT.append(float(row['SpTnum']))
                        d_pc.append(float(row['d_pc']))
                        mag0.append(float(row[y_axis+'_mag']))


absmag = [mag0[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag0))]
absmag_ = [mag0_[i] - 5*np.log10(d_pc_[i]) + 5 for i in range(len(mag0_))]

# Kiman 2019
# Table 4: Mean Gaia DR2 Absolute Magnitudes and Colors for MLSDSS-GaiaDR2.
# Table 6: Best fit parameters for SDSS and 2MASS magnitudes as a function of spectral type.

MG_kim19 = [[0, 8.13], [1, 8.78], [2, 9.35], [3, 9.97], [4, 10.77], [5, 11.86], [6, 12.92], [7, 13.54], [
    8, 14.6], [9, 15.26], [10, 16.11], [11, 16.82], [12, 17.11], [13, 17.89], [14, 18.51], [16, 18.92]]


def MJ_kim19(SpT):
    C = [-0.007, 0.65, 5.7]
    Cerr = [0.002, 0.04, 0.1]
    MJ = C[0]*SpT**2 + C[1]*SpT + C[2]
    return MJ


xp = np.linspace(-2, 18, 1000)

# Hawley et al. 2002
# Equation 2, page 3423.


def MJ_haw02(x):
    if(-1 <= x <= 3):
        return 6.46 + 0.26*x  # K5 to M3
    if(x == 4):
        return 8.34  # M4
    if(5 <= x <= 7):
        return 5.73 + 0.74*x  # M5 to M7
    if(8 <= x <= 15):
        return 8.83 + 0.29*x  # M8 to L5
    if(16 <= x <= 28):
        return 12.13 + 0.14*x  # L6 to T8

# Variables


x_cif_ = SpT_
y_cif_ = absmag_

x_cif = SpT
y_cif = absmag
z_cif = SpT

x_haw02 = xp
y_haw02 = []
for i in range(len(x_haw02)):
    y_haw02.append(MJ_haw02(x_haw02[i]))

x_kim19 = xp
y_kim19 = MJ_kim19(xp)


# PLOTTING

# Labels

xlabel = r'Spectral type'
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
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

ax.plot(x_kim19, y_kim19, c='magenta', ls='--', lw=2)
ax.plot(x_haw02, y_haw02, c='blue', ls='--', lw=2)
ax.plot(4, 8.34, 'o', c='blue')  # Single value for Hawley+02

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# Axes: ticks
# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

ax.set_xticks(np.arange(-2, len(SpT_name)-1, 2))
ax.set_xticklabels(SpT_half)
ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)

# Axes: range & scale

ax.set_ylim(3.5, 15.5)
ax.invert_yaxis()

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

plt.savefig(filename+'.png', bbox_inches='tight')
plt.show()
