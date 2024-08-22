
# %%

# PLOT: GAIA DR2 BP-RP EXCESS
# Cifuentes et al. 2020
# See also Evans et al. 2018

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
filename = 'cif01_plot_colour_excess'

booleans = ['Bool_delta', 'Bool_young',
            'Bool_dphot', 'Bool_RUWE', 'Bool_excess']

x_axis = ['BP', 'RP']
y_axis = 'phot_bp_rp_excess_factor_1'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    d_pc = []
    d_pc_ = []
    mag0 = []
    mag0_ = []
    mag1 = []
    mag1_ = []
    bp_rp_excess = []
    bp_rp_excess_ = []
    for row in reader:
        if row[y_axis] != '':
            if row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
                if row['d_pc'] != '':
                    SpT_.append(float(row['SpTnum']))
                    d_pc_.append(float(row['d_pc']))
                    mag0_.append(float(row[x_axis[0]+'_mag']))
                    mag1_.append(float(row[x_axis[1]+'_mag']))
                    bp_rp_excess_.append(float(row[y_axis]))
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        if row['Qf_'+x_axis[0]] == 'false':
                            SpT.append(float(row['SpTnum']))
                            d_pc.append(float(row['d_pc']))
                            mag0.append(float(row[x_axis[0]+'_mag']))
                            mag1.append(float(row[x_axis[1]+'_mag']))
                            bp_rp_excess.append(float(row[y_axis]))

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_x_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]

# Variables

x_cif_ = colour_x_
y_cif_ = bp_rp_excess_

x_cif = colour_x
y_cif = bp_rp_excess
z_cif = SpT

# Models: see Evans et al. 2018
# If excess_p_rp = C then
# C < f1 AND C > f2

xp = np.linspace(0, 6, 100)


def f1(x):
    f = 1.3 + 0.060 * x**2
    return f


def f2(x):
    f = 1.0 + 0.015 * x**2
    return f


# PLOTTING

# Labels

xlabel = r'$G_{B_P}-G_{R_P}$ [mag]'
ylabel = r'phot_bp_rp_excess_factor'
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

# Plots: Models (Evans+18)

plt.plot(xp, f1(xp), 'b--')
plt.plot(xp, f2(xp), 'b--')

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

plt.xlim(1.1, 5.8)
plt.ylim(1.2, 2.8)


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
