# %%

# PLOT: COLOUR - SPT (W/META)
# Cifuentes et al. 2020
# G - J vs. SpT

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
filename = 'cif20_plot_GJ_SpT_meta'

booleans = ['Bool_delta', 'Bool_young']

x_axis = 'SpTnum'
y_axis = ['GG', 'J']


# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    Teff = []
    Teff_ = []
    mag0 = []
    mag0_ = []
    mag1 = []
    mag1_ = []
    FeH_lit = []
    eFeH_lit = []
    for row in reader:
        if row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '':
            if row['Teff'] != '':
                SpT_.append(float(row['SpTnum']))
                Teff_.append(float(row['Teff']))
                mag0_.append(float(row[y_axis[0]+'_mag']))
                mag1_.append(float(row[y_axis[1]+'_mag']))
                if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == 'false':
                    if row[booleans[0]] == row[booleans[1]] == 'false':
                        if row['FeH_lit'] != '':
                            SpT.append(float(row['SpTnum']))
                            Teff.append(float(row['Teff']))
                            mag0.append(float(row[y_axis[0]+'_mag']))
                            mag1.append(float(row[y_axis[1]+'_mag']))
                            FeH_lit.append(float(row['FeH_lit']))
                            eFeH_lit.append(float(row['eFeH_lit']))

colour_y = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_y_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]

# Average colours

SpT_num_bis = np.arange(0, 18.5, 0.5)
SpT_num = [-2.0, -1.0]
for i in range(len(SpT_num_bis)):
    SpT_num.append(np.ndarray.tolist(np.arange(0, 18.5, 0.5))[i])

avg_mag = [[] for _ in range(len(SpT_num))]  # a column for each filter

colour_mean = []
colour_stdev = []
colour_size = []

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    for row in Mother:
        for n in range(len(SpT_num)):
            if row['SpTnum'] == str(SpT_num[n]):
                if row[booleans[0]] == row[booleans[1]] == 'false':
                    if row['Qf_'+y_axis[0]] == 'false' and row['Qf_'+y_axis[1]] == 'false':
                        avg_mag[n].append(
                            float(row[y_axis[0]+'_mag'])-float(row[y_axis[1]+'_mag']))

for m in range(len(SpT_num)):
    colour_mean.append(round(np.mean(avg_mag[m]), 2))
    colour_stdev.append(round(np.std(avg_mag[m], ddof=1), 2))
    colour_size.append(len(avg_mag[m]))

# Variables

x_cif = SpT
y_cif = colour_y
z_cif = FeH_lit

x_cif_ = SpT_
y_cif_ = colour_y_
z_cif_ = Teff_

x_mean = SpT_num
y_mean = colour_mean
z_mean = colour_size
xerr = 0
yerr = [colour_stdev[i] for i in range(len(colour_stdev))]

# PLOTTING

# Labels

xlabel = r'Spectral type'
ylabel = r'$G-J$[mag]'

cbarlabel = r'[Fe/H]'
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
cmap = plt.get_cmap('copper')

# Plots: distribution

ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='none',
           edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, c=z_cif, cmap=cmap,
                s=pointsize, marker='o', label='This work')

# Plots: mean values and stdev

for i in range(len(colour_stdev)-8):  # limit to L4
    ax.scatter(x_mean[i], y_mean[i], s=z_mean[i],
               facecolors='none', edgecolors='dimgrey')
    ax.errorbar(x_mean[i], y_mean[i], xerr=0,
                yerr=yerr[i], c='dimgrey', capsize=5)

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

ax.set_ylim(1.2, 5.8)  # GJ

# Colorbar

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax)  # Colorbar
sc.set_clim(vmin=-1.0, vmax=0.6)  # Do not comment
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
