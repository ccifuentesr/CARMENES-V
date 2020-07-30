# %%

# PLOT: LITERATURE > COLOUR - SPECTRAL TYPE
# r - i vs. SpT
# Cifuentes et al. 2020
# Hawley et al. 2002
# Covey et al. 2007
# Bochanski et al. 2008

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
filename = 'cif20_literature_ri_SpT'

booleans = ['Bool_delta', 'Bool_young']

y_axis = ['r', 'i']
x_axis = 'SpTnum'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    Teff = []
    mag0 = []
    mag1 = []
    SpT_ = []
    Teff_ = []
    mag0_ = []
    mag1_ = []
    for row in reader:
        if row[x_axis] != '' and row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '':
            if row['Teff'] != '':
                SpT_.append(float(row[x_axis]))
                Teff_.append(float(row['Teff']))
                mag0_.append(float(row[y_axis[0]+'_mag']))
                mag1_.append(float(row[y_axis[1]+'_mag']))
                if row[booleans[0]] == row[booleans[1]] == 'false':
                    if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == 'false':
                        SpT.append(float(row[x_axis]))
                        Teff.append(float(row['Teff']))
                        mag0.append(float(row[y_axis[0]+'_mag']))
                        mag1.append(float(row[y_axis[1]+'_mag']))

colour_y_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]
colour_y = [mag0[i] - mag1[i] for i in range(len(mag0))]

# Average colours (Mother)

SpT_num = [-1.0, -2.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5,
           6.0, 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 11.5, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]

avg_mag = [[] for _ in range(len(SpT_num))]
avg_SpT = [[] for _ in range(len(SpT_num))]

n_colour = []
mean_colour = []
mean_SpT = []
stdev = []
items = []

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    for row in Mother:
        for n in range(len(SpT_num)):
            if row['SpTnum'] == str(SpT_num[n]):
                if row[booleans[0]] == row[booleans[1]] == 'false':
                    if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == 'false':
                        if float(row['SpTnum']) < 12.5:
                            avg_mag[n].append(
                                float(row[y_axis[0]+'_mag'])-float(row[y_axis[1]+'_mag']))
                            avg_SpT[n].append(float(row['SpTnum']))


for m in range(len(SpT_num)):
    n_colour.append(len(avg_mag[m]))
    mean_colour.append(round(np.mean(avg_mag[m]), 2))
    mean_SpT.append(round(np.mean(avg_SpT[m]), 2))
    stdev.append(round(np.std(avg_mag[m], ddof=1), 2))
    items.append(len(avg_mag[m]))

# Hawley et al. 2002

paper_haw02 = 'Literature/hawley02.colours.csv'

colour_haw02 = []
ecolour_haw02 = []
SpT_haw02 = []

with open(paper_haw02, 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    for row in reader:
        if row['r-i'] != '' and row[x_axis] != '':
            SpT_haw02.append(float(row[x_axis]))
            colour_haw02.append(float(row['r-i']))
            ecolour_haw02.append(float(row['e_r-i']))

# Covey et al. 2007

paper_cov07 = 'Literature/covey07.SpT_colours_alt.csv'

colour_cov07 = []
SpT_cov07 = []

with open(paper_cov07, 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    for row in reader:
        if row['r-i'] != '' and row[x_axis] != '':
            SpT_cov07.append(float(row[x_axis]))
            colour_cov07.append(float(row['r-i']))

# Bochanski et al. 2007
# Table 6: Template Colors
# Active and Inactive colours
# [M0, ..., L0] in steps of 1.0

SpT_boc07 = list(np.arange(0, 8, 1))

# r-i and
ri_boc07_inactive = [0.65, 0.80, 1.04, 1.28, 1.42, 1.72, 1.99, 2.35]
eri_boc07_inactive = [0.12, 0.11, 0.18, 0.19, 0.16, 0.21, 0.12, 0.18]
ri_boc07_active = [0.68, 0.77, 1.07, 1.30,
                   1.49, 1.74, 1.9, 2.33, 2.76, 2.83, 2.54]
eri_boc07_active = [0.13, 0.16, 0.03, 0.15,
                    0.28, 0.21, 0.11, 0.19, 0.12, 0.07, 0.08]

# Sample size
size_boc07_active = [7, 4, 1, 6, 25, 171, 1132, 400, 16, 5, 4]
size_boc07_inactive = [570, 275, 66, 186, 137, 235, 899, 150]

# West et al. 2008
# Table 1: Colors of Late-Type Dwarfs

SpT_wes08 = list(np.arange(0, 11, 1))
colour_wes08 = [0.66, 0.82, 1.00, 1.21,
                1.46, 1.91, 2.11, 2.50, 2.73, 2.81, 2.49]
ecolour_wes08 = [0.12, 0.18, 0.13, 0.16,
                 0.15, 0.13, 0.14, 0.18, 0.25, 0.21, 0.09]

# Variables

x_cif = SpT
y_cif = colour_y
z_cif = SpT

x_cif_ = SpT_
y_cif_ = colour_y_

dx = 0.1  # offset in SpT

x_haw02 = SpT_haw02
y_haw02 = colour_haw02
ey_haw02 = ecolour_haw02

x_cov07 = [SpT_cov07[i] - dx for i in range(len(SpT_cov07))]
y_cov07 = colour_cov07

x_boc07 = SpT_boc07
y_boc07 = ri_boc07_inactive
ey_boc07 = eri_boc07_inactive

x_wes08 = [SpT_wes08[i] + dx for i in range(len(SpT_wes08))]
y_wes08 = colour_wes08
ey_wes08 = ecolour_wes08

# PLOTTING

# Labels

xlabel = r'Spectral type'
ylabel = r'$r-i$ [mag]'
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

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap, zorder=1)
ax.scatter(mean_SpT, mean_colour, s=n_colour, facecolors='black',
           edgecolors='black')

ax.scatter(x_haw02, y_haw02, s=pointsize, facecolors='',
           edgecolors='blue')
for i in range(len(x_haw02)):
    ax.errorbar(x_haw02[i], y_haw02[i], yerr=ey_haw02[i],
                xerr=0, fmt='o', mfc='none', color='blue', capthick=1)

# ax.scatter(x_cov07, y_cov07, s=pointsize, facecolors='',
#            edgecolors='green')
ax.scatter(x_boc07, y_boc07, s=pointsize, facecolors='',
           edgecolors='green')
for i in range(len(x_boc07)):
    ax.errorbar(x_boc07[i], y_boc07[i], yerr=ey_boc07[i],
                xerr=0, fmt='o', mfc='none', color='green', capthick=1)

ax.scatter(x_wes08, y_wes08, s=pointsize, facecolors='',
           edgecolors='red')
for i in range(len(x_wes08)):
    ax.errorbar(x_wes08[i], y_wes08[i], yerr=ey_wes08[i],
                xerr=0, fmt='o', mfc='none', color='red', capthick=1)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

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

ax.set_xlim(-2.5, 19.5)
ax.set_ylim(0, 3.6)

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT.
# Uncomment only (#) to hide colorbar keeping colour consistency.

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
