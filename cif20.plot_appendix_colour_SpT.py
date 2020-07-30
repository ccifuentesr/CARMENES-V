# %%

# PLOT: PLOT > COLOUR - SPECTRAL TYPE
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
# Uncomment to select colour against SpT.

Mother_version = '01'
booleans = ['Bool_delta', 'Bool_young']
x_axis = ['SpT']

SpT_range_pre = np.arange(0, 18.5, 0.5)
SpT_range = [-2, -1]
for i in range(len(SpT_range_pre)):
    SpT_range.append(np.ndarray.tolist(np.arange(0, 18.5, 0.5))[i])

# y_axis = ['NUV', 'RP']
# filename = 'cif20_plot_NUVRP_SpT'

# y_axis = ['B', 'V']
# filename = 'cif20_plot_BV_SpT'

# y_axis = ['GG', 'W3']
# filename = 'cif20_plot_GW3_SpT'

# y_axis = ['r', 'i']
# filename = 'cif20_plot_ri_SpT'

# y_axis = ['RP', 'W1']
# filename = 'cif20_plot_RPW1_SpT'

y_axis = ['J', 'W2']
filename = 'cif20_plot_JW2_SpT'

# Mother
# Appending values and obtaining mean, stdev and sample size.

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    Teff = []
    Teff_ = []
    mag2 = []
    mag2_ = []
    mag3 = []
    mag3_ = []
    for row in reader:
        if row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '':
            if row['Teff'] != '':
                SpT_.append(float(row['SpTnum']))
                Teff_.append(float(row['Teff']))
                mag2_.append(float(row[y_axis[0]+'_mag']))
                mag3_.append(float(row[y_axis[1]+'_mag']))
                if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == 'false':
                    if row[booleans[0]] == row[booleans[1]] == 'false':
                        SpT.append(float(row['SpTnum']))
                        Teff.append(float(row['Teff']))
                        mag2.append(float(row[y_axis[0]+'_mag']))
                        mag3.append(float(row[y_axis[1]+'_mag']))

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

colour_y = [mag2[i] - mag3[i] for i in range(len(mag2))]
colour_y_ = [mag2_[i] - mag3_[i] for i in range(len(mag2_))]

x_cif_ = SpT_
y_cif_ = colour_y_

x_cif = SpT
y_cif = colour_y
z_cif = Teff

x_mean = SpT_num
y_mean = colour_mean
y_mean_err = colour_stdev
z_mean = colour_size

# Statistics

print('STATISTICS\n-----------')
print('Sample size: ', colour_size)
print('Average colours: ', colour_mean)
print('Mean standard deviations: ', colour_stdev)
print('Mean stdev: ', np.round(np.nanmean(colour_stdev), 2))

# PLOTTING

# Labels
# Uncomment when plotting different colour_y

# ylabel = r'$NUV-G_{RP}$ [mag]'
# ylabel = r'$B-V$ [mag]'
# ylabel = r'$G-W3$ [mag]'
# ylabel = r'$r-i$ [mag]'
# ylabel = r'$G_{RP}-W1$ [mag]'
ylabel = r'$J-W2$ [mag]'

xlabel = r'Spectral type'
cbarlabel = r'T$_{\rm eff}$'
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
cmap = plt.get_cmap('coolwarm_r')

# Plots: distribution

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

# Plots: statistics
# Use range(14) to limit to M6V in NUV-RP and B-V

for i in range(len(x_mean)):
    ax.errorbar(x_mean[i], y_mean[i], xerr=0,
                yerr=y_mean_err[i], c='dimgrey', capsize=5)
    ax.scatter(x_mean[i], y_mean[i], s=z_mean[i],
               facecolors='none', edgecolors='dimgrey')

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

# Axes: ticks
# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

# ax.set_xticks(np.arange(-2, len(SpT_name)-1, 1)) #
ax.set_xticks(np.arange(-2, len(SpT_name)-1, 2))
# ax.set_xticklabels(SpT_name)
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
# Uncomment when plotting different colour_y

# plt.ylim(6.5, 13.5)  # NUV-RP
# plt.ylim(0.5, 2.8)  # B-V
# plt.ylim(1.9,10.1) # G-W3
# plt.ylim(-0.1, 3.2)  # r-i
# plt.ylim(1.2,7.0) # RP-W1
plt.ylim(0.4, 3.4)  # J-W2

# Colorbar
# Uncomment to fix the colour range

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax) # Colorbar
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=1200, vmax=4601)  # Keep uncommented to hide the colorbar
# cbar.set_ticks(np.arange(1200, 4601, 300))
# cbar.ax.invert_yaxis()
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
