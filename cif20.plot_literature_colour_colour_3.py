# %%

# PLOT: LITERATURE > COLOUR - COLOUR
# g - r vs. r - i
# Cifuentes et al. 2020
# Bochanski et al. 2007
# Davenport et al. 2014

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams.update({'errorbar.capsize': 4})

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_literature_gr_ri'

booleans = ['Bool_delta', 'Bool_young']

x_axis = ['r', 'i']
y_axis = ['g', 'r']

x_axis_alt = 'r-i'
y_axis_alt = 'g-r'
z_axis_alt = 'Nstars'

SpT_range_pre = np.arange(0, 18.5, 0.5)  # M0.0 to L8.0
SpT_range = [-2, -1]  # K5 and K7
for i in range(len(SpT_range_pre)):
    SpT_range.append(np.ndarray.tolist(
        np.arange(0, 18.5, 0.5))[i])  # all range

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    d_pc = []
    mag0 = []
    mag1 = []
    mag2 = []
    mag3 = []
    SpT_ = []
    d_pc_ = []
    mag0_ = []
    mag1_ = []
    mag2_ = []
    mag3_ = []
    mag0_mean = [[] for i in range(len(SpT_range))]
    mag1_mean = [[] for i in range(len(SpT_range))]
    mag2_mean = [[] for i in range(len(SpT_range))]
    mag3_mean = [[] for i in range(len(SpT_range))]
    for row in reader:
        if row['d_pc'] != '':
            if row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
                SpT_.append(float(row['SpTnum']))
                d_pc_.append(float(row['d_pc']))
                mag0_.append(float(row[x_axis[0]+'_mag']))
                mag1_.append(float(row[x_axis[1]+'_mag']))
                mag2_.append(float(row[y_axis[0]+'_mag']))
                mag3_.append(float(row[y_axis[1]+'_mag']))
                if row[booleans[0]] == row[booleans[1]] == 'false':
                    if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                        SpT.append(float(row['SpTnum']))
                        d_pc.append(float(row['d_pc']))
                        mag0.append(float(row[x_axis[0]+'_mag']))
                        mag1.append(float(row[x_axis[1]+'_mag']))
                        mag2.append(float(row[y_axis[0]+'_mag']))
                        mag3.append(float(row[y_axis[1]+'_mag']))
                        for m in range(len(SpT_range)):
                            if float(row['SpTnum']) == SpT_range[m]:  # Choose Spectral type
                                mag0_mean[m].append(
                                    float(row[x_axis[0]+'_mag']))
                                mag1_mean[m].append(
                                    float(row[x_axis[1]+'_mag']))
                                mag2_mean[m].append(
                                    float(row[y_axis[0]+'_mag']))
                                mag3_mean[m].append(
                                    float(row[y_axis[1]+'_mag']))

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_x_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]
colour_y = [mag2[i] - mag3[i] for i in range(len(mag2))]
colour_y_ = [mag2_[i] - mag3_[i] for i in range(len(mag2_))]

colour_x_mean = []
colour_y_mean = []
colour_size = []
for i in range(len(SpT_range)):
    colour_x_mean.append(
        np.mean([mag0_mean[i][j]-mag1_mean[i][j] for j in range(len(mag0_mean[i])-12)]))  # L2 or less
    colour_y_mean.append(
        np.mean([mag2_mean[i][j]-mag3_mean[i][j] for j in range(len(mag2_mean[i])-12)]))  # L2 or less
    colour_size.append(len(mag0_mean[i]))

# Literature
# Introduce the literature source as [['filename.csv', 'alias']]
# Covey07 and Davenport14

papers = [['Literature/davenport14.colours1.csv', 'Dav14']]

colour_x_lit = [[] for _ in range(len(papers))]
colour_y_lit = [[] for _ in range(len(papers))]
nstars_lit = [[] for _ in range(len(papers))]
z_lit = [[] for _ in range(len(papers))]

for i in range(len(papers)):
    with open(papers[i][0], 'r') as mycsv:
        reader = csv.DictReader(mycsv)
        for row in reader:
            if row[x_axis_alt] != '' and row[y_axis_alt] != '':
                colour_x_lit[i].append(float(row[x_axis_alt]))
                colour_y_lit[i].append(float(row[y_axis_alt]))
                nstars_lit[i].append(float(row[z_axis_alt]))
        for j in range(len(nstars_lit[i])):
            z_lit[i].append(nstars_lit[i][j]/100)

# Bochanski et al. 2007
# "Inactive" colours, except from M8 included.
# [M0, ..., L0] in steps of 1.0

gr_boc07 = [1.40, 1.47, 1.60, 1.59, 1.55, 1.57, 1.55, 1.55, 1.80, 1.74, 2.50]
ri_boc07 = [0.65, 0.80, 1.04, 1.28, 1.42, 1.72, 1.99, 2.35, 2.76, 2.83, 2.54]
SpT_boc07 = list(np.arange(0, 11, 1))

# Variables

x_cif = colour_x
y_cif = colour_y
z_cif = SpT

x_cif_ = colour_x_
y_cif_ = colour_y_

x_mean = colour_x_mean
y_mean = colour_y_mean
z_mean = colour_size

# PLOTTING

# Labels

xlabel = r'$r-i$ [mag]'
ylabel = r'$g-r$ [mag]'
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
plot_colors = ['blue']

# Plots: distribution

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, c=z_cif, cmap=cmap, s=pointsize,
                marker='o', label='This work', zorder=0)
for i in range(len(papers)):
    ax.scatter(colour_x_lit[i], colour_y_lit[i], s=z_lit[i],
               facecolors='', edgecolors=plot_colors[i], label=papers[i][1])

# Plots: Colour locus

ax.scatter(x_mean, y_mean, s=z_mean,
           facecolors='black', edgecolors='black')

# Plots: models

ax.scatter(ri_boc07, gr_boc07, s=pointsize/10,
           facecolors='', edgecolors='green')
ax.plot(ri_boc07, gr_boc07, c='green')

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
plt.gca().xaxis.set_major_locator(MaxNLocator(prune='lower'))

# Axes: range & scale

ax.set_xlim(0.0, 3.1)
ax.set_ylim(0.0, 3.5)

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

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
