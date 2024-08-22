# %%

# PLOT: LITERATURE > COLOUR - COLOR
# J - H vs. g - i
# Cifuentes et al. 2020
# Covey et al. 2007
# Davenport et al. 2014

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
filename = 'cif01_literature_JH_gi'

booleans = ['Bool_delta', 'Bool_young']

x_axis = ['g', 'i']
y_axis = ['J', 'H']

x_axis_alt = 'g-i'
y_axis_alt = 'J-H'
z_axis_alt = 'Nstars'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    d_pc = []
    mag0 = []
    mag1 = []
    mag2 = []
    mag3 = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == 'false':
            if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                if row['d_pc'] != '':
                    if row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
                        SpT.append(float(row['SpTnum']))
                        d_pc.append(float(row['d_pc']))
                        mag0.append(float(row[x_axis[0]+'_mag']))
                        mag1.append(float(row[x_axis[1]+'_mag']))
                        mag2.append(float(row[y_axis[0]+'_mag']))
                        mag3.append(float(row[y_axis[1]+'_mag']))

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_y = [mag2[i] - mag3[i] for i in range(len(mag2))]

# Literature
# Introduce the literature source as [['filename.csv', 'alias']]
# Covey07 and Davenport14

papers = [['Literature/davenport14.colours1.csv', 'Dav14'],
          ['Literature/covey07.colours.csv', 'Cov07']]

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

# Variables

x_cif = colour_x
y_cif = colour_y
z_cif = SpT

# PLOTTING

# Labels

xlabel = r'$g-i$ [mag]'
ylabel = r'$J-H$ [mag]'
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
plot_colors = ['green', 'blue']

# Plots: distribution

sc = ax.scatter(x_cif, y_cif, c=z_cif, cmap=cmap, s=pointsize,
                marker='o', label='This work', zorder=0)
for i in range(len(papers)):
    ax.scatter(colour_x_lit[i], colour_y_lit[i], s=z_lit[i],
               facecolors='', edgecolors=plot_colors[i], label=papers[i][1])

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
# plt.legend(loc='best', prop={'size': legendsize})

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)

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

plt.savefig(filename+'.png', bbox_inches='tight')
plt.show()
