# %%

# PLOT: LITERATURE > ABSOLUTE MAGNITUDE - COLOUR
# MJ vs. J - Ks
# Cifuentes et al. 2020
# Lepine et al. 2011, 2013
# Knapp et al. 2004

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
filename = 'cif20_literature_MJ_JKs'
booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

y_axis = 'J'
x_axis = ['J', 'K']
x_axis_alt = ['J', 'Ks']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    d_pc = []
    mag0 = []
    mag1 = []
    mag2 = []
    SpT_ = []
    d_pc_ = []
    mag0_ = []
    mag1_ = []
    mag2_ = []
    for row in reader:
        if row['d_pc'] != '':
            if row[y_axis+'_mag'] != '' and row[x_axis_alt[0]+'_mag'] != '' and row[x_axis_alt[1]+'_mag'] != '':
                SpT_.append(float(row['SpTnum']))
                d_pc_.append(float(row['d_pc']))
                mag0_.append(float(row[x_axis_alt[0]+'_mag']))
                mag1_.append(float(row[x_axis_alt[1]+'_mag']))
                mag2_.append(float(row[y_axis+'_mag']))
                if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                    if row['Qf_'+y_axis] == row['Qf_'+x_axis_alt[0]] == row['Qf_'+x_axis_alt[1]] == 'false':
                        SpT.append(float(row['SpTnum']))
                        d_pc.append(float(row['d_pc']))
                        mag0.append(float(row[x_axis_alt[0]+'_mag']))
                        mag1.append(float(row[x_axis_alt[1]+'_mag']))
                        mag2.append(float(row[y_axis+'_mag']))

colour = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]
Mabs = [mag2[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag0))]
Mabs_ = [mag2_[i] - 5*np.log10(d_pc_[i]) + 5 for i in range(len(mag0_))]

# Lepine et al. 2013
# Introduce additional literature sources as [['filename.csv', 'alias']].
# Lepine13 'plx1' is Trigonometric parallax. Name normalised to 'plx'.

papers = [['Literature/lepine13.csv', 'Lep13']]

mag0_lit = [[] for _ in range(len(papers))]
mag1_lit = [[] for _ in range(len(papers))]
mag2_lit = [[] for _ in range(len(papers))]
plx_lit = [[] for _ in range(len(papers))]
colour_lit = [[] for _ in range(len(papers))]
Mabs_lit = [[] for _ in range(len(papers))]

for i in range(len(papers)):
    with open(papers[i][0], 'r') as mycsv:
        reader = csv.DictReader(mycsv)
        for row in reader:
            if row[x_axis[0]+'mag'] != '' and row[x_axis[1]+'mag'] != '' and row[y_axis+'mag'] != '':
                if row['plx'] != '':
                    mag0_lit[i].append(float(row[x_axis[0]+'mag']))
                    mag1_lit[i].append(float(row[x_axis[1]+'mag']))
                    mag2_lit[i].append(float(row[y_axis+'mag']))
                    plx_lit[i].append(float(row['plx']))
        colour_lit[i].append([mag0_lit[i][j] - mag1_lit[i][j]
                              for j in range(len(mag0_lit[i]))])
        Mabs_lit[i].append([mag2_lit[i][j] - 5*np.log10(1 /
                                                        plx_lit[i][j]) + 5 for j in range(len(mag0_lit[i]))])

# Knapp et al. 2004
# 'photometry_c_L' has had T objects removed, only L left.

with open('Literature/knapp04.photometry_c_L.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Mabs_Kna04 = []
    colour_Kna04 = []
    for row in reader:
        if row['JMAG'] != '':
            Mabs_Kna04.append(float(row['JMAG']))
            colour_Kna04.append(float(row['J-K']))

# Variables

x_cif = colour
y_cif = Mabs
z_cif = SpT

x_cif_ = colour_
y_cif_ = Mabs_

x_Kna04 = colour_Kna04
y_Kna04 = Mabs_Kna04

# PLOTTING

# Labels

xlabel = r'$J-Ks$ [mag]'
ylabel = r'M$_J$ [mag]'
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
    ax.scatter(colour_lit[i], Mabs_lit[i], s=pointsize,
               facecolors='', edgecolors=plot_colors[i], label=papers[i][1])

ax.scatter(x_Kna04, y_Kna04, s=pointsize, marker='o',
           facecolors='', edgecolors='green', label='Kna04')

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

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

ax.set_xlim(0.5, 2.3)
ax.set_ylim(3.5, 15.5)
ax.invert_yaxis()

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax)  # Colorbar
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)  # Keep only this uncomment to hide colorbar.
# cbar.set_ticks(np.arange(-2, 19, 2))
# cbar.ax.set_yticklabels(SpT_half)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
