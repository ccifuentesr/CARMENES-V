# %%

# PLOT: Plx_A ~ Plx_B FOR BINARIES
# Cifuentes et al. 2020

import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_binaries_distance'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

y_axis = 'parallax'
x_axis = 'Plx'

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    plx_A_ = []
    eplx_A_ = []
    plx_B_ = []
    eplx_B_ = []
    plx_A = []
    eplx_A = []
    plx_B = []
    eplx_B = []
    for row in reader:
        if row['Binarity'] == 'Known':
            if row['parallax_1'] != '' and row['parallax_2'] != '':
                SpT_.append(float(row['SpTnum']))
                plx_A_.append(float(row['parallax_1']))
                eplx_A_.append(float(row['parallax_error_1']))
                plx_B_.append(float(row['parallax_2']))
                eplx_B_.append(float(row['parallax_error_2']))
        if row['Binarity'] == 'New':
            if row['parallax_1'] != '' and row['parallax_2'] != '':
                SpT.append(float(row['SpTnum']))
                plx_A.append(float(row['parallax_1']))
                eplx_A.append(float(row['parallax_error_1']))
                plx_B.append(float(row['parallax_2']))
                eplx_B.append(float(row['parallax_error_2']))


# Variables

x_cif = plx_A
xerr_cif = eplx_A
y_cif = plx_B
yerr_cif = eplx_B
z_cif = SpT

x_cif_ = plx_A_
xerr_cif_ = eplx_A_
y_cif_ = plx_B_
yerr_cif_ = eplx_B_

x_cif_res = plx_A
y_cif_res = [(plx_B[i]-plx_A[i])/plx_A[i] for i in range(len(plx_A))]
z_cif_res = SpT

x_cif_res_ = plx_A_
y_cif_res_ = [(plx_B_[i]-plx_A_[i])/plx_A_[i] for i in range(len(plx_A_))]

xp = np.linspace(0, 100, 1000)

# PLOTTING

# Labels

xlabel = r'$\varpi_A$ [mas]'
ylabel = r'$\varpi_B$ [mas]'
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

sc = ax.scatter(x_cif, y_cif, c=z_cif, s=pointsize,
                cmap=cmap, marker='o', edgecolor='k')
sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize,
                 facecolors='none', edgecolors='grey', zorder=0)

# Plots: models

ax.plot(xp, xp, '--', c='grey', lw=2, zorder=0)
ax.plot(xp, xp*1.05, '--', c='grey', lw=1.5, zorder=0)
ax.plot(xp, xp*0.95, '--', c='grey', lw=1.5, zorder=0)

# Aesthetics

# Axes: range & scale

ax.set_xlim(5, 100)
ax.set_ylim(5, 100)
ax.set_xscale('log')
ax.set_yscale('log')
# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
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
ax.set_xticks((10, 100))
ax.set_yticks((10, 100))

# Labels: outliers

labels_out = ['PM J15380+3224', 'PM J11585+4626']
coords_out = [[16.5, 15.1], [16.2, 14.5]]  # ,[6.769,6.424]]
coords_out_ = [[23, 14.5], [23, 12.5]]  # ,[7.0,6.0]]

plt.arrow(22.5, 15.1, -5, 0, fc="k", ec="k", head_width=0.5, head_length=0.5)
plt.arrow(22.5, 13, -5.2, 1.1, fc="k", ec="k", head_width=0.5, head_length=0.5)

for i, type in enumerate(labels_out):
    plt.text(coords_out_[i][0], coords_out_[i][1], type, fontsize=18)

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

# %%
# PLOTTING (W/RESIDUALS)

# Canvas & Colours

fig, ax = plt.subplots(2, 1, sharex='col', gridspec_kw={
                       'hspace': 0.1, 'wspace': 0.4, 'height_ratios': [10, 2], 'hspace': 0.03}, figsize=figsize)

# Plot: distribution

sc0 = ax[0].scatter(x_cif, y_cif, c=z_cif, s=pointsize,
                    cmap=cmap, marker='o', edgecolor='k')
sc0_ = ax[0].scatter(x_cif_, y_cif_, s=pointsize,
                     facecolors='none', edgecolors='grey', zorder=0)

sc1 = ax[1].scatter(x_cif_res, y_cif_res, c=z_cif_res, s=pointsize,
                    cmap=cmap, marker='o', edgecolor='k')
sc1_ = ax[1].scatter(x_cif_res_, y_cif_res_, s=pointsize,
                     facecolors='none', edgecolors='grey', zorder=0)

# Plots: lines

ax[0].plot(xp, xp, '--', c='grey', lw=2, zorder=0)
ax[0].plot(xp, xp*1.05, '--', c='grey', lw=1.5, zorder=0)
ax[0].plot(xp, xp*0.95, '--', c='grey', lw=1.5, zorder=0)
ax[1].axhline(y=0, color='grey', linestyle='--', lw=2, zorder=0)
ax[1].axhline(y=0+0.05, color='grey', linestyle='--', lw=2, zorder=0)
ax[1].axhline(y=0-0.05, color='grey', linestyle='--', lw=2, zorder=0)

# Plots: labels

for i, type in enumerate(labels_out):
    ax[0].text(coords_out_[i][0], coords_out_[i][1], type, fontsize=18)

ax[0].arrow(22.5, 15.1, -5, 0, fc="k", ec="k", head_width=0.5, head_length=0.5)
ax[0].arrow(22.5, 13, -5.2, 1.1, fc="k", ec="k",
            head_width=0.5, head_length=0.5)

# Aesthetics

# Plot: empty-faced

# sc1.set_facecolor('none')

# Axes: labels & legend

ax[0].set_ylabel(ylabel, size=labelsize)
ax[1].set_ylabel(r'$\nabla\varpi$', size=labelsize*0.8)
ax[1].set_xlabel(xlabel, size=labelsize)
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# Axes: ticks

ax[0].tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=False)
ax[0].tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax[1].tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both')
ax[1].tick_params(axis='y', labelsize=tickssize*0.8, direction='in',
                  right=True, labelright=False, which='both')
ax[0].tick_params('both', length=10, width=1, which='major')
ax[0].tick_params('both', length=5, width=1, which='minor')
ax[1].tick_params('both', length=10, width=1, which='major')
ax[1].tick_params('both', length=5, width=1, which='minor')
ax[0].xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax[0].minorticks_on()
ax[1].xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Axes: range & scale

ax[0].set_xlim(5, 100)
ax[0].set_ylim(5, 100)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_ylim([0-12E-2, 0+12E-2])

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

cbar = fig.colorbar(sc0, ax=ax.ravel().tolist(), pad=0.01, aspect=50)
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=45)
sc0.set_clim(vmin=-2, vmax=18)
sc1.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
