# %%

# PLOT: LITERATURE > LUMINOSITY (W/METALLICITY)
# Cifuentes et al. 2020

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
filename = 'cif20_literature_Lbol_meta'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    Teff = []
    Teff_meta = []
    Lbol = []
    eLbol = []
    Lbol_meta = []
    eLbol_meta = []
    FeH_meta = []
    eFeH_meta = []
    FeH_lit = []
    eFeH_lit = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
            if row['Teff'] != '' and row['Teff_meta'] != '':
                if row['FeH_lit'] != '':
                    SpT.append(float(row['SpTnum']))
                    Teff.append(float(row['Teff']))
                    Teff_meta.append(float(row['Teff_meta']))
                    Lbol.append(float(row['Lbol']))
                    eLbol.append(float(row['Lberr']))
                    Lbol_meta.append(float(row['Lbol_meta']))
                    eLbol_meta.append(float(row['Lberr_meta']))
                    FeH_meta.append(float(row['Meta_meta']))
                    FeH_lit.append(float(row['FeH_lit']))

# Variables

x = Lbol_meta
y = Lbol
z = FeH_lit

xerr = [np.log10(eLbol_meta[i]) for i in range(len(Lbol_meta))]
yerr = [np.log10(eLbol[i]) for i in range(len(Lbol))]

xp = np.logspace(-3, 0, 100)
# yp = [0 for i in range(len(xp))]
yp = xp


Delta_Lbol = [(Lbol_meta[i]-Lbol[i])/Lbol[i] for i in range(len(Lbol))]
eDelta_Lbol = [np.sqrt((Lbol[i]*eLbol_meta[i]/Lbol_meta[i]) **
                       2 + (-eLbol[i]/Lbol_meta[i])**2) for i in range(len(Lbol))]
# Delta_Lbol = [np.log10(y[i]) - np.log10(x[i]) for i in range(len(Lbol))]

x_res = x
y_res = Delta_Lbol
z_res = z

# PLOTTING

# Labels

xlabel = r'$L$ [$L_{\rm sol}$] $_{\rm BT-Settl}$'
ylabel = r'$L$ [$L_{\rm sol}$] $_{\rm BT-Settl~CIFIST}$'
# ylabel = r'$\Delta L$'
cbarlabel = r'[Fe/H]'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (12, 12)
pointsize = 40
tickssize = 22
labelsize = 22
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('copper')

# Plots: distribution
sc = ax.scatter(x, y, c=z, cmap=cmap, s=pointsize,
                marker='o', zorder=1)

# Loop over each data point to plot
for xr, yr, ex, ey in zip(x, y, xerr, yerr):
    plt.errorbar(xr, yr,  xerr=10**ex, yerr=10**ey, lw=1,
                 capsize=3, color='gainsboro', zorder=0)


# Plots: models
# 1:1 line

plt.plot(xp, yp, 'k--', lw=2, zorder=0)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

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
# Limited in x by Sch19

ax.set_xlim(1.2E-3, 2E-1)
ax.set_ylim(1.2E-3, 2E-1)
# ax.set_ylim(-0.1, 0.1)
ax.set_xscale('log')
ax.set_yscale('log')
# ax.invert_xaxis()
# ax.invert_yaxis()
# ax.xaxis.set_major_locator(MaxNLocator(prune='lower', nbins= 5))

# Colorbar

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)  # Colorbar
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
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

sc0 = ax[0].scatter(x, y, c=z, s=pointsize,
                    cmap=cmap, marker='o', edgecolor='')

sc1 = ax[1].scatter(x_res, y_res, c=z_res, s=pointsize,
                    cmap=cmap, marker='o', edgecolor='', alpha=0.5)

# Loop over each data point to plot
for xr, yr, ex, ey in zip(x, y, xerr, yerr):
    ax[0].errorbar(xr, yr,  xerr=10**ex, yerr=10**ey, lw=1,
                   capsize=3, color='gainsboro', zorder=0)

# Plots: lines

ax[0].plot(xp, yp, color='grey', linestyle='--', lw=2, zorder=0)
ax[1].axhline(y=0, color='grey', linestyle='--', lw=2, zorder=0)
# ax[1].axhline(y=0+0.05, color='grey', linestyle='--', lw=2, zorder=0)
# ax[1].axhline(y=0-0.05, color='grey', linestyle='--', lw=2, zorder=0)

# Aesthetics

# Axes: labels & legend

ax[0].set_ylabel(ylabel, size=labelsize)
ax[1].set_ylabel(r'$\nabla L$', size=labelsize*0.8)
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

ax[0].set_xlim(1.2E-3, 2E-1)
ax[0].set_ylim(1.2E-3, 2E-1)
ax[0].set_xscale('log')
ax[0].set_yscale('log')

# ax[1].set_ylim([0-55E-3, 0+55E-3])
ax[1].set_ylim([0-12E-2, 0+12E-2])

# Colorbar

# cbar = fig.colorbar(sc0, ax=ax.ravel().tolist(), pad=0.01, aspect=50)
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=45)
sc0.set_clim(vmin=-1, vmax=0.6)
sc1.set_clim(vmin=-1, vmax=0.6)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
