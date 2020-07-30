# %%

# PLOT: LITERATURE > METALLICITY
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
filename = 'cif20_literature_meta_meta'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    FeH_meta = []
    eFeH_meta = []
    FeH_lit = []
    eFeH_lit = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
            if row['Teff'] != '' and row['Teff_meta'] != '':
                if row['FeH_lit'] != '':
                    SpT.append(float(row['SpTnum']))
                    FeH_meta.append(float(row['Meta_meta']))
                    FeH_lit.append(float(row['FeH_lit']))
                    eFeH_lit.append(float(row['eFeH_lit']))

# Mother (Bis)

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Meta_lit15 = []
    Meta_lit10 = []
    Meta_lit05 = []
    Meta_lit00 = []
    Meta_lit003 = []
    Meta_lit005 = []
    Meta_VOSA15 = []
    Meta_VOSA10 = []
    Meta_VOSA05 = []
    Meta_VOSA00 = []
    Meta_VOSA003 = []
    Meta_VOSA005 = []
    Teff_lit15 = []
    Teff_lit10 = []
    Teff_lit05 = []
    Teff_lit00 = []
    Teff_lit003 = []
    Teff_lit005 = []
    for row in reader:
        if row['FeH_lit'] != '':
            if row['Meta_meta'] == '-1.5':
                Teff_lit15.append(float(row['Teff']))
                Meta_lit15.append(float(row['FeH_lit']))
                Meta_VOSA15.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '-1.0':
                Teff_lit10.append(float(row['Teff']))
                Meta_lit10.append(float(row['FeH_lit']))
                Meta_VOSA10.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '-0.5':
                Teff_lit05.append(float(row['Teff']))
                Meta_lit05.append(float(row['FeH_lit']))
                Meta_VOSA05.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '0.0':
                Teff_lit00.append(float(row['Teff']))
                Meta_lit00.append(float(row['FeH_lit']))
                Meta_VOSA00.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '0.3':
                Teff_lit003.append(float(row['Teff']))
                Meta_lit003.append(float(row['FeH_lit']))
                Meta_VOSA003.append(float(row['Meta_meta']))
            if row['Meta_meta'] == '0.5':
                Teff_lit005.append(float(row['Teff']))
                Meta_lit005.append(float(row['FeH_lit']))
                Meta_VOSA005.append(float(row['Meta_meta']))

# Variables

x_all = FeH_meta
y_all = FeH_lit
z_all = FeH_lit

y = [Meta_lit15, Meta_lit10, Meta_lit05, Meta_lit00, Meta_lit003, Meta_lit005]
x = [Meta_VOSA15, Meta_VOSA10, Meta_VOSA05,
     Meta_VOSA00, Meta_VOSA003, Meta_VOSA005]
z = [Teff_lit15, Teff_lit10, Teff_lit05, Teff_lit00, Teff_lit003, Teff_lit005]

xp = np.linspace(-2, 1, 100)
yp = xp

# PLOTTING

# Labels

xlabel = r'[Fe/H]$_{\rm BT-Settl}$'
ylabel = r'[Fe/H]$_{\rm Literature}$'
cbarlabel = r'T$_{\rm eff}$ [K]'

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('coolwarm_r')
cmap_bis = plt.get_cmap('copper')

x_lines = [-1.75, -1.25, -0.75, -0.25, 0.15, 0.40, 0.50]
y_lines = [np.median(y[i]) for i in range(len(y))]

x0 = np.linspace(x_lines[0], x_lines[0+1], 100)
y0 = [np.median(y[0]) for i in range(len(x0))]

xr = [[] for _ in range(len(x_lines))]
yr = [[] for _ in range(len(x_lines))]

for i in range(len(x_lines)-1):
    xr[i].append((np.linspace(x_lines[i], x_lines[i+1], 100)))
    yr[i].append([np.median(y[i]) for _ in range(100)])
    sc = ax.scatter(x[i], y[i], c=z[i], cmap=cmap,
                    s=pointsize*2, marker='o', zorder=1)
    plt.plot(xr[i][0], yr[i][0], c=cmap_bis(i/6), lw=3)
    sc.set_facecolor('none')

# Plots: equality

plt.plot(xp, yp, c='grey', linestyle='--', lw=2, zorder=0)

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
ax.xaxis.set_tick_params(which='minor', bottom=False, top=False)
ax.yaxis.set_tick_params(which='minor', left=False, right=False)
ax.set_xticks((-1.5, -1.0, -0.5, 0.0, 0.3, 0.5))
ax.set_yticks((-1.5, -1.0, -0.5, 0.0, 0.3, 0.5))

# Axes: range & scale

ax.set_xlim(-1.55, 0.65)
ax.set_ylim(-1.55, 0.65)

# Colorbar

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)  # Colorbar
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=2499, vmax=4001)
cbar.set_ticks(np.arange(2500, 4001, 150))
cbar.ax.invert_yaxis()
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
