# %%

# PLOT: LITERATURE > BOLOMETRIC LUMINOSITIES
# Cifuentes et al. 2020
# Pecaut & Mamajek 2013
# Newton et al. 2015
# Faherty et al. 2016

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
filename = 'cif20_literature_Lbol_Teff'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

x_axis = 'Teff'
y_axis = 'Lbol'
y_axis_alt = 'logL'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Lbol = []
    Teff = []
    SpT = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
            if row['Lbol'] != '' and row['Teff'] != '':
                Lbol.append(float(row['Lbol']))
                Teff.append(float(row['Teff']))
                SpT.append(float(row['SpTnum']))


# Newton et al. 2015, Faherty et al. 2016
# Faherty et al. 2016 has 'Lbol' but is logarithmic, changed to logL.
# Newton et al. 2015: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJ/800/85.

papers_alt = [['Literature/newton15_2.csv', 'New15_2'],
              ['Literature/faherty16.young2.csv', 'Fah16']]

Lbol_lit_log = [[] for _ in range(len(papers_alt))]
Lbol_lit_alt = [[] for _ in range(len(papers_alt))]
Teff_lit_alt = [[] for _ in range(len(papers_alt))]

for i in range(len(papers_alt)):
    with open(papers_alt[i][0], 'r') as mycsv:
        reader = csv.DictReader(mycsv)
        for row in reader:
            if row[x_axis] != '' and row[y_axis_alt] != '' and row[y_axis_alt] != '...':
                Teff_lit_alt[i].append(float(row[x_axis]))
                Lbol_lit_log[i].append(float(row[y_axis_alt]))
    for j in range(len(Lbol_lit_log[i])):
        Lbol_lit_alt[i].append(10**Lbol_lit_log[i][j])


# Pecaut & Mamajek 2013
# http://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/ApJS/208/9/table9

with open('Literature/pec13_SED.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Lbol1_pec13 = []
    Lbol2_pec13 = []
    Lbol1_log_pec13 = []
    Lbol2_log_pec13 = []
    Teff1_pec13 = []
    Teff2_pec13 = []
    for row in reader:
        if row['logL1'] != '':
            Teff1_pec13.append(float(row['Teff1']))
            Lbol1_log_pec13.append(float(row['logL1']))
        if row['logL2'] != '':
            Teff2_pec13.append(float(row['Teff2']))
            Lbol2_log_pec13.append(float(row['logL2']))
    for i in range(len(Lbol1_log_pec13)):
        Lbol1_pec13.append(10**Lbol1_log_pec13[i])
    for i in range(len(Lbol2_log_pec13)):
        Lbol2_pec13.append(10**Lbol2_log_pec13[i])

# Variables

x_cif = Teff
y_cif = Lbol
z_cif = SpT

x_bonus = [Teff1_pec13]
y_bonus = [Lbol1_pec13]

# PLOTTING

# Labels

xlabel = r'$T_{\rm eff}$ [K]'
ylabel = r'$L$ [L$_\odot$]'
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
plot_colors = ['red']
plot_colors_alt = ['blue', 'magenta']
bonus_colours = ['green']

# Plots: distribution

sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif,
                cmap=cmap, zorder=0)

for i in range(len(papers_alt)):
    ax.scatter(Teff_lit_alt[i], Lbol_lit_alt[i], s=pointsize,
               facecolors='', edgecolors=plot_colors_alt[i], label=papers_alt[i][1], zorder=i+1, alpha=1)

for i in range(len(x_bonus)):
    ax.scatter(x_bonus[i], y_bonus[i], s=pointsize, facecolors='',
               edgecolors=bonus_colours[i], label='Pec13', alpha=1)

# Aesthetics

# Axes: range & scale

ax.set_xlim(900, 5000)
ax.set_ylim(1E-5, 1)
ax.set_yscale('log')
ax.invert_xaxis()

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)

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

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
