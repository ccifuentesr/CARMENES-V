# %%

# PLOT: LITERATURE > EFFECTIVE TEMPERATURE - COLOUR
# Teff vs. V - J
# Cifuentes et al. 2020
# Casagrande et al 2008
# Pecaut & Mamajek 2013

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
filename = 'cif20_literature_Teff_colour'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

y_axis = 'Teff'
x_axis = ['V', 'J']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_ = []
    Teff_ = []
    SpT = []
    Teff = []
    mag0 = []
    mag1 = []
    mag0_ = []
    mag1_ = []
    for row in reader:
        if row[y_axis] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
            SpT_.append(float(row['SpTnum']))
            Teff_.append(float(row[y_axis]))
            mag0_.append(float(row[x_axis[0]+'_mag']))
            mag1_.append(float(row[x_axis[1]+'_mag']))
            if row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                    SpT.append(float(row['SpTnum']))
                    Teff.append(float(row[y_axis]))
                    mag0.append(float(row[x_axis[0]+'_mag']))
                    mag1.append(float(row[x_axis[1]+'_mag']))

colour = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]

# Average colours based on Teff

T = np.arange(min(Teff), max(Teff)+100, 100)
avg_cols = [[] for _ in range(len(T))]

for i in range(len(Teff)):
    for j in range(len(T)):
        if Teff[i] == T[j]:
            avg_cols[j].append(mag0[i] - mag1[i])

avg_col = []
n_col = []
for n in range(len(avg_cols)):
    avg_col.append(np.mean(avg_cols[n]))
    n_col.append(len(avg_cols[n]))

# Literature
# Introduce additional sources as [['filename.csv', 'alias']]
# Passeger+18,19 filter names need to have '_' removed for correct data appending.

papers = [['Literature/casagrande08.csv', 'Cas08']]

mag0_lit = [[] for _ in range(len(papers))]
mag1_lit = [[] for _ in range(len(papers))]
colour_lit = [[] for _ in range(len(papers))]
Teff_lit = [[] for _ in range(len(papers))]

for i in range(len(papers)):
    with open(papers[i][0], 'r') as mycsv:
        reader = csv.DictReader(mycsv)
        for row in reader:
            if row[x_axis[0]+'mag'] != '' and row[x_axis[1]+'mag'] != '' and row[y_axis] != '':
                Teff_lit[i].append(float(row[y_axis]))
                mag0_lit[i].append(float(row[x_axis[0]+'mag']))
                mag1_lit[i].append(float(row[x_axis[1]+'mag']))
    colour_lit[i].append([mag0_lit[i][j] - mag1_lit[i][j]
                          for j in range(len(mag0_lit[i]))])

# Pecaut & Mamajek 2013

with open('Literature/mamajek_privcom.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Teff_mamajek = []
    VK_mamajek = []
    JH_mamajek = []
    HKs_mamajek = []
    for row in reader:
        if float(row['Teff']) > 2000 and float(row['Teff']) <= 4600:
            Teff_mamajek.append(float(row['Teff']))
            VK_mamajek.append(float(row['V-Ks']))
            JH_mamajek.append(float(row['J-H']))
            HKs_mamajek.append(float(row['H-Ks']))

VJ_mamajek = [VK_mamajek[i] - JH_mamajek[i] - HKs_mamajek[i]
              for i in range(len(VK_mamajek))]

# Casagrande et al. 2008
# Equation 6: Teff = a + b*X + c*X**2 + d*X**3 (X = V - J)


def Teff_cas08(x):
    [a, b, c, d] = [0.1926, 0.5738, -0.0726, 0.0042]
    Teff_cas08 = 5040/(a + b*x + c*x**2 + d*x**3)
    return Teff_cas08

# Variables


x_cif = colour
y_cif = Teff
z_cif = SpT

x_cif_ = colour_
y_cif_ = Teff_

x_mamajek = VJ_mamajek
y_mamajek = Teff_mamajek

x_lit_model = np.linspace(2.260, 7.231, 1000)  # Table 3 in Cas08
y_lit_model = Teff_cas08(x_lit_model)

# PLOTTING

# Labels

xlabel = r'$V-J$ [mag]'
ylabel = r'$T_{\rm eff}$ [K]'
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
cmap = plt.get_cmap('magma_r')
plot_colors = ['blue', 'green', 'red', 'orange']

# Plots: distribution

for i in range(len(papers)):
    ax.scatter(colour_lit[i], Teff_lit[i], s=pointsize,
               facecolors='', edgecolors=plot_colors[i], label=papers[i][1], alpha=1, zorder=2)

ax.scatter(x_cif_, y_cif_, s=pointsize, marker='o', facecolors='',
           edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap, zorder=1)

# Plots: literature

ax.scatter(x_mamajek, y_mamajek, s=pointsize/10,
           facecolors='none', edgecolors='green', alpha=1)
ax.plot(x_mamajek, y_mamajek, c='green')

# Plots: average colours

for n in range(len(avg_col)):
    ax.scatter(avg_col[n], T[n], s=n_col[n],
               facecolors='black', edgecolors='black')
ax.scatter(avg_col[3], T[3], s=n_col[3], facecolors='none',
           edgecolors='black', label='This work')

# Aesthetics

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

# Axes: range & scale

ax.set_xlim(1.8, 8.6)

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
