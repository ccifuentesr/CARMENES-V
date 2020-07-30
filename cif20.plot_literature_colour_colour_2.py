# %%

# PLOT: LITERATURE > COLOUR - COLOR
# J - H vs. g - i
# Cifuentes et al. 2020
# Ansdell et al. 2015

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

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif01_plot_NUVKs_VJ'

booleans = ['Bool_delta', 'Bool_young']

x_axis = ['V', 'J']
y_axis = ['NUV', 'Ks']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Karmn = []
    Karmn_ = []
    Karmn_young = []
    SpT = []
    SpT_ = []
    SpT_young = []
    mag0 = []
    mag0_ = []
    mag0_young = []
    mag1 = []
    mag1_ = []
    mag1_young = []
    mag2 = []
    mag2_ = []
    mag2_young = []
    mag3 = []
    mag3_ = []
    mag3_young = []
    for row in reader:
        if row[y_axis[0]+'_mag'] != '' and row[y_axis[1]+'_mag'] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
            Karmn_.append(str(row['Karmn']))
            SpT_.append(float(row['SpTnum']))
            mag0_.append(float(row[x_axis[0]+'_mag']))
            mag1_.append(float(row[x_axis[1]+'_mag']))
            mag2_.append(float(row[y_axis[0]+'_mag']))
            mag3_.append(float(row[y_axis[1]+'_mag']))
            if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                if row[booleans[0]] == row[booleans[1]] == 'false':
                    Karmn.append(str(row['Karmn']))
                    SpT.append(float(row['SpTnum']))
                    mag0.append(float(row[x_axis[0]+'_mag']))
                    mag1.append(float(row[x_axis[1]+'_mag']))
                    mag2.append(float(row[y_axis[0]+'_mag']))
                    mag3.append(float(row[y_axis[1]+'_mag']))
            if row['Qf_'+y_axis[0]] == row['Qf_'+y_axis[1]] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                if row[booleans[1]] == 'true':
                    Karmn_young.append(str(row['Karmn']))
                    SpT_young.append(float(row['SpTnum']))
                    mag0_young.append(float(row[x_axis[0]+'_mag']))
                    mag1_young.append(float(row[x_axis[1]+'_mag']))
                    mag2_young.append(float(row[y_axis[0]+'_mag']))
                    mag3_young.append(float(row[y_axis[1]+'_mag']))

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_y = [mag2[i] - mag3[i] for i in range(len(mag2))]
colour_x_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]
colour_y_ = [mag2_[i] - mag3_[i] for i in range(len(mag2_))]
colour_x_young = [mag0_young[i] - mag1_young[i]
                  for i in range(len(mag0_young))]
colour_y_young = [mag2_young[i] - mag3_young[i]
                  for i in range(len(mag2_young))]

# Ansdell et al. 2015
# Equation 1: NUV − Ks = 7.72 + 1.66(V − J)
# "Basal NUV locus"


def ansdell_locus(x):
    return (7.72 + 1.66*x)

# Variables


x_cif = colour_x
y_cif = colour_y
z_cif = SpT

x_cif_ = colour_x_
y_cif_ = colour_y_
z_cif_ = SpT_

x_cif_young = colour_x_young
y_cif_young = colour_y_young
z_cif_young = SpT_young

xp = np.linspace(0, 6, 100)

# Shaded region
avg_col = 4.15  # Average V-J for M4V
x_shadow_a = np.linspace(0, avg_col, 100)
x_shadow_b = np.linspace(avg_col, 10, 100)

# Young stars
# According to Ansdell et al. 2015

x_young = []
y_young = []
z_young = []
karmn_young_index = []

for i in range(len(x_cif)):
    try:
        x_young.append(float(np.where(x_cif[i] < avg_col and y_cif[i] < (
            ansdell_locus(x_cif[i])-1.12), x_cif[i], '')))
        y_young.append(float(np.where(x_cif[i] < avg_col and y_cif[i] < (
            ansdell_locus(x_cif[i])-1.12), y_cif[i], '')))
        z_young.append(float(np.where(x_cif[i] < avg_col and y_cif[i] < (
            ansdell_locus(x_cif[i])-1.12), z_cif[i], '')))
        karmn_young_index.append(Karmn.index((np.where(
            x_cif[i] < avg_col and y_cif[i] < (ansdell_locus(x_cif[i])-1.12), Karmn[i], ''))))
    except Exception as e:
        pass

# Not Young stars
# Difference of samples

indices = list(range(0, len(x_cif)))
karmn_noyoung_index = list(set(indices) - set(karmn_young_index))
x_noyoung = []
y_noyoung = []
z_noyoung = []

for i in karmn_noyoung_index:
    x_noyoung.append(x_cif[i])
    y_noyoung.append(y_cif[i])
    z_noyoung.append(z_cif[i])

# PLOTTING

# Labels

xlabel = r'$V-J$ [mag]'
ylabel = r'$NUV-Ks$ [mag]'
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

sc_ = ax.scatter(x_noyoung, y_noyoung, s=pointsize,
                 c=z_noyoung, cmap=cmap, zorder=0)
# sc_.set_facecolor('none')
ax.scatter(x_cif_young, y_cif_young, s=pointsize*2,
           marker='*', facecolors='black', edgecolors='black')
ax.scatter(x_young, y_young, s=pointsize, facecolors='none', edgecolors='blue')

# Plots: Lines

ax.plot(xp, ansdell_locus(xp), 'k-.')
ax.plot(xp, ansdell_locus(xp)-1.12, 'k--', alpha=1)
ax.plot(xp, ansdell_locus(xp)+1.12, 'k--', alpha=1)

# Shaded regions

ax.fill_betweenx(xp, ansdell_locus(xp), color='grey', alpha='0.15')
ax.fill_between(x_shadow_b, 9, 20, alpha=0.4, color='gainsboro')
ax.fill_between(x_shadow_a, 7.72 + 1.66*x_shadow_a -
                1.12, 20, alpha=0.4, color='gainsboro')
plt.axvline(x=avg_col, color='grey', linestyle='--')

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
ax.text(avg_col-0.32, 15.85, r'M3$\,$V', fontsize=18, c='grey')

# Axes: range & scale

plt.xlim(1.8, 5.8)
plt.ylim(9.5, 16)
ax.invert_yaxis()

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc_, cax=cax)  # Colorbar
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc_.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.png', bbox_inches='tight')
plt.show()
