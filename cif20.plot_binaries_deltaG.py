# %%

# PLOT: ∆G ~ rho FOR BINARIES
# Cifuentes et al. 2020

import csv
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rc
from matplotlib.ticker import MaxNLocator

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_binaries_deltag'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

x_axis = 'rho'
y_axis = ['GG_mag', 'phot_g_mean_mag_2']

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    G_A_new = []  # New candidates
    G_B_new = []
    rho_new = []
    SpT_new = []
    G_A_new_ = []  # New candidates without Plx in the secondary
    G_B_new_ = []
    rho_new_ = []
    SpT_new_ = []
    G_A_known = []  # Known binaries
    G_B_known = []
    rho_known = []
    SpT_known = []
    G_A_back = []  # Background sources
    G_B_back = []
    rho_back = []
    SpT_back = []
    G_A_delta = []  # 'Delta' sources
    G_B_delta = []
    rho_delta = []
    for row in reader:
        try:
            if row['Binarity'] == 'New' and row['parallax_2'] != '':
                rho_new.append(float(row[x_axis]))
                SpT_new.append(float(row['SpTnum']))
                G_A_new.append(float(row[y_axis[0]]))
                G_B_new.append(float(row[y_axis[1]]))
            if row['Binarity'] == 'New' and row['parallax_2'] == '':
                rho_new_.append(float(row[x_axis]))
                SpT_new_.append(float(row['SpTnum']))
                G_A_new_.append(float(row[y_axis[0]]))
                G_B_new_.append(float(row[y_axis[1]]))
            if row['Binarity'] == 'Known':
                rho_known.append(float(row[x_axis]))
                SpT_known.append(float(row['SpTnum']))
                G_A_known.append(float(row[y_axis[0]]))
                G_B_known.append(float(row[y_axis[1]]))
            if row['Binarity'] == 'Background':
                if row[y_axis[0]] != '' and row[y_axis[1]] != '':
                    rho_back.append(float(row[x_axis]))
                    SpT_back.append(float(row['SpTnum']))
                    G_A_back.append(float(row[y_axis[0]]))
                    G_B_back.append(float(row[y_axis[1]]))
            if row['Class'] == 'Delta':
                rho_delta.append(float(row['rho']))
                G_A_delta.append(float(row[y_axis[0]]))
                G_B_delta.append(float(row[y_axis[1]]))
        except:
            pass

# Variables

G_AB_new = [abs(G_A_new[i] - G_B_new[i]) for i in range(len(G_A_new))]
G_AB_new_ = [abs(G_A_new_[i] - G_B_new_[i]) for i in range(len(G_A_new_))]
G_AB_known = [abs(G_A_known[i] - G_B_known[i]) for i in range(len(G_A_known))]
G_AB_back = [abs(G_A_back[i] - G_B_back[i]) for i in range(len(G_A_back))]
G_AB_delta = [abs(G_A_delta[i] - G_B_delta[i]) for i in range(len(G_A_delta))]

x0 = rho_new
y0 = G_AB_new
z0 = SpT_new
x1 = rho_new_
y1 = G_AB_new_
z1 = SpT_new_
x2 = rho_known
y2 = G_AB_known
z2 = SpT_known
x3 = rho_back
y3 = G_AB_back
z3 = SpT_back

xp = np.linspace(0, 5, 100)

# PLOTTING

# Labels

xlabel = r'$\rho$ [arcsec]'
ylabel = r'$|G_A-G_B|$ [mag]'
cbarlabel = r'Spectral type'
SpTypes = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
           'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpTypes_bis = [SpTypes[2*i] for i in range(0, int(round(len(SpTypes)/2, 0))+1)]

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)

# Plots: distribution

plt.scatter(x0, y0, s=pointsize, color='blue', zorder=3, label='New')
plt.scatter(x1, y1, s=pointsize, color='blue',
            facecolor='', zorder=3, label='New?')
plt.scatter(x2, y2, s=pointsize, color='silver', zorder=2, label='Known')
plt.scatter(x3, y3, s=pointsize, color='silver',
            marker='x', zorder=1, label='Background')

# Plots: Hortizontal lines & filling

plt.axhline(y=0, color='grey', linestyle='--', linewidth=2.5)
plt.axhline(y=5, color='red', linestyle='--', linewidth=2.5)
ax.fill_betweenx(xp, 5, color='gainsboro', zorder=0)

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
ax.xaxis.set_major_locator(MaxNLocator(prune='lower', nbins=5))

# Axes: range & scale

plt.xlim(0, 5)
plt.ylim(0, 10)
ax.invert_yaxis()

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
