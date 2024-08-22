# %%

# PLOT: BC - COLOUR
# Cifuentes et al. 2020
# Plot all BC against G-J colour in a single plot.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from scipy import stats
from scipy import optimize
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
filename = 'cif20_plot_BC_GJ'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']


x_axis = ['GG', 'J']  # filters[10] and filters[13] in the list below
y_axis = ['GG', 'J']

filters = ['FUV', 'NUV', 'u', 'BT', 'B', 'g', 'BP', 'VT', 'V',
           'r', 'GG', 'i', 'RP', 'J', 'H', 'Ks', 'W1', 'W2', 'W3', 'W4']
filters_lab = ['$FUV$', '$NUV$', '$u$', '$B_T$', '$B$', '$g$',
               '$G_{BP}$', '$V_T$', '$V$', '$r$', '$G$', '$i$', '$G_{RP}$', '$J$', '$H$', '$Ks$', '$W1$', '$W2$', '$W3$', '$W4$']
Lsun = 3.828E26  # in Watts (IAU B2 Resolution)

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = [[] for x in range(len(filters))]
    d_pc = [[] for x in range(len(filters))]
    Lbol = [[] for x in range(len(filters))]
    mags = [[] for x in range(len(filters))]
    mag0 = [[] for x in range(len(filters))]
    mag1 = [[] for x in range(len(filters))]
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == 'false':
            if row['d_pc'] != '' and row['Lbol'] != '':
                if row['GG_mag'] != '' and row['J_mag'] != '':
                    if row['Qf_J'] == 'false' and row['Qf_GG'] == 'false':
                        for j in range(len(filters)):
                            if row[filters[j]+'_mag'] != '' and row['Qf_'+filters[j]] == 'false':
                                SpT[j].append(float(row['SpTnum']))
                                d_pc[j].append(float(row['d_pc']))
                                Lbol[j].append(float(row['Lbol']))
                                mags[j].append(float(row[filters[j]+'_mag']))
                                mag0[j].append(float(row[x_axis[0]+'_mag']))
                                mag1[j].append(float(row[x_axis[1]+'_mag']))

colour = []
Mabs = []
Mbol = []
BC = []

for j in range(len(filters)):
    colour.append([mag0[j][i] - mag1[j][i] for i in range(len(mag0[j]))])
    Mabs.append([mags[j][i] - 5*np.log10(d_pc[j][i]) +
                 5 for i in range(len(mags[j]))])
    Mbol.append([71.197425-2.5*np.log10(Lbol[j][i]*Lsun)
                 for i in range(len(Lbol[j]))])  # (IAU B2 Resolution)
    BC.append([np.array(- Mabs[j][i]) + np.array(Mbol[j][i])
               for i in range(len(Mabs[j]))])

# PLOTTING

# Labels

xlabel = r'$G-J$ [mag]'
ylabel = r'BC$_\lambda$ [mag]'
cbarlabel = r'Spectral type'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (12, 18)
pointsize = 20
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('rainbow')

# Plots: distribution
# Avoid warning using [] in c argument.

for i in range(len(filters)):
    ax.scatter(colour[i], BC[i], s=pointsize, facecolors=[cmap(i/len(filters))],
               edgecolors='none', label=filters_lab[i])

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# for i in xrange(5):
#     ax.plot(x, i * x, label='$y = %ix$'%i)

# Shrink current axis by 20%
box = ax.get_position()
handles, labels = ax.get_legend_handles_labels()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis (reversed)
ax.legend(reversed(handles), reversed(labels), handletextpad=0.01, loc='center left',
          bbox_to_anchor=(.98, 0.5), fontsize=16, frameon=False)

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

plt.ylim(-15.5, 7.0)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
