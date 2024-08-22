# %%

# CALCULATOR: COMPLETENESS
# Cifuentes et al. 2020
# Sample volume completeness

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

Mother_version = '01'
filename = 'cif01_plot_completeness'

SpTmax = 5.0
distances = np.arange(7, 100, 1)
volumes = [4/3*np.pi*distances[i]**3 for i in range(len(distances))]

# Mother

Karmn = [[] for x in range(len(distances))]

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    for row in reader:
        for i in range(len(distances)):
            if row['d_pc'] != '':
                if float(row['d_pc']) <= distances[i] and float(row['SpTnum']) <= SpTmax and float(row['SpTnum']) >= 0.0:
                    if float(row['DE_J2000']) > -23:
                        Karmn[i].append(str(row['Karmn']))

# stars per cubic parsec (assuming complete up to 7 pc)
Nstars = [len(Karmn[i]) for i in range(len(distances))]
Nstars_volume_ref = Nstars[0]/volumes[0]
Nstars_volume = [Nstars[i]/volumes[i] for i in range(len(distances))]

# Variables

x = distances+0.5  # offline
y = Nstars_volume/Nstars_volume_ref*100

# Plotting

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)

# Labels

xlabel = r'$d$ [pc]'
ylabel = r'Completeness'
cbarlabel = r'Spectral type'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (12, 10)
pointsize = 30
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Plots: distribution

ax.bar(x, y, width=1)
plt.axvline(7, color='black', linestyle='dashed', linewidth=2)

# Aesthetics

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

# Axes: range & scale

ax.set_xlim(1, 100)

# Show & Save

plt.savefig(filename+'.png', bbox_inches='tight')
plt.show()

# Print out

print('From M0.0 to', SpTmax, '\nDistance | Completeness')
for i in range(len(distances)):
    print(distances[i], '|', y[i])
