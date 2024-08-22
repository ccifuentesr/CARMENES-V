# %%

# PLOT: HISTOGRAM > PHOTOMETRY COMPLETENESS
# Cifuentes et al. 2020

# Completeness of the photometric sample in each passband.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data

Mother_version = '01'
filename = 'cif20_hist_filters'

# Mother
# Appending two variables for each filter.

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    FUV = []
    tFUV = []
    NUV = []
    tNUV = []
    u = []
    tu = []
    BT = []
    tBT = []
    B = []
    tB = []
    g = []
    tg = []
    BP = []
    tBP = []
    G = []
    tG = []
    VT = []
    tVT = []
    V = []
    tV = []
    r = []
    tr = []
    ii = []
    tii = []
    RP = []
    tRP = []
    J = []
    tJ = []
    H = []
    tH = []
    K = []
    tK = []
    W1 = []
    tW1 = []
    W2 = []
    tW2 = []
    W3 = []
    tW3 = []
    W4 = []
    tW4 = []
    for row in reader:
        if row['Qf_FUV'] != '':
            FUV.append(str(row['Qf_FUV']))
        if row['Qf_NUV'] != '':
            NUV.append(str(row['Qf_NUV']))
        if row['Qf_u'] != '':
            u.append(str(row['Qf_u']))
        if row['Qf_BT'] != '':
            BT.append(str(row['Qf_BT']))
        if row['Qf_B'] != '':
            B.append(str(row['Qf_B']))
        if row['Qf_g'] != '':
            g.append(str(row['Qf_g']))
        if row['Qf_BP'] != '':
            BP.append(str(row['Qf_BP']))
        if row['Qf_VT'] != '':
            VT.append(str(row['Qf_VT']))
        if row['Qf_V'] != '':
            V.append(str(row['Qf_V']))
        if row['Qf_GG'] != '':
            G.append(str(row['Qf_GG']))
        if row['Qf_r'] != '':
            r.append(str(row['Qf_r']))
        if row['Qf_i'] != '':
            ii.append(str(row['Qf_i']))
        if row['Qf_RP'] != '':
            RP.append(str(row['Qf_RP']))
        if row['Qf_J'] != '':
            J.append(str(row['Qf_J']))
        if row['Qf_H'] != '':
            H.append(str(row['Qf_H']))
        if row['Qf_Ks'] != '':
            K.append(str(row['Qf_Ks']))
        if row['Qf_W1'] != '':
            W1.append(str(row['Qf_W1']))
        if row['Qf_W2'] != '':
            W2.append(str(row['Qf_W2']))
        if row['Qf_W3'] != '':
            W3.append(str(row['Qf_W3']))
        if row['Qf_W4'] != '':
            W4.append(str(row['Qf_W4']))

# Variables

x = [FUV, NUV, u, BT, B, g, BP, VT, V, r, G, ii, RP, J, H, K, W1, W2, W3, W4]

false = []
true = []
for i in range(len(x)):
    false.append(x[i].count('false'))
    true.append(x[i].count('true'))

xlab = ['$F$', '$N$', '$u$', '$B_T$', '$B$', '$g$', '$G_{BP}$', '$V_T$', '$V$', '$r$', '$G$', '$i$', '$G_{RP}$',
        '$J$', '$H$', '$K$', '$W1$', '$W2$', '$W3$', '$W4$']
index = np.arange(len(x))    # the x locations for the groups
width = 1.0       # the width of the bars: can also be len(x) sequence

xhline = np.arange(-2, 21, 1)
yhline = [len(J) for i in range(len(xhline))]  # J_mag is the most complete

colorset = []
for i in range(len(x)):
    colorset.append(cm.rainbow(i/len(x)))

# PLOTTING

# Labels

xlabel = r'Passband'
ylabel = r'Number of stars'
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

# Plots: stacked bars

p1 = plt.bar(index, false, width, color=colorset)
p2 = plt.bar(index, true, width, bottom=false, color=colorset, alpha=0.5)

# Horizontal lines

plt.plot(xhline, yhline, '--', color='0.5')


# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize-2, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
plt.rcParams["axes.linewidth"] = 1
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=False, top=False)
ax.xaxis.set_major_locator(MaxNLocator(prune='lower', nbins=5))


# Axes: range & scale

# plt.xlim(0,5)
# plt.ylim(0,10)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.invert_xaxis()
# ax.invert_yaxis()

plt.xticks(index, xlab, size='22')
ax.tick_params(axis='x', labelsize=22, direction='in',
               top=False, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=23, direction='in',
               right=True, labelright=False, which='both')
plt.xlim(-0.5, len(x)-0.5)  # removes blank columns after and before
plt.rc('axes', linewidth=2)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()

# Stats

# Used in VOSA
# Magnitudes markes as `Nofit` in VOSA. Taken from `Nfit` and `Ntot` data delivered by VOSA.
# This is for non-binary subsample, not the total sample of stars.

print('STATISTICS\n----------')

print(sum(true), 'Nofit,', sum(false),
      'Fitted, and', sum(false)+sum(true), 'total.')

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Nfit = []
    Ntot = []
    for row in reader:
        if row['Ntot'] != '':
            Nfit.append(int(row['Nfit']))
            Ntot.append(int(row['Ntot']))

print(sum(Nfit), 'Nfit,', sum(Ntot), 'Ntot,',
      round(sum(Nfit)/sum(Ntot)*100, 2), '% ratio.')
