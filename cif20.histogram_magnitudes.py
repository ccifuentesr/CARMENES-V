# %%

# PLOT: HISTOGRAM > MAGNITUDES
# Cifuentes et al. 2020

import csv
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_hist_magnitudes'

# Mother
# Appending data for each passband.

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    FUV = []
    eFUV = []
    NUV = []
    eNUV = []
    u = []
    eu = []
    BT = []
    eBT = []
    B = []
    eB = []
    g = []
    eg = []
    BP = []
    eBP = []
    GG = []
    eGG = []
    VT = []
    eVT = []
    V = []
    eV = []
    r = []
    er = []
    ii = []
    eii = []
    RP = []
    eRP = []
    J = []
    eJ = []
    H = []
    eH = []
    K = []
    eK = []
    W1 = []
    eW1 = []
    W2 = []
    eW2 = []
    W3 = []
    eW3 = []
    W4 = []
    eW4 = []
    for row in reader:
        if row['FUV_mag'] != '' and row['eFUV_mag'] != '':
            FUV.append(float(row['FUV_mag']))
            eFUV.append(float(row['eFUV_mag']))
        if row['NUV_mag'] != '' and row['eNUV_mag'] != '':
            NUV.append(float(row['NUV_mag']))
            eNUV.append(float(row['eNUV_mag']))
        if row['u_mag'] != '' and row['eu_mag'] != '':
            u.append(float(row['u_mag']))
            eu.append(float(row['eu_mag']))
        if row['GG_mag'] != '' and row['eGG_mag'] != '':
            GG.append(float(row['GG_mag']))
            eGG.append(float(row['eGG_mag']))
        if row['BT_mag'] != '' and row['eBT_mag'] != '':
            BT.append(float(row['BT_mag']))
            eBT.append(float(row['eBT_mag']))
        if row['B_mag'] != '' and row['eB_mag'] != '':
            B.append(float(row['B_mag']))
            eB.append(float(row['eB_mag']))
        if row['g_mag'] != '' and row['eg_mag'] != '':
            g.append(float(row['g_mag']))
            eg.append(float(row['eg_mag']))
        if row['BP_mag'] != '' and row['eBP_mag'] != '':
            BP.append(float(row['BP_mag']))
            eBP.append(float(row['eBP_mag']))
        if row['VT_mag'] != '' and row['eVT_mag'] != '':
            VT.append(float(row['VT_mag']))
            eVT.append(float(row['eVT_mag']))
        if row['V_mag'] != '' and row['eV_mag'] != '':
            V.append(float(row['V_mag']))
            eV.append(float(row['eV_mag']))
        if row['r_mag'] != '' and row['er_mag'] != '':
            r.append(float(row['r_mag']))
            er.append(float(row['er_mag']))
        if row['i_mag'] != '' and row['ei_mag'] != '':
            ii.append(float(row['i_mag']))
            eii.append(float(row['ei_mag']))
        if row['RP_mag'] != '' and row['eRP_mag'] != '':
            RP.append(float(row['RP_mag']))
            eRP.append(float(row['eRP_mag']))
        if row['J_mag'] != '' and row['eJ_mag'] != '':
            J.append(float(row['J_mag']))
            eJ.append(float(row['eJ_mag']))
        if row['H_mag'] != '' and row['eH_mag'] != '':
            H.append(float(row['H_mag']))
            eH.append(float(row['eH_mag']))
        if row['Ks_mag'] != '' and row['eKs_mag'] != '':
            K.append(float(row['Ks_mag']))
            eK.append(float(row['eKs_mag']))
        if row['W1_mag'] != '' and row['eW1_mag'] != '':
            W1.append(float(row['W1_mag']))
            eW1.append(float(row['eW1_mag']))
        if row['W2_mag'] != '' and row['eW2_mag'] != '':
            W2.append(float(row['W2_mag']))
            eW2.append(float(row['eW2_mag']))
        if row['W3_mag'] != '' and row['eW3_mag'] != '':
            W3.append(float(row['W3_mag']))
            eW3.append(float(row['eW3_mag']))
        if row['W4_mag'] != '' and row['eW4_mag'] != '':
            W4.append(float(row['W4_mag']))
            eW4.append(float(row['eW4_mag']))

# PLOTTING

# Labels

x = [FUV, NUV, u, BT, B, g, BP, VT, V, r, GG, ii, RP, J, H, K, W1, W2, W3, W4]
xlab = ['$FUV$', '$NUV$', '$u$', '$B_T$', '$B$', '$g$', '$G_{BP}$', '$V_T$', '$V$', '$r$', '$G$', '$i$', '$G_{RP}$',
        '$J$', '$H$', '$K$', '$W1$', '$W2$', '$W3$', '$W4$']
cmap = plt.get_cmap('rainbow')

# Sizes

figsize = (5, 15)
tickssize = 22
labelsize = 10
legendsize = 18
cblabsize = 18
linewidth = 3

# Canvas & Colours

fig, ax = plt.subplots(len(x), 1, sharex=True, figsize=figsize)
fig.subplots_adjust(hspace=0)

for i in range(len(x)):
    legend = []
    ax[i].hist(x[i], bins='fd', color=cmap(i/len(x)), label=xlab[i])
    if i in [0, 1, 3, 7]:  # Reduce y-tick value for smaller samples
        ax[i].set_yticks([0, 50])
        ax[i].tick_params(axis='x', labelsize=labelsize, direction='in')
        ax[i].tick_params(axis='y', labelsize=labelsize, direction='in')
        legend = ax[i].legend(loc="upper left", markerscale=0, frameon=False)
        plt.setp(legend.get_texts(), color=cmap(i/len(x)))
        for item in legend.legendHandles:  # Remove coloured square
            item.set_visible(False)
    else:
        ax[i].set_yticks([0, 100])
        ax[i].tick_params(axis='x', labelsize=labelsize, direction='in')
        ax[i].tick_params(axis='y', labelsize=labelsize, direction='in')
        legend = ax[i].legend(loc="upper left", markerscale=0, frameon=False)
        plt.setp(legend.get_texts(), color=cmap(i/len(x)))
        for item in legend.legendHandles:  # Remove coloured square
            item.set_visible(False)

# Aesthetics

# Axes: range & scale

plt.xlim(min(x[-1]), max(x[0]))

# Axes: labels & legend

plt.xlabel(r'Magnitude [mag]', size='10')
plt.rc('axes', linewidth=.1)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
