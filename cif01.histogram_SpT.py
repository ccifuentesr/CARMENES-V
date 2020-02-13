#!/usr/bin/env python
# coding: utf-8

# In[49]:


"""
SpT
Histogram
=========
Cifuentes et al. 2019
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
import random

rc('font',**{'family':'serif','serif':['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

"""
Data
----
"""

## Global

Mother_version = '367'

filename = 'cif01_hist_SpT'

with open('cif01.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    SpTnum = []
    for row in Mother:
        try:
            SpTnum.append(float(row['SpTnum']))
        except ValueError:
            next(Mother)

SpTypes = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
           'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']

"""
Plotting
--------
"""

## Labels

xlabel = r'Spectral type'
ylabel = r'Number of stars'

## Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18
bins = len(SpTypes)-1

## Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)

## Plots: histogram

plt.hist(SpTnum, bins=bins, density=False, histtype='step', color='blue', linewidth=2)

"""
Aesthetics
----------
"""

## Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.set_xticklabels(SpTypes)

## Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize-2, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=False, top=False)
plt.xticks(np.arange(-2, 21, 1))
plt.setp(ax.xaxis.get_majorticklabels(), rotation=-45, ha="left")

## Axes: range & scale

plt.xlim(-2,19)
# plt.ylim(,)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.invert_xaxis()
# ax.invert_yaxis()

"""
Show & Save
----------
"""

plt.savefig(filename+'.png', bbox_inches='tight')
plt.show()

# %%
