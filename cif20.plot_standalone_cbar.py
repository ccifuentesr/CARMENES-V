# %%

# PLOT: CUSTOM STANDALONE COLOURBAR
# Cifuentes et al. 2020

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Uncomment for SpT/Teff/[Fe/H] colourbar.

# filename = 'cif20_colorbar_SpT.eps'
filename = 'cif20_colorbar_Teff.eps'
# filename = 'cif20_colorbar_FeH.eps'


a = np.array([[-2, 19]])
SpTypes = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
           'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']

# PLOTTING

# Canvas & Colours
# Uncomment for SpT/Teff/[Fe/H] colourbar.

plt.figure(figsize=(16, 0.5))
# img = plt.imshow(a, cmap="magma_r")  # SpT
img = plt.imshow(a, cmap="coolwarm_r")  # Teff
# img = plt.imshow(a, cmap="copper")  # [Fe/H]
plt.gca().set_visible(False)

# Colourbar
# Uncomment for SpT/Teff/[Fe/H] colourbar.

cax = plt.axes([0.1, 0.2, 0.8, 0.4])
cbar = plt.colorbar(orientation='horizontal', cax=cax)
cbar.ax.tick_params(labelsize=16)  # SpT and Teff
# cbar.ax.tick_params(labelsize=18) # [Fe/H]

# SpT
# cbar.set_ticks(range(-2, 19))
# plt.clim(-2, 18)
# cbar.ax.set_xticklabels(SpTypes)
# cbar.outline.set_visible(False)

# Teff
plt.clim(1300, 4601)
cbar.set_ticks(np.arange(1300, 4601, 300))
cbar.outline.set_visible(False)
plt.gca().invert_xaxis()

# [Fe/H]
# plt.clim(-1.0, 0.6)
# cbar.set_ticks(np.arange(-1.0, 0.61, 0.1))
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename, bbox_inches='tight')
plt.show()
