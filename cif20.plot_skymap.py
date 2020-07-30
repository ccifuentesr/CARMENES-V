# %%

# PLOT: SKY MAP (EQUATORIAL/GALACTIC)
# Cifuentes et al. 2020

from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import SkyCoord
from astropy.io import ascii
import astropy.units as u
import astropy.coordinates as coord
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import csv
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Choose between equatorial and galactic coordinates

Mother_version = '01'

# filename = 'cif20_plot_equatorial.eps'
filename = 'cif20_plot_galactic.eps'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    RAJ2000 = []
    DEJ2000 = []
    SpT = []
    for row in reader:
        RAJ2000.append(float(row['RA_J2000']))
        DEJ2000.append(float(row['DE_J2000']))
        SpT.append(float(row['SpTnum']))

# Variables

equatorial = SkyCoord(RAJ2000[:], DEJ2000[:],
                      frame='icrs', unit=(u.degree, u.degree))
galactic = equatorial.galactic

# Equatorial:
xeq = [equatorial[i].ra.wrap_at(180*u.deg) for i in range(len(equatorial))]
yeq = [equatorial[i].dec.wrap_at(180*u.deg) for i in range(len(equatorial))]
xeq_rad = [xeq[i].radian for i in range(len(xeq))]
yeq_rad = [yeq[i].radian for i in range(len(yeq))]

# Galactic:
xga = [galactic[i].l.wrap_at(180*u.deg) for i in range(len(galactic))]
yga = [galactic[i].b.wrap_at(180*u.deg) for i in range(len(galactic))]
xga_rad = [xga[i].radian for i in range(len(xga))]
yga_rad = [yga[i].radian for i in range(len(yga))]

# Uncomment for equatorial

# x = xeq_rad
# y = yeq_rad
# z = SpT

# Uncomment for galactic

x = xga_rad
y = yga_rad
z = SpT

# PLOTTING

# Sizes

figsize = (12, 10)
pointsize = 30
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111, projection="mollweide")
plt.grid(True)

cmap = plt.get_cmap('magma_r')

# Plots: distribution

ax.scatter(x, y, s=pointsize, c=z, cmap=cmap)

# Labels

# Uncomment for equatorial

# xlabel = r'$RA$ [deg]'
# ylabel = r'$DE$ [deg]'

# Uncomment for galactic

xlabel = r'$l$ [deg]'
ylabel = r'$b$ [deg]'

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize*0.8)
ax.tick_params(axis='y', labelsize=tickssize)

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax)  # Colorbar
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)  # Keep uncommented to hide the colorbar
# cbar.set_ticks(np.arange(-2, 19, 2))
# cbar.ax.set_yticklabels(SpT_half)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename, bbox_inches='tight')
plt.show()
