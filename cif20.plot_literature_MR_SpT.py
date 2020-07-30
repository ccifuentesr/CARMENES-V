# %%

# PLOT: LITERATURE > MASS/RADIUS - SPECTRAL TYPE
# Cifuentes et al. 2020
# Pecaut & Mamajek 2013
# Mann et al. 2015

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
# Uncomment to select mass or radius.

Mother_version = '01'

# Radius
filename = 'cif20_literature_radius_SpT'
y_axis = 'R'

# Mass
# filename = 'cif20_literature_mass_SpT'
# y_axis = 'M'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']
x_axis = 'SpTnum'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    Teff = []
    Lbol = []
    eLbol = []
    for row in reader:
        for n in range(len(booleans)):
            if row[booleans[n]] == 'false':
                if row['Lbol'] != '' and row['Teff'] != '':
                    Lbol.append(float(row['Lbol']))
                    Teff.append(float(row['Teff']))
                    SpT.append(float(row[x_axis]))

# Mann et al. 2015

SpT_man15 = []
Mass_man15 = []

with open('Literature/mann15.stars.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    for row in reader:
        if row[x_axis] != '':
            Mass_man15.append(float(row[y_axis]))
            SpT_man15.append(float(row[x_axis]))

# Pecaut & Mamajek 2013

with open('Literature/mamajek_privcom.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_mam = []
    Radius_mam = []
    Mass_mam = []
    for row in reader:
        SpT_mam.append(float(row['SpTnum']))
        if row['R_Rsun'] != '':
            Radius_mam.append(float(row['R_Rsun']))
        if row['M_Msun'] != '':
            Mass_mam.append(float(row['M_Msun']))


# Radius and Mass


def Radius_SB(Lbol, Teff):
    """Stellar radius from the Stefan–Boltzmann law under the black body approximation.

    Args:
        Lbol (float): Bolometric luminosity in solar units.
        Teff (float): Effective temperature in Kelvin.

    Returns:
        float: Stellar radius in solar units.

    Nominal solar values from the IAU B.3 resolution
    on recommended nominal conversion constants for selected solar and planetary properties:
    https://www.iau.org/static/resolutions/IAU2015_English.pdf

    Nominal solar luminosity: 3.828 x 10+26 W (exact)
    Nominal solar radius: 6.957 x 10+8 m (exact)

    Stefan-Boltzmann constant value from 2018 CODATA recommended values:
    https://physics.nist.gov/cuu/pdf/wall_2018.pdf

    Stefan-Boltzman constant, k: 5.670 374 419 x 10-8 W m-2 K-4 (exact)

    """
    Lsun = 3.828*1e26
    Rsun = 6.957*1e8
    sigma = 5.670374419*1e-8
    R = 1/Rsun * np.sqrt(Lbol * Lsun/(4 * np.pi * sigma * Teff**4))
    return R


def Mass_sch19(Radius):
    """Stellar mass from the empirical relation by Schweitzer et al. 2019 
    (2019A&A...625A..68S), based on masses and radii of eclipsing binaries.

    Args:
        Radius (float): Stellar radius in solar units.

    Returns:
        float: Stellar mass in solar units.

    (See Equation 6 in Schweitzer et al. 2019 and references therein).
    """
    a = -0.0240
    b = 1.055
    a_err = 0.0076
    b_err = 0.017
    M = a + b * Radius
    return M


# Variables
# Uncomment to select mass or radius.

Radius = [Radius_SB(Lbol[i], Teff[i]) for i in range(len(Lbol))]
Mass = [Mass_sch19(Radius[i]) for i in range(len(Radius))]

x_cif = SpT
x_mam = SpT_mam
z_cif = SpT

# Radius
y_cif = Radius
y_mam = Radius_mam

# Mass
# y_cif = Mass
# y_mam = Mass_mam

# PLOTTING

# Labels
# Uncomment to select mass or radius.

ylabel = r'$\mathcal{R}$ [$\mathcal{R}_\odot$]'
# ylabel = r'$\mathcal{M}$ [$\mathcal{M}_\odot$]'

xlabel = r'Spectral type'
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

sc = ax.scatter(x_cif, y_cif, s=pointsize, marker='o',
                c=z_cif, cmap=cmap, zorder=0)

ax.scatter(x_mam, y_mam, s=pointsize*0.25, marker='o', color='green', zorder=1)
ax.plot(x_mam, y_mam, lw=2, color='green', zorder=1)

ax.scatter(SpT_man15, Mass_man15, s=pointsize,
           facecolors='', edgecolors='blue', label='Man15')  # Uncomment for Mass only

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)

# Axes: ticks
# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

ax.set_xticks(np.arange(-2, len(SpT_name)-1, 2))
ax.set_xticklabels(SpT_half)  #
ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)

# Axes: range & scale

ax.set_xlim(-2.5, 19)
ax.set_ylim(6e-2, 1.3e0)
ax.set_yscale('log')

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
