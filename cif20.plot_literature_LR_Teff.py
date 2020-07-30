# %%

# PLOT: LITERATURE > RADIUS - EFFECTIVE TEMPERATURE
# PLOT: LITERATURE > RADIUS - BOLOMETRIC LUMINOSITY
# Cifuentes et al. 2020
# Rabus et al. 2019

from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.ticker import FormatStrFormatter
import csv
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import sys
sys.path


rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_literature_Radius_Teff'
# filename = 'cif20_literature_Radius_Lbol'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

y_axis = 'R_Rsol'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    d_pc = []
    d_pc_ = []
    Teff = []
    Teff_ = []
    Lbol = []
    Lbol_ = []
    eLbol = []
    eTeff = []
    for row in reader:
        # and row[y_axis] != '':
        if row['d_pc'] != '' and row['Lbol'] != '' and row['Teff'] != '':
            SpT_.append(float(row['SpTnum']))
            d_pc_.append(float(row['d_pc']))
            Lbol_.append(float(row['Lbol']))
            Teff_.append(float(row['Teff']))
            if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                SpT.append(float(row['SpTnum']))
                d_pc.append(float(row['d_pc']))
                Lbol.append(float(row['Lbol']))
                eLbol.append(float(row['Lberr']))
                Teff.append(float(row['Teff']))


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


R_ = [Radius_SB(Lbol_[i], Teff_[i]) for i in range(len(Lbol_))]
R = [Radius_SB(Lbol[i], Teff[i]) for i in range(len(Lbol))]

# Rabus et al. 2019
# Page 5: "Interestingly, in Figure 4 we identified a discontinuous behaviour between 3200 and 3340 K (gray shaded area)"
# Equation 7: Ranges should be in mass (~0.23 Msol), but we use Teff because mass determination would come from Sch19.

Teff_rab19 = [3200, 3340]  #
Tsol = 5777  # K


def R_rab19_1(x):
    C = [-0.277, 0.869]
    eC = [0.060, 0.113]
    return C[0] + C[1]*x/Tsol


def R_rab19_2(x):
    C = [-1.223, 2.700]
    eC = [0.085, 0.0462]
    return C[0] + C[1]*x/Tsol


xp_rab19_1 = np.linspace(min(Teff), Teff_rab19[1], 1000)
xp_rab19_2 = np.linspace(Teff_rab19[0], max(Teff), 1000)
yp_rab19_1 = []
yp_rab19_2 = []
for i in range(len(xp_rab19_1)):
    yp_rab19_1.append(R_rab19_1(xp_rab19_1[i]))
    yp_rab19_2.append(R_rab19_2(xp_rab19_2[i]))

# Baraffe's Isochrones BCAH98
# http://perso.ens-lyon.fr/isabelle.baraffe/
# BCAH98_iso1 with log(t) = 9.6, 9.0 and 9.9 (1, 4.6, an 10 Ga) (remember, t_Sun = 4.6 Ga = 10^9.6 Ga)

isoc_bcah98 = ['BCAH98_iso1_logt_90.dat',
               'BCAH98_iso1_logt_99.dat', 'BCAH98_iso1_logt_96.dat']
isoc_dusty = ['DUSTY00_t_1.dat', 'DUSTY00_t_5.dat', 'DUSTY00_t_10.dat']


def Radius_logg(M, logg):
    # Units from IAU B.3 resolution: https://www.iau.org/static/resolutions/IAU2015_English.pdf
    G = 6.67408e-11  # m3 kg-1 s-2
    Msun = 1.9891*1e30  # kg
    Rsun = 6.957*1e8  # m
    R = np.sqrt(G*M*Msun/(1e-2*10**logg))  # cgs to mks
    return R/Rsun


Mass_bcah98 = [[] for i in range(len(isoc_bcah98))]
Teff_bcah98 = [[] for i in range(len(isoc_bcah98))]
logg_bcah98 = [[] for i in range(len(isoc_bcah98))]
Lbol_bcah98 = [[] for i in range(len(isoc_bcah98))]
Radius_bcah98 = [[] for i in range(len(isoc_bcah98))]

for i in range(len(isoc_bcah98)):
    with open('Baraffe/'+isoc_bcah98[i], 'r') as mycsv:
        #  m, Teff, g, logL, Mv,Mr, Mi, Mj, Mh, Mk, Ml', Mm
        dat_reader = csv.reader(mycsv)
        dat_ref = list(dat_reader)
        for j in range(0, len(dat_ref)-1):
            Mass_bcah98[i].append(float(dat_ref[j][0].split()[0]))  # Msol
            Teff_bcah98[i].append(float(dat_ref[j][0].split()[1]))  # K
            logg_bcah98[i].append(float(dat_ref[j][0].split()[2]))  # cm s-2
            Lbol_bcah98[i].append(10**float(dat_ref[j][0].split()[3]))  # Lsun

for i in range(len(isoc_bcah98)):
    for j in range(len(Mass_bcah98[i])):
        Radius_bcah98[i].append(Radius_logg(
            Mass_bcah98[i][j], logg_bcah98[i][j]))

# Baraffe's Isochrones DUSTY00
# http://perso.ens-lyon.fr/isabelle.baraffe/
# DUSTY00_models with t = 1, 5 and 10 Ga

Mass_dusty = [[] for i in range(len(isoc_dusty))]
Teff_dusty = [[] for i in range(len(isoc_dusty))]
logg_dusty = [[] for i in range(len(isoc_dusty))]
Lbol_dusty = [[] for i in range(len(isoc_dusty))]
Radius_dusty = [[] for i in range(len(isoc_dusty))]

for i in range(len(isoc_dusty)):
    with open('Baraffe/'+isoc_dusty[i], 'r') as mycsv:
        #  m, Teff, L, g, R, Mv, Li/L0, Mv, Mr, Mi, Mj, Mh, Mk, Ml', Mm
        dat_reader = csv.reader(mycsv)
        dat_ref = list(dat_reader)
        for j in range(0, len(dat_ref)-1):
            Mass_dusty[i].append(float(dat_ref[j][0].split()[0]))  # Msol
            Teff_dusty[i].append(float(dat_ref[j][0].split()[1]))  # K
            Lbol_dusty[i].append(10**float(dat_ref[j][0].split()[2]))  # Lsun
            logg_dusty[i].append(float(dat_ref[j][0].split()[3]))  # cm s-2
            Radius_dusty[i].append(float(dat_ref[j][0].split()[4]))  # Msol

# Variables
# Uncomment to select R or Lbol

# Teff

x_cif_ = Teff_
y_cif_ = R_
x_cif = Teff
y_cif = R
z_cif = SpT
x_bcah98 = Teff_bcah98
y_bcah98 = Radius_bcah98
x_dusty = Teff_dusty
y_dusty = Radius_dusty

# Lbol

# x_cif_ = Lbol_
# y_cif_ = R_
# x_cif = Lbol
# y_cif = R
# z_cif = SpT
# x_bcah98 = Lbol_bcah98
# y_bcah98 = Radius_bcah98
# x_dusty = Teff_dusty
# x_dusty = Lbol_dusty

# Shaded region
# Use for Teff only.

x_shadow = np.linspace(Teff_rab19[0], Teff_rab19[1], 100)

# PLOTTING

# Labels

xlabel = r'T$_{\rm eff}$ [K]'
# xlabel = r'$L$ [$L_{\rm sol}$]'

ylabel = r'$\mathcal{R}$ [$\mathcal{R}_\odot$]'
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

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma_r')

# Plot: distribution

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap, zorder=1)

# Plot: Rabus et al. 2019
# Plot: Isochrones by Baraffe
# R vs. Teff only

ax.plot(xp_rab19_1, yp_rab19_1, 'b--', lw=2, zorder=2)
ax.plot(xp_rab19_2, yp_rab19_2, 'b--', lw=2, zorder=2)
plt.fill_between(x_shadow, 0, 1, color='gainsboro', zorder=0)

# BCAH98

colour_bcah98 = ['grey', 'grey', 'black']

for i in range(len(isoc_bcah98)):
    ax.scatter(x_bcah98[i], y_bcah98[i], s=pointsize*0.25,
               facecolors='', edgecolors=colour_bcah98[i], zorder=1)
    ax.plot(x_bcah98[i], y_bcah98[i], c=colour_bcah98[i], lw=2, zorder=i)

# DUSTY00

colour_dusty = ['grey', 'grey', 'black']

for i in range(len(isoc_dusty)):
    ax.scatter(x_dusty[i], y_dusty[i], s=pointsize*0.25,
               facecolors='', edgecolors=colour_dusty[i], zorder=1)
    ax.plot(x_dusty[i], y_dusty[i], c=colour_bcah98[i], lw=2, zorder=i)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

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

plt.ylim(0.07, 0.75)  # Teff or Lbol
plt.xlim(2150, max(x_cif)+0)  # Teff (Avoid R of the coolest objects)
# plt.xlim(2e-4, 2e-1)  # Lbol
# ax.set_xscale('log')  # Lbol

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax)  # Colorbar
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)  # Keep uncommented to hide colorbar.
# cbar.set_ticks(np.arange(-2, 19, 2))
# cbar.ax.set_yticklabels(SpT_half)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
