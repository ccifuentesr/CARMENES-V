# %%

# PLOT: LITERATURE > BOLOMETRIC LUMINOSITIES
# PLOT: LITERATURE > EFFECTIVE TEMPERATURES
# PLOT: LITERATURE > MASSES
# PLOT: LITERATURE > RADII
# Cifuentes et al. 2020

# Comparison with literature values from many authors
# (see complete list in Cifuentes et al. 2020)

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams.update({'errorbar.capsize': 4})

# Data
# Booleans allow to select only clean data.
# Uncomment for L / T / M / R diagram.

Mother_version = '01'
booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# filename = 'cif20_literature_luminosities'
# filename = 'cif20_literature_Teff'
filename = 'cif20_literature_masses'
# filename = 'cif20_literature_radii'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    Lbol = []
    Lbol_lit = []
    Lberr = []
    Lberr_lit = []
    Teff = []
    Teff_lit = []
    eTeff_lit = []
    for row in reader:
        if row['Lbol'] != '' and row['Lbol_lit'] != '' and row['Lberr_lit'] != '':  # Lbol
            if row['Teff_lit'] != '':
                SpT.append(float(row['SpTnum']))
                Lbol_lit.append(float(row['Lbol_lit']))
                Lberr_lit.append(float(row['Lberr_lit']))
                Lbol.append(float(row['Lbol']))
                Lberr.append(float(row['Lberr']))
                Teff.append(float(row['Teff']))
                Teff_lit.append(float(row['Teff_lit']))
                eTeff_lit.append(float(row['eTeff_lit']))

# Exclusive append for M and R plots

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Karmn = []
    SpT_ = []
    Lbol_ = []
    Lbol_lit_ = []
    Lberr_ = []
    Lberr_lit_ = []
    Teff_ = []
    Teff_lit_ = []
    eTeff_lit_ = []
    R_lit = []
    eR_lit = []
    M_lit = []
    eM_lit = []
    for row in reader:
        if row['Lbol'] != '' and row['Lbol_lit'] != '' and row['Lberr_lit'] != '':  # Lbol
            if row['R_lit'] != '' and row['eR_lit'] != '':  # Radius
                if row['M_lit'] != '' and row['eM_lit'] != '':  # Mass
                    if row['Teff'] != '' and row['Teff_lit'] != '' and row['eTeff_lit'] != '':  # Teff
                        SpT_.append(float(row['SpTnum']))
                        Lbol_lit_.append(float(row['Lbol_lit']))
                        Lberr_lit_.append(float(row['Lberr_lit']))
                        Lbol_.append(float(row['Lbol']))
                        Lberr_.append(float(row['Lberr']))
                        Teff_.append(float(row['Teff']))
                        Teff_lit_.append(float(row['Teff_lit']))
                        eTeff_lit_.append(float(row['eTeff_lit']))
                        R_lit.append(float(row['R_lit']))
                        eR_lit.append(float(row['eR_lit']))
                        M_lit.append(float(row['M_lit']))
                        eM_lit.append(float(row['eM_lit']))

# Gaia DR2 values

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_gaia_cif = []
    Lbol_gaia = []
    Lbol_gaia_cif = []
    eLbol_gaia_cif = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
            if row['Teff'] != '':
                if row['lum_val_1'] != '':
                    SpT_gaia_cif.append(float(row['SpTnum']))
                    Lbol_gaia_cif.append(float(row['Lbol']))
                    eLbol_gaia_cif.append(float(row['Lberr']))
                    Lbol_gaia.append(float(row['lum_val_1']))

# Passegger et al. 19

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_pas19 = []
    Teff_pas19 = []
    TeffVN_pas19 = []
    eTeffVN_pas19 = []
    for row in reader:
        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
            if row['Teff'] != '':
                if row['Teff_Pas19'] != '':
                    SpT_pas19.append(float(row['SpTnum']))
                    Teff_pas19.append(float(row['Teff']))
                    TeffVN_pas19.append(float(row['Teff_Pas19']))
                    eTeffVN_pas19.append(float(row['eTeff_Pas19']))

# Literature: Schweitzer+19
# http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/625/A68
# VOSA values end with _cif

with open('Literature/schweitzer19.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_sch19 = []
    Teff_sch19 = []
    eTeff_sch19 = []
    R_sch19 = []
    eR_sch19 = []
    M_sch19 = []
    eM_sch19 = []
    Lbol_sch19_cif = []
    Lberr_sch19_cif = []
    Teff_sch19_cif = []
    for row in reader:
        if row['Lbol'] != '':
            if row['Com'] != 'Out':
                SpT_sch19.append(float(row['SpTnum']))
                Teff_sch19.append(float(row['Teff']))
                eTeff_sch19.append(float(row['e_Teff']))
                R_sch19.append(float(row['Rad']))
                eR_sch19.append(float(row['e_Rad']))
                M_sch19.append(float(row['MassMR']))
                eM_sch19.append(float(row['e_MassMR']))
                Lbol_sch19_cif.append(float(row['Lbol']))
                Lberr_sch19_cif.append(float(row['Lberr']))
                Teff_sch19_cif.append(float(row['Teff_cif']))


Teff_sch19_cif = [Teff_sch19_cif[i] +
                  delta_T for i in range(len(Teff_sch19_cif))]

# Radius: Stefan-Boltzmann (SB)


def Radius_SB(Lbol, Lberr, Teff, eTeff):
    """Stellar radius and its error from the Stefan–Boltzmann law under the black body approximation.

    Args:
        Lbol (float): Bolometric luminosity in solar units.
        Lberr (float): Bolometric luminosity uncertainty in solar units.
        Teff (float): Effective temperature in Kelvin.
        eTeff (float): Effective temperature uncertainty in Kelvin.

    Returns:
        float: Stellar radius in solar units.
        float: Stellar radius error in solar units.

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
    sigma = 5.670367*1e-8
    a = (Lbol*Lsun)/(4*np.pi*sigma*Teff**4*Rsun**2)
    R = 1/Rsun * np.sqrt(Lbol*Lsun/(4*np.pi*sigma*Teff**4))
    eR = np.sqrt(a * ((Lberr/(2*Lbol))**2 + (-2*eTeff/Teff)**2))
    return R, eR

# Masses: R-M relation for eclipsing binaries (Schweitzer et al. 2019)


def Mass_sch19(Radius, eRadius):
    """Stellar mass and its error from the empirical relation by Schweitzer et al. 2019 
    (2019A&A...625A..68S), based on masses and radii of eclipsing binaries.

    Args:
        Radius (float): Stellar radius in solar units.
        eRadius (float): Stellar radius uncertainty in solar units.

    Returns:
        float: Stellar mass in solar units.
        float: Stellar mass error in solar units.

    (See Equation 6 in Schweitzer et al. 2019 and references therein).
    """
    a = -0.024048024
    b = 1.0552427
    a_err = 0.007592668
    b_err = 0.017044148
    M = a + b * Radius
    eM = np.sqrt((a_err)**2 + (Radius * b_err)**2 + (b * eRadius)**2)
    return M, eM


# Variables

# Preliminar computations

eTeff = [50 for i in range(len(Teff))]
eTeff_ = [50 for i in range(len(Teff_))]
eTeff_pas19 = [50 for i in range(len(Teff_pas19))]
eTeff_sch19_cif = [50 for i in range(len(Teff_sch19_cif))]

R = [Radius_SB(Lbol_[i], Lberr_[i], Teff_[i], eTeff_[i])[0]
     for i in range(len(Lbol_))]
eR = [Radius_SB(Lbol_[i], Lberr_[i], Teff_[i], eTeff_[i])[1]
      for i in range(len(Lbol_))]
R_sch19_cif = [Radius_SB(Lbol_sch19_cif[i], Lberr_sch19_cif[i], Teff_sch19_cif[i], eTeff_sch19_cif[i])[0]
               for i in range(len(Lbol_sch19_cif))]
eR_sch19_cif = [Radius_SB(Lbol_sch19_cif[i], Lberr_sch19_cif[i], Teff_sch19_cif[i], eTeff_sch19_cif[i])[1]
                for i in range(len(Lbol_sch19_cif))]

M = [Mass_sch19(R[i], eR[i])[0]
     for i in range(len(R))]
eM = [Mass_sch19(R[i], eR[i])[1]
      for i in range(len(R))]
M_sch19_cif = [Mass_sch19(R_sch19_cif[i], eR_sch19_cif[i])[0]
               for i in range(len(R_sch19_cif))]
eM_sch19_cif = [Mass_sch19(R_sch19_cif[i], eR_sch19_cif[i])[1]
                for i in range(len(R_sch19_cif))]

# Uncomment for L / T / M / R diagram.

# Lbol

# x = Lbol
# y = Lbol_lit
# z = SpT
# xerr = Lberr
# yerr = Lberr_lit

# x_gaia = Lbol_gaia_cif
# y_gaia = Lbol_gaia
# z_gaia = SpT_gaia_cif
# xerr_gaia = eLbol_gaia_cif
# yerr_gaia = 0

# xp = np.linspace(0, 1, 100)
# yp = xp

# Uncomment for Teff

# x = Teff
# y = Teff_lit
# z = SpT
# xerr = eTeff
# yerr = eTeff_lit

# x_pas19 = Teff_pas19
# y_pas19 = TeffVN_pas19
# z_pas19 = SpT_pas19
# xerr_pas19 = eTeff_pas19
# yerr_pas19 = eTeffVN_pas19

# xp = np.linspace(1000, 5000, 1000)
# yp = xp

# Uncomment for Radius

# x = R
# y = R_lit
# z = SpT_
# xerr = eR
# yerr = eR_lit

# x_sch19 = R_sch19_cif
# y_sch19 = R_sch19
# z_sch19 = SpT_sch19
# xerr_sch19 = eR_sch19_cif
# yerr_sch19 = eR_sch19

# xp = np.linspace(0, 1, 100)
# yp = xp

# Uncomment for Mass

x = M
y = M_lit
z = SpT_
xerr = eM
yerr = eM_lit

x_sch19 = M_sch19_cif
y_sch19 = M_sch19
z_sch19 = SpT_sch19
xerr_sch19 = M_sch19_cif
yerr_sch19 = eM_sch19

xp = np.linspace(0, 1, 100)
yp = xp

# PLOTTING

# Labels
# Uncomment for L / T / M / R diagram.

# Lbol
# xlabel = r'$L$ [$L_{\rm sol}$] $_{\rm This~work}$'
# ylabel = r'$L$ [$L_{\rm sol}$] $_{\rm Literature}$'

# Teff
# xlabel = r'$T_{\rm eff}$ [K] $_{\rm (This\,work)}$'
# ylabel = r'$T_{\rm eff}$ [K] $_{\rm (Literature)}$'

# Radius
# xlabel = r'$\mathcal{R}$ [$\mathcal{R}_\odot$] $_{\rm This~work}$'
# ylabel = r'$\mathcal{R}$ [$\mathcal{R}_\odot$] $_{\rm Literature}$'

# Mass
xlabel = r'$\mathcal{M}$ [$\mathcal{M}_\odot$] $_{\rm This~work}$'
ylabel = r'$\mathcal{M}$ [$\mathcal{M}_\odot$] $_{\rm Literature}$'

cbarlabel = r'Spectral type'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma_r')

# Plot: distribution

sc = ax.scatter(x, y, c=z, cmap=cmap, s=pointsize,
                marker='o', edgecolor='none', zorder=1)
plt.plot(xp, yp, 'k--', lw=2)

# Uncomment to plot grey errorbars

ax.errorbar(x, y, xerr=xerr, yerr=yerr, lw=1,
            ls='none', capsize=3, color='grey', zorder=0)

# Plots: literature
# Uncomment the following for L / T / M / R diagram.

# Uncomment only for Lbol
# sc_ = ax.scatter(x_gaia, y_gaia, c=z_gaia, cmap=cmap, s=pointsize,
#                  marker='o', edgecolor='none', zorder=1)
# ax.errorbar(x_gaia, y_gaia, xerr=xerr_gaia, yerr=yerr_gaia, lw=1,
#             ls='none', capsize=3, color='grey', zorder=0)

# Uncomment only for Teff

# sc = ax.scatter(x_pas19, y_pas19, c=z_pas19, cmap=cmap, s=pointsize,
#                 marker='o', edgecolor='none', zorder=1)
# ax.errorbar(x_pas19, y_pas19,  xerr=xerr_pas19, yerr=yerr_pas19, lw=1,
#             ls='none', capsize=3, color='grey', zorder=0)

# Uncomment only for R or M
sc = ax.scatter(x_sch19, y_sch19, c=z_sch19, cmap=cmap, s=pointsize,
                marker='o', edgecolor='none', zorder=1)

# Aesthetics

# Axes: labels & legend
# Uncomment the following for L / T / M / R diagram.

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))  # L, M, R
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))  # L, M, R
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))  # T
# ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))  # T

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
# Uncomment the following for L / T / M / R diagram.

# Uncomment for L
# ax.set_xlim(4E-5, 2E-1)
# ax.set_ylim(4E-5, 2E-1)
# ax.set_xscale('log')
# ax.set_yscale('log')

# Uncomment for T
# dx = 150
# ax.set_xlim(min(x)-dx, max(x)+dx)
# ax.set_ylim(min(x)-dx, max(x)+dx)

# Uncomment for R or M
ax.set_xlim(0.06, 0.75)
ax.set_ylim(0.06, 0.75)

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)  # Colorbar
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)
sc_.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()


# %%
# To csv
filename = 'cif20.MR_literature.csv'
with open(filename, mode='w') as mycsv:
    writer = csv.writer(mycsv, delimiter='\n')
    header = []
    rows = []
    rows.append("Karmn, SpTnum, R, eR, M, eM, R_lit, eR_lit, M_lit, eM_lit")
    for i in range(len(Karmn)):
        rows.append("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(
            Karmn[i], SpT[i], np.round(R[i], 4), np.round(eR[i], 4), np.round(M[i], 4), np.round(eM[i], 4), R_lit[i], eR_lit[i], M_lit[i], eM_lit[i]))
    writer.writerow(rows)
