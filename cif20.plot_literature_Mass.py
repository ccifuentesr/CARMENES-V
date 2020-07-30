# %%

# PLOT: LITERATURE > MASS - MASS
# Cifuentes et al. 2020
# Delfosse et al. 2000
# Benedict et al. 2016
# Mann et al. 2019

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
# Uncomment to choose one author to compare.

Mother_version = '01'
# filename = 'cif20_literature_mass_del00' #
# filename = 'cif20_literature_mass_ben16' #
filename = 'cif20_literature_mass_man19'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

mags = ['Ks']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    d_pc = []
    ed_pc = []
    d_pc_ = []
    ed_pc_ = []
    Lbol = []
    Lbol_ = []
    Teff = []
    Teff_ = []
    Ks_mag = []
    Ks_mag_ = []
    for row in reader:
        if row['Qf_'+mags[0]] == 'false' and row[mags[0]+'_mag'] != '':
            if row['d_pc'] != '' and row['Teff'] != '':
                if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                    d_pc_.append(float(row['d_pc']))
                    ed_pc_.append(float(row['ed_pc']))
                    Lbol_.append(float(row['Lbol']))
                    Teff_.append(float(row['Teff']))
                    Ks_mag_.append(float(row[mags[0]+'_mag']))

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


# Delfosse et al. 2000

def M_del00(Mabs):
    """Stellar mass from the empirical relation MK - Mass by Delfosse et al. 2000 (2000A&A...364..217D).
    Absolute magnitude in 2MASS Ks needs to be converted to CIT K magnitude before the transformation.

    Args:
        Mabs (float): Absolute magnitude in K(CIT).

    Returns:
        float: Stellar mass in solar units.

    (See page 220 in Astron. Astrophys. 364, 217–224, 2000).
    """
    C = [1.8, 6.12, 13.205, -6.2315, 0.37529]
    M = 10**((C[0] + C[1]*Mabs + C[2]*Mabs**2 +
              C[3]*Mabs**3 + C[4]*Mabs**4)*1E-3)
    return M


def JHK_CIT(J_2MASS, H_2MASS, K_2MASS):
    """Transformation from 2MASS JHKs to CIT JHK systems from Carpenter et al. 2000 (2001AJ....121.2851C).
    The Caltech (CIT) photometric system was described by Frogel et al. (1978) and defined by the 
    standard-star observations from Elias et al. (1982).
    Valid in the range MK ∈ [4.5, 9.5] mag.

    Args:
        J_2MASS (float): magnitude in 2MASS J passband.
        H_2MASS (float): magnitude in 2MASS H passband.
        K_2MASS (float): magnitude in 2MASS Ks passband.

    Returns:
        float: magnitude in CIT J passband. 
        float: magnitude in CIT H passband. 
        float: magnitude in CIT K passband. 

    (See §4.3 in Carpenter et al. 2000).
    """
    K_CIT = K_2MASS + 0.024  # Masses by Del00 only make use of this transformation
    J_CIT = 1/1.056 * (J_2MASS - K_2MASS + 0.013) + K_CIT
    H_CIT = 1/1.026 * (H_2MASS - K_2MASS + 0.028) + K_CIT
    return J_CIT, H_CIT, K_CIT

# Literature: Benedict et al. 2016


def M_ben16(Mabs):
    """Stellar mass from the empirical relation MKs - Mass by Benedict et al. 2016 (2016AJ....152..141B).
    Valid in the range MK <= 10 mag.

    Args:
        Mabs (float): Absolute magnitude in Ks.

    Returns:
        float: Stellar mass in solar units.

    (See Equation 11 and Table 13 in Benedict et al. 2016).
    """
    C = [0.2311, -0.1352, 0.0400, 0.0038, -0.0032]
    Cerr = [0.0004, 0.0007, 0.0005, 0.0002, 0.0001]
    x0 = 7.5
    M = C[0] + C[1]*(Mabs - x0) + C[2]*(Mabs - x0)**2 + \
        C[3]*(Mabs - x0)**3 + C[4]*(Mabs - x0)**4
    return M


# Mann et al. 2019

def M_man19(Mabs):
    """Stellar mass from the empirical relation MKs - Mass by Mann et al. 2019 (2019ApJ...871...63M).
    Valid in the range MKs ∈ [0.4, 11.0] mag ('safe' [0.45, 10.5] mag).

    Args:
        Mabs (float): Absolute magnitude in Ks.

    Returns:
        float: Stellar mass in solar units.

    (See Table 6 in Mann et al. 2019; n = 5 fit is preferred).
    """
    C = [-0.642, -0.208, -8.43*1e-4, 7.87*1e-3, 1.42*1e-4, -2.13*1e-4]
    zp = 7.5
    M = 10**(C[0] + C[1]*(Mabs - zp) + C[2]*(Mabs - zp)**2 + C[3] *
             (Mabs - zp)**3 + C[4]*(Mabs - zp)**4 + C[5]*(Mabs - zp)**5)
    return M


# Variables
# Define the validity limits for each case.
# Uncomment only one case.

# This work
# Keep for ALL
Teff_ = [Teff_[i] + 95 for i in range(len(Teff_))]

Radius_MR_ = [Radius_SB(Lbol_[i], Teff_[i]) for i in range(len(Lbol_))]
Mass_MR_ = [Mass_sch19(Radius_MR_[i]) for i in range(len(Radius_MR_))]

xp = np.linspace(0, 1, 100)
yp = xp

# Delfosse et al. 2000

# # ALL
# Mabs_del00_ = [(Ks_mag_[i] + 0.024) - 5 *
#                np.log10(d_pc_[i]) + 5 for i in range(len(Ks_mag_))]  # K_CIT = K_2MASS + 0.024
# Mass_del00_ = [M_del00(Mabs_del00_[i]) for i in range(len(Mabs_del00_))]

# # VALID in MK ∈ [4.5, 9.5]
# Mabs_del00_indices = [i for i in range(len(Mabs_del00_)) if Mabs_del00_[
#     i] > 4.5 and Mabs_del00_[i] < 9.5]
# Mabs_del00 = [Mabs_del00_[Mabs_del00_indices[i]]
#               for i in range(len(Mabs_del00_indices))]
# Mass_del00 = [M_del00(Mabs_del00[i]) for i in range(len(Mabs_del00))]

# # Masses: This work
# Lbol_del00 = [Lbol_[Mabs_del00_indices[i]]
#               for i in range(len(Mabs_del00_indices))]
# Teff_del00 = [Teff_[Mabs_del00_indices[i]]
#               for i in range(len(Mabs_del00_indices))]
# Radius_MR_del00 = [Radius_SB(Lbol_del00[i], Teff_del00[i])
#                    for i in range(len(Lbol_del00))]
# Mass_MR_del00 = [Mass_sch19(Radius_MR_del00[i])
#                  for i in range(len(Radius_MR_del00))]

# # Plotting
# x_ = Mass_MR_
# y_ = Mass_del00_
# z_ = Mabs_del00_
# x = Mass_MR_del00
# y = Mass_del00
# z = Mabs_del00

# Benedict et al. 2016

# # ALL
# Mabs_ben16_ = [(Ks_mag_[i]) - 5 *
#                np.log10(d_pc_[i]) + 5 for i in range(len(Ks_mag_))]
# Mass_ben16_ = [M_ben16(Mabs_ben16_[i]) for i in range(len(Mabs_ben16_))]

# # VALID (Mabs <= 10)
# Mabs_ben16_indices = [i for i in range(
#     len(Mabs_ben16_)) if Mabs_ben16_[i] <= 10]
# Mabs_ben16 = [Mabs_ben16_[Mabs_ben16_indices[i]]
#               for i in range(len(Mabs_ben16_indices))]
# Mass_ben16 = [M_ben16(Mabs_ben16[i]) for i in range(len(Mabs_ben16))]

# # Masses: This work
# Lbol_ben16 = [Lbol_[Mabs_ben16_indices[i]]
#               for i in range(len(Mabs_ben16_indices))]
# Teff_ben16 = [Teff_[Mabs_ben16_indices[i]]
#               for i in range(len(Mabs_ben16_indices))]
# Radius_MR_ben16 = [Radius_SB(Lbol_ben16[i], Teff_ben16[i])
#                    for i in range(len(Lbol_ben16))]
# Mass_MR_ben16 = [Mass_sch19(Radius_MR_ben16[i])
#                  for i in range(len(Radius_MR_ben16))]

# # Plotting
# x_ = Mass_MR_
# y_ = Mass_ben16_
# z_ = Mabs_ben16_
# x = Mass_MR_ben16
# y = Mass_ben16
# z = Mabs_ben16

# Mann et al. 2019

# ALL
Mabs_man19_ = [Ks_mag_[i] - 5 *
               np.log10(d_pc_[i]) + 5 for i in range(len(Ks_mag_))]
Mass_man19_ = [M_man19(Mabs_man19_[i]) for i in range(len(Mabs_man19_))]

# VALID (MKs ∈ [0.4, 11.0])
Mabs_man19_indices = [i for i in range(len(Mabs_man19_)) if Mabs_man19_[
    i] > 4.5 and Mabs_man19_[i] < 10.5]
Mabs_man19 = [Mabs_man19_[Mabs_man19_indices[i]]
              for i in range(len(Mabs_man19_indices))]
Mass_man19 = [M_man19(Mabs_man19[i]) for i in range(len(Mabs_man19))]

# Masses: This work
Lbol_man19 = [Lbol_[Mabs_man19_indices[i]]
              for i in range(len(Mabs_man19_indices))]
Teff_man19 = [Teff_[Mabs_man19_indices[i]]
              for i in range(len(Mabs_man19_indices))]
Radius_MR_man19 = [Radius_SB(Lbol_man19[i], Teff_man19[i])
                   for i in range(len(Lbol_man19))]
Mass_MR_man19 = [Mass_sch19(Radius_MR_man19[i])
                 for i in range(len(Radius_MR_man19))]

# Plotting
x_ = Mass_MR_
y_ = Mass_man19_
z_ = Mabs_man19_
x = Mass_MR_man19
y = Mass_man19
z = Mabs_man19

# PLOTTING

# Labels
# Uncomment to select the appropriate authors.

# ylabel = r'$\mathcal{M}$ [$\mathcal{M}_\odot$] $_{\rm Del00}$' #
# ylabel = r'$\mathcal{M}$ [$\mathcal{M}_\odot$] $_{\rm Ben16}$'
ylabel = r'$\mathcal{M}$ [$\mathcal{M}_\odot$] $_{\rm Man19}$'

xlabel = r'$\mathcal{M}$ [$\mathcal{M}_\odot$] $_{\rm This~work}$'
cbarlabel = r'M$_{Ks}$ [mag]'

# Sizes

figsize = (12, 10)
pointsize = 40
tickssize = 22
labelsize = 22
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('viridis_r')

# Plots: distribution

sc_ = plt.scatter(x_, y_, s=pointsize, edgecolors='gainsboro',
                  facecolors='', zorder=0)  # All
sc = plt.scatter(x, y, s=pointsize, c=z, cmap=cmap, zorder=1)  # Valid

# Plots: equality

plt.plot(xp, yp, 'k--', lw=2, zorder=2)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
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

ax.set_xlim(0.06, 0.75)
ax.set_ylim(0.06, 0.75)

# Colorbar

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
cbar.ax.tick_params(labelsize=labelsize)
cbar.ax.invert_yaxis()
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()

# Stats

res = [(y[i]-x[i]) for i in range(len(x))]
res_per = [(y[i]-x[i])/x[i]*100 for i in range(len(x))]

res_mean = np.mean(res)
res_std = np.std(res)
res_per_mean = np.mean(res_per)
res_per_std = np.std(res_per)

print('Mean:', np.round(res_mean, 4), 'Stdev:', np.round(res_std, 4))
print('Mean (%):', np.round(res_per_mean, 4),
      'Stdev (%):', np.round(res_per_std, 4))
