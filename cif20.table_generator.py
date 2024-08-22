# %%

# GENERATOR: PARAMETERS
# PLOT: MEDIAN VALUES
# Cifuentes et al. 2020 (Table 6)
# Average astrophysical parameters for K5 V to L2.0 objects.

# Read https://stackoverflow.com/questions/34050491/standard-deviation-in-numpy
# for 'ddof=1' in numpy std.
# Output ready for pasting into LaTeX tabular environment.

# Example of data access:
# mags[0] means K5 V, all filters, all stars;
# mags[0][0] means K5 V, FUV filter, all stars;
# mags[0][0][0] means K5 V, FUV filter, first star.
# Teff[0] means all Teffs for K5 V;
# Teff[0][0] means all Teffs for K5 V and first star.

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import pyperclip
import math
from astropy.stats import sigma_clip
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

Mother_version = '01'
booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Filter and spectral sequence lists

filters = ['FUV', 'NUV', 'u', 'BT', 'B', 'g', 'BP', 'VT', 'V',
           'GG', 'r', 'i', 'RP', 'J', 'H', 'Ks', 'W1', 'W2', 'W3', 'W4']
SpT_name = ['K5\,V', 'K7\,V', 'M0.0\,V', 'M0.5\,V', 'M1.0\,V', 'M1.5\,V', 'M2.0\,V', 'M2.5\,V', 'M3.0\,V', 'M3.5\,V', 'M4.0\,V', 'M4.5\,V', 'M5.0\,V', 'M5.5\,V', 'M6.0\,V', 'M6.5\,V', 'M7.0\,V',
            'M7.5\,V', 'M8.0\,V', 'M8.5\,V', 'M9.0\,V', 'M9.5\,V', 'L0', 'L0.5', 'L1', 'L1.5', 'L2', 'L2.5', 'L3', 'L3.5', 'L4', 'L4.5', 'L5', 'L5.5', 'L6', 'L6.5', 'L7', 'L7.5', 'L8']
SpT_range_pre = np.arange(0, 18.5, 0.5)  # M0.0 to L8.0
SpT_range = [-2, -1]  # K5 and K7
for i in range(len(SpT_range_pre)):
    SpT_range.append(np.ndarray.tolist(
        np.arange(0, 18.5, 0.5))[i])  # all range

alt_size = 0.8  # To resize all points proportionally

# Variables
# One column, number of filters long, for each spectral type.

mags = [[[] for i in range(len(filters))]
        for j in range(len(SpT_range))]
Mbol = [[[] for i in range(len(filters))]
        for j in range(len(SpT_range))]
BC_all = [[[] for i in range(len(filters))]
          for j in range(len(SpT_range))]
Lbol_all = [[[] for i in range(len(filters))]
            for j in range(len(SpT_range))]
Teff_all = [[[] for i in range(len(filters))]
            for j in range(len(SpT_range))]
Radius_all = [[[] for i in range(len(filters))]
              for j in range(len(SpT_range))]
Mass_all = [[[] for i in range(len(filters))]
            for j in range(len(SpT_range))]

mags_mean = [[] for i in range(len(SpT_range))]
mags_median = [[] for i in range(len(SpT_range))]
mags_std = [[] for i in range(len(SpT_range))]
mags_size = [[] for i in range(len(SpT_range))]

# Parameters

SpTnum = [[[] for i in range(len(filters))] for j in range(len(SpT_range))]
d_pc = [[[] for i in range(len(filters))] for j in range(len(SpT_range))]
Lbol = [[] for i in range(len(SpT_range))]
Teff = [[] for i in range(len(SpT_range))]
logg = [[] for i in range(len(SpT_range))]

Radius_all = [[] for i in range(len(SpT_range))]
Mass_all = [[] for i in range(len(SpT_range))]

# Mother: all values

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_cif = []
    SpT_cif_ = []
    Teff_cif = []
    Teff_cif_ = []
    Lbol_cif = []
    eLbol_cif = []
    Lbol_cif_ = []
    eLbol_cif_ = []
    for row in reader:
        if row['Teff'] != '' and row['Lbol'] != '':
            SpT_cif_.append(float(row['SpTnum']))
            Teff_cif_.append(float(row['Teff']))
            Lbol_cif_.append(float(row['Lbol']))
            eLbol_cif_.append(float(row['Lberr']))
            if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                SpT_cif.append(float(row['SpTnum']))
                Teff_cif.append(float(row['Teff']))
                Lbol_cif.append(float(row['Lbol']))
                eLbol_cif.append(float(row['Lberr']))

eTeff_cif_ = [50 for i in range(len(Teff_cif_))]
eTeff_cif = [50 for i in range(len(Teff_cif))]

# Mother: mean values

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    for row in Mother:
        for m in range(len(SpT_range)-1):
            for n in range(len(filters)):
                if float(row['SpTnum']) == SpT_range[m]:  # Choose Spectral type
                    # Only clean data
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        # Only valid data
                        if row[filters[n]+'_mag'] != '' and row['e'+filters[n]+'_mag'] != '':
                            if row['Qf_'+filters[n]] == 'false':  # Clean photometry
                                if row['Lbol'] != '':  # To compute BC
                                    SpTnum[m][n].append(float(row['SpTnum']))
                                    d_pc[m][n].append(float(row['d_pc']))
                                    mags[m][n].append(
                                        float(row[filters[n]+'_mag']))
                                    Lbol_all[m][n].append(float(row['Lbol']))
                                    Teff_all[m][n].append(float(row['Teff']))

# Mother: VOSA parameters

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    for row in Mother:
        for m in range(len(SpT_range)):  # only stars with VOSA
            if float(row['SpTnum']) == SpT_range[m]:  # Choose Spectral type
                if row['Teff'] != '':
                    # Only clean stars
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        Lbol[m].append(float(row['Lbol']))
                        Teff[m].append(float(row['Teff']))
                        logg[m].append(float(row['logg']))

# Parameters statistics
# Produce combined lists for rolling median
# K5, K7, K7 + M0.0 + M0.5, M0.0 + M0.5 + M1.0, etc.

Lbol_roll = []
for i in np.arange(2, len(SpT_range)-1):
    Lbol_roll.append(Lbol[i-1] + Lbol[i] + Lbol[i+1])

Lbol_roll.insert(0, Lbol[0])
Lbol_roll.insert(1, Lbol[1])
Lbol_roll.insert(len(Lbol_roll), Lbol[-1])

Teff_roll = []
for i in np.arange(2, len(SpT_range)-1):
    Teff_roll.append(Teff[i-1] + Teff[i] + Teff[i+1])

Teff_roll.insert(0, Teff[0])
Teff_roll.insert(1, Teff[1])
Teff_roll.insert(len(Teff_roll), Teff[-1])

# Lbol
# Sigma-clipped median values

sigma = 2.5

Lbol_median = []
Lbol_std = []
Lbol_size = []

for i in range(len(SpT_range)):
    Lbol_median.append(np.ma.median(sigma_clip(
        Lbol_roll[i], sigma=sigma, cenfunc='median', masked='False')))
    Lbol_std.append(np.std(sigma_clip(
        Lbol_roll[i], sigma=sigma, cenfunc='median', masked='False')))
    Lbol_size.append(len(sigma_clip(
        Lbol_roll[i], sigma=sigma, cenfunc='median', masked='False')))

Lbol_cif_clipped = []
for i in range(len(SpT_range)):
    Lbol_cif_clipped.append(sigma_clip(
        Lbol[i], sigma=sigma, cenfunc='median', masked='False'))

SpT_cif_clipped = []
for i in range(len(SpT_range)):
    SpT_cif_clipped.append(len(Lbol_cif_clipped[i])*[SpT_range[i]])

# Lbol for BC calculations

Lbol_all_roll = [[] for i in range(len(SpT_range))]
for m in np.arange(2, len(SpT_range)-1):
    for n in range(len(filters)):
        Lbol_all_roll[m].append(
            Lbol_all[m-1][n] + Lbol_all[m][n] + Lbol_all[m+1][n])

for n in reversed(range(len(filters))):
    Lbol_all_roll[0].insert(0, Lbol_all[0][n])
    Lbol_all_roll[1].insert(0, Lbol_all[1][n])
    Lbol_all_roll[-1].insert(0, Lbol_all[-1][n])

# Teff

Teff_median = []
Teff_std = []
Teff_size = []

for i in range(len(SpT_range)):
    Teff_median.append(np.ma.median(sigma_clip(
        Teff_roll[i], sigma=sigma, cenfunc='median', masked='False')))
    Teff_std.append(np.std(sigma_clip(
        Teff_roll[i], sigma=sigma, cenfunc='median', masked='False')))
    Teff_size.append(len(sigma_clip(
        Teff_roll[i], sigma=sigma, cenfunc='median', masked='False')))

# Photometry

for m in range(len(SpT_range)):
    for n in range(len(filters)):
        mags_median[m].append(np.median(mags[m][n]))
        mags_std[m].append(np.std(mags[m][n]))
        mags_size[m].append(len(mags[m][n]))

# Absolute magnitudes

Mabs = [[] for i in range(len(SpT_range))]
Mabs_median = [[] for i in range(len(SpT_range))]
Mabs_std = [[] for i in range(len(SpT_range))]

for m in range(len(SpT_range)):
    for n in range(len(filters)):
        Mabs[m].append(((mags[m][n]) - 5 * np.log10(d_pc[m][n]) + 5).tolist())

Mabs_roll = [[] for i in range(len(SpT_range))]
for m in np.arange(2, len(SpT_range)-1):
    for n in range(len(filters)):
        Mabs_roll[m].append(Mabs[m-1][n] + Mabs[m][n] + Mabs[m+1][n])

for n in reversed(range(len(filters))):
    Mabs_roll[0].insert(0, Mabs[0][n])
    Mabs_roll[1].insert(0, Mabs[1][n])
    Mabs_roll[-1].insert(0, Mabs[-1][n])

for m in np.arange(len(SpT_range)):
    for n in range(len(filters)):
        Mabs_median[m].append(np.ma.median(sigma_clip(
            Mabs_roll[m][n], sigma=sigma, cenfunc='median', masked='False')))
        Mabs_std[m].append(np.std(sigma_clip(
            Mabs_roll[m][n], sigma=sigma, cenfunc='median', masked='False')))

# Absolute magnitudes: individual filters

Mabs_B = []
eMabs_B = []
Mabs_g = []
eMabs_g = []
Mabs_BP = []
eMabs_BP = []
Mabs_V = []
eMabs_V = []
Mabs_r = []
eMabs_r = []
Mabs_GG = []
eMabs_GG = []
Mabs_i = []
eMabs_i = []
Mabs_RP = []
eMabs_RP = []
Mabs_J = []
eMabs_J = []
Mabs_H = []
eMabs_H = []
Mabs_Ks = []
eMabs_Ks = []
Mabs_W1 = []
eMabs_W1 = []
Mabs_W2 = []
eMabs_W2 = []
Mabs_W3 = []
eMabs_W3 = []

for i in range(len(Mabs_median)):
    Mabs_B.append(np.median(Mabs_median[i][4]))
    Mabs_g.append(np.median(Mabs_median[i][5]))
    Mabs_BP.append(np.median(Mabs_median[i][6]))
    Mabs_V.append(np.median(Mabs_median[i][8]))
    Mabs_GG.append(np.median(Mabs_median[i][9]))
    Mabs_r.append(np.median(Mabs_median[i][10]))
    Mabs_i.append(np.median(Mabs_median[i][11]))
    Mabs_RP.append(np.median(Mabs_median[i][12]))
    Mabs_J.append(np.median(Mabs_median[i][13]))
    Mabs_H.append(np.median(Mabs_median[i][14]))
    Mabs_Ks.append(np.median(Mabs_median[i][15]))
    Mabs_W1.append(np.median(Mabs_median[i][16]))
    Mabs_W2.append(np.median(Mabs_median[i][17]))
    Mabs_W3.append(np.median(Mabs_median[i][18]))

for i in range(len(Mabs_std)-1):
    eMabs_B.append(Mabs_std[i][4])
    eMabs_g.append(Mabs_std[i][5])
    eMabs_BP.append(Mabs_std[i][6])
    eMabs_V.append(Mabs_std[i][8])
    eMabs_GG.append(Mabs_std[i][9])
    eMabs_r.append(Mabs_std[i][10])
    eMabs_i.append(Mabs_std[i][11])
    eMabs_RP.append(Mabs_std[i][12])
    eMabs_J.append(Mabs_std[i][13])
    eMabs_H.append(Mabs_std[i][14])
    eMabs_Ks.append(Mabs_std[i][15])
    eMabs_W1.append(Mabs_std[i][16])
    eMabs_W2.append(Mabs_std[i][17])
    eMabs_W3.append(Mabs_std[i][18])

# Bolometric corrections
# BC = Mbol - Mabs

for m in range(len(SpT_range)):
    for n in range(len(filters)):
        for i in range(len(Lbol_all_roll[m][n])):
            Lsun = 3.8275*1e26  # W
            BC_all[m][n].append(
                (71.197425-2.5*np.log10(Lbol_all_roll[m][n][i]*Lsun)) - Mabs_roll[m][n][i])

BC_median = [[] for i in range(len(SpT_range))]
BC_std = [[] for i in range(len(SpT_range))]
BC_size = [[] for i in range(len(SpT_range))]

for m in range(len(SpT_range)):
    for n in range(len(filters)):
        BC_median[m].append(np.ma.median(sigma_clip(
            BC_all[m][n], sigma=sigma, cenfunc='median', masked='False')))
        BC_std[m].append(np.std(sigma_clip(
            BC_all[m][n], sigma=sigma, cenfunc='median', masked='False')))
        BC_size[m].append(len(sigma_clip(
            BC_all[m][n], sigma=sigma, cenfunc='median', masked='False')))

BC_GG = []
BC_J = []
eBC_GG = []
eBC_J = []

for i in range(len(SpT_range)):
    BC_GG.append(BC_median[i][9])
    BC_J.append(BC_median[i][13])
    eBC_GG.append(BC_std[i][9])
    eBC_J.append(BC_std[i][13])

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


for m in range(len(SpT_range)):
    for n in range(len(Lbol_roll[m])):
        Radius_all[m].append(Radius_SB(Lbol_roll[m][n], Teff_roll[m][n]))
        Mass_all[m].append(Mass_sch19(Radius_all[m][n]))

Radius_median = []
Radius_std = []
Mass_median = []
Mass_std = []

for i in range(len(SpT_range)):
    Radius_median.append(np.median(sigma_clip(
        Radius_all[i], sigma=sigma, cenfunc='median', masked='False')))
    Radius_std.append(np.std(sigma_clip(
        Radius_all[i], sigma=sigma, cenfunc='median', masked='False')))
    Mass_median.append(np.median(sigma_clip(
        Mass_all[i], sigma=sigma, cenfunc='median', masked='False')))
    Mass_std.append(np.std(sigma_clip(
        Mass_all[i], sigma=sigma, cenfunc='median', masked='False')))

# %%

# WRITE OUT
# Uncomment (#) for Tables 6, 7 and 8, respectively.

filename = 'cif01.table_parameters.csv'
with open(filename, mode='w') as mycsv:
    writer = csv.writer(mycsv, delimiter='\n')
    header = []
    rows = []
    for i in range(len(SpT_name)-12):  # Maximum L2
        # Table 6 (#)
        rows.append(("{} & {:.3f} $\pm$ {:.3f} & {:.3f} $\pm$ {:.3f} & {:g} $\pm$ {:g} & {} $\pm$ {} & {:.3f} $\pm$ {:.3f}  & {:.3f} $\pm$ {:.3f} & {} \\\\").format(
            SpT_name[i], BC_GG[i], eBC_GG[i], BC_J[i], eBC_J[i],
            float('{:.3g}'.format(Lbol_median[i]*1e4)), float('{:.3g}'.format(
                Lbol_std[i]*1e4)), int(Teff_median[i]), int(Teff_std[i]),
            Radius_median[i], Radius_std[i],
            Mass_median[i], Mass_std[i],
            len(Lbol[i])))
        # Table 7/8 (#)
        # rows.append(("{} & {:.2f} $\pm$ {:.2f} & {:.2f} $\pm$ {:.2f} & {:.2f} $\pm$ {:.2f} & {:.2f} $\pm$ {:.2f} & {:.2f} $\pm$ {:.2f} & {:.2f} $\pm$ {:.2f} & {:.2f} $\pm$ {:.2f} \\\\").format(
        # SpT_name[i], Mabs_B[i], eMabs_B[i], Mabs_g[i], eMabs_g[i], Mabs_BP[i], eMabs_BP[i], Mabs_V[i], eMabs_V[i], Mabs_r[i], eMabs_r[i], Mabs_GG[i], eMabs_GG[i], Mabs_i[i], eMabs_i[i]))
        # SpT_name[i], Mabs_RP[i], eMabs_RP[i], Mabs_J[i], eMabs_J[i], Mabs_H[i], eMabs_H[i], Mabs_Ks[i], eMabs_Ks[i], Mabs_W1[i], eMabs_W1[i], Mabs_W2[i], eMabs_W2[i], Mabs_W3[i], eMabs_W3[i]))
    writer.writerow(rows)

# Reopen file and make-up / Copy into clipboard / Ready to paste in LaTeX table.

text = open(filename, 'r')
text = ''.join([i for i in text]).replace("nan $\pm$ nan", " \ldots ")
text = ''.join([i for i in text]).replace("0.0 $\pm$ nan", " \ldots ")
text = ''.join([i for i in text]).replace("-", "--")
x = open(filename, 'w')
x.writelines(text)
pyperclip.copy(text)
x.close()

# %%

# PLOTTING

# Data appending
# Literature: Pecaut+13 (E. Mamajek, priv. comm.)

with open('Literature/Mamajek_privcom.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_mamajek = []
    Teff_mamajek = []
    Lbol_mamajek = []
    R_mamajek = []
    for row in reader:
        SpT_mamajek.append(float(row['SpTnum']))
        Teff_mamajek.append(float(row['Teff']))
        Lbol_mamajek.append(float(row['logL']))
        R_mamajek.append(float(row['R_Rsun']))

Lbol_mamajek = [10**Lbol_mamajek[i] for i in range(len(Lbol_mamajek))]

# Literature: Mann+15

with open('Literature/mann15.stars.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT_man15 = []
    M_man15 = []
    for row in reader:
        try:
            SpT_man15.append(float(row['SpTnum']))
            M_man15.append(float(row['M']))
        except:
            pass

# Literature: Passegger+19

SpT_pas19 = [[] for i in range(len(SpT_range))]
TeffVN_pas19 = [[] for i in range(len(SpT_range))]
eTeffVN_pas19 = [[] for i in range(len(SpT_range))]

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    for row in reader:
        for m in range(len(SpT_range)):
            if float(row['SpTnum']) == SpT_range[m]:  # Choose Spectral type
                if row[booleans[0]] == row[booleans[1]] == 'false':
                    if row['Teff'] != '' and row['Teff_Pas19'] != '':
                        SpT_pas19[m].append(float(row['SpTnum']))
                        TeffVN_pas19[m].append(float(row['Teff_Pas19']))
                        eTeffVN_pas19[m].append(float(row['eTeff_Pas19']))

TeffVN_pas19_median = []
TeffVN_pas19_std = []
TeffVN_pas19_size = []

for i in range(len(SpT_range)):
    TeffVN_pas19_median.append(np.median(TeffVN_pas19[i]))
    TeffVN_pas19_std.append(np.std(TeffVN_pas19[i]))
    TeffVN_pas19_size.append(len(TeffVN_pas19[i]))

# Literature: Rajpurohit+18

SpT_raj18 = [[] for i in range(len(SpT_range))]
Teff_raj18 = [[] for i in range(len(SpT_range))]
eTeff_raj18 = [[] for i in range(len(SpT_range))]

with open('Literature/rajpurohit18.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    for row in reader:
        for m in range(len(SpT_range)):
            if float(row['SpTnum']) == SpT_range[m]:  # Choose Spectral type
                if row['Teff'] != '':
                    SpT_raj18[m].append(float(row['SpTnum']))
                    Teff_raj18[m].append(float(row['Teff']))

Teff_raj18_median = []
Teff_raj18_std = []
Teff_raj18_size = []

for i in range(len(SpT_range)):
    Teff_raj18_median.append(np.median(Teff_raj18[i]))
    Teff_raj18_std.append(np.std(Teff_raj18[i]))
    Teff_raj18_size.append(len(Teff_raj18[i]))


# Variables
# Uncomment the corresponding section (#) to plot

# Teff vs. SpT (#)

filename = 'cif20_plot_Teff_mean'
ylabel = r'$T_{\rm eff}\,[K]$'

x = [SpT_range[i] for i in range(len(SpT_range))]
y = [Teff_median[i] for i in range(len(SpT_range))]
z = [Teff_size[i] * alt_size for i in range(len(Teff_size))]

x_cif = SpT_cif
y_cif = Teff_cif
z_cif = SpT_cif

x_mamajek = SpT_mamajek
y_mamajek = Teff_mamajek
z_mamajek = SpT_mamajek

x_pas19 = SpT_range
y_pas19 = TeffVN_pas19_median
z_pas19 = TeffVN_pas19_size

x_raj18 = SpT_range
y_raj18 = Teff_raj18_median
z_raj18 = Teff_raj18_size

# Lbol vs. SpT (#)

# filename = 'cif20_plot_Lbol_mean'
# ylabel = r'$L\,[L_{\rm sol}]$'

# x = [SpT_range[i] for i in range(22)]  # only <=M9.5
# y = [Lbol_median[i] for i in range(22)]
# z = [Lbol_size[i] * alt_size for i in range(22)]

# x_cif = SpT_cif
# y_cif = Lbol_cif
# z_cif = SpT_cif

# x_mamajek = SpT_mamajek
# y_mamajek = Lbol_mamajek
# z_mamajek = SpT_mamajek

# R vs. SpT (#)

# filename = 'cif20_plot_radius_mean'
# ylabel = r'$\mathcal{R}$ [$\mathcal{R}_\odot$]'

# R_cif = [Radius_SB(Lbol_cif[i], Teff_cif[i]) for i in range(len(Lbol_cif))]
# M_cif = [Mass_sch19(R_cif[i]) for i in range(len(R_cif))]

# x = [SpT_range[i] for i in range(22)]  # only <=M9.5
# y = [Radius_median[i] for i in range(22)]
# z = [Lbol_size[i] for i in range(22)]

# x_cif = SpT_cif
# y_cif = R_cif
# z_cif = SpT_cif

# x_mamajek = SpT_mamajek
# y_mamajek = R_mamajek
# z_mamajek = SpT_mamajek

# x_man15 = SpT_man15
# y_man15 = M_man15


# Labels

xlabel = r'Spectral type'
cbarlabel = r'Spectral type'
SpT_name_ = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
             'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name_[2*i]
            for i in range(0, int(round(len(SpT_name_)/2, 0))+1)]

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

# Plots: all stars

sc = ax.scatter(x_cif, y_cif, c=z_cif, cmap=cmap,
                s=pointsize, marker='o', zorder=0)

# Plots: median values

plt.scatter(x, y, c='k', s=z, zorder=4)

# Plots: Pecaut+13

plt.plot(x_mamajek, y_mamajek, c='green', lw=2, zorder=1)
plt.scatter(x_mamajek, y_mamajek, c='green', s=pointsize*0.25, zorder=2)

# Plots: Mann et al. 2015 (#)
# Only R vs. SpT

# plt.scatter(x_man15, y_man15, edgecolors='blue',
#             facecolors='', s=pointsize, zorder=1)

# Plots: Passegger+19 and Rajopurohit+18 (#)
# Only Teff vs. SpT

plt.plot(x_raj18, y_raj18, c='red', lw=2, zorder=2)
plt.plot(x_pas19, y_pas19, c='blue', lw=2, zorder=3)

# for i in range(len(SpT_range)):
#     plt.scatter(x_raj18[i], y_raj18[i], c='red',
#                 s=z_raj18[i]*alt_size, zorder=2)
#     plt.scatter(x_pas19[i], y_pas19[i], c='blue',
#                 s=z_pas19[i]*alt_size, zorder=3)


# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

# Axes: ticks
# Use 1/2 to show all/half SpT.
# Use SpT_name/SpT_half to show all/half SpT.

ax.set_xticks(np.arange(-2, len(SpT_name_)-1, 2))
ax.set_xticklabels(SpT_half)
ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()
# ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)

# Axes: range & scale
# Limit on x-axis only.
# Uncomment depending on the variables represented.


ax.set_xlim(-2.5, 12.5)  # All
ax.set_ylim(1850, max(Teff_median)+100)  # Teff
# ax.set_ylim(1e-4, .4e0)  # Lbol
# ax.set_ylim(0.07, 1.2)  # R or M

# ax.set_xscale('log')
# ax.set_yscale('log')  # Lbol or R
# ax.invert_xaxis()
# ax.invert_yaxis()

# ax.set_yticks((0.1, 1))  # R or M
# ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
# ax.get_yaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())

# Colorbar
# Use 1/2 to show all/half SpT.
# Use SpT_name/SpT_half to show all/half SpT.

# divider = make_axes_locatable(plt.gca())
# cax = divider.append_axes("right", "2%", pad="1%")
# cbar = plt.colorbar(sc, cax=cax)  # Colorbar
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
# Only uncommnt this to keep colouring consistent.
sc.set_clim(vmin=-2, vmax=18)
# cbar.set_ticks(np.arange(-2, 19, 2))
# cbar.ax.set_yticklabels(SpT_half)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)


# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()


# %%
# REFERENCE
# Filter and spectral type number correspondence.

for i in range(len(filters)):
    print("Filter {}: {}".format(i, filters[i]))

for i in range(len(SpT_range)):
    print("SpTnum {}: {}".format(i, SpT_range[i]))
