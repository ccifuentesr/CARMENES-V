# %%

# CALCULATOR: MASSES & RADII
# Cifuentes et al. 2020

# Stefan-Boltzmann L-R
# Schweitzer et al. 2019 M-R

import numpy as np

Mother_version = '01'

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

# %%

# WRITE OUT


Karmn = []
Lbol = []
Lberr = []
Teff = []

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    for row in Mother:
        if row['Lbol'] != '':
            Karmn.append(str(row['Karmn']))
            Lbol.append(float(row['Lbol']))
            Lberr.append(float(row['Lberr']))
            Teff.append(float(row['Teff']))

eTeff = [50 for i in range(len(Teff))]

Radius = []
Mass = []
for i in range(len(Lbol)):
    Radius.append(Radius_SB(Lbol[i], Lberr[i], Teff[i], eTeff[i]))
    Mass.append(Mass_sch19(Radius[i][0], Radius[i][1]))

filename = 'cif20.Mother.v'+Mother_version+'_MR.csv'
with open(filename, mode='w') as mycsv:
    writer = csv.writer(mycsv, delimiter='\n')
    header = []
    rows = []
    rows.append("Karmn,Radius,eRadius,Mass,eMass")
    for i in range(len(Lbol)):
        rows.append(("{},{:.4f},{:.4f},{:.4f},{:.4f}").format(
            Karmn[i], Radius[i][0], Radius[i][1], Mass[i][0], Mass[i][1]))
    writer.writerow(rows)
