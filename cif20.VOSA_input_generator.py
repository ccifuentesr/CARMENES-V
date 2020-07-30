# %%

# PLOT: GENERATOR > VOSA INPUT
# Cifuentes et al. 2020

# Produces a file ready to be digested by VOSA.
# See http://svo2.cab.inta-csic.es/theory/vosa/.

import numpy as np
import pandas as pd
import pyperclip
import csv

# Data

Mother_version = '01'

filters = ['BT', 'B', 'g', 'BP', 'VT', 'V',
           'GG', 'r', 'i', 'RP', 'J', 'H', 'Ks', 'W1', 'W2', 'W3', 'W4']

VOSA_filters = ['TYCHO/TYCHO.B', 'Misc/UCAC.B', 'Misc/UCAC.sdss_g', 'GAIA/GAIA2r.Gbp', 'TYCHO/TYCHO.V',
                'Misc/UCAC.V', 'GAIA/GAIA2r.G', 'Misc/UCAC.sdss_r', 'Misc/UCAC.sdss_i', 'GAIA/GAIA2r.Grp', '2MASS/2MASS.J', '2MASS/2MASS.H',
                '2MASS/2MASS.Ks', 'WISE/WISE.W1', 'WISE/WISE.W2', 'WISE/WISE.W3', 'WISE/WISE.W4']

# Mother
# Uncomment to select G1, G2 or G3 group.
# G1: K5 V to K7 V.
# G2: M0.0 V to M9.5 V.
# G3: L0.0 to L8.0 objects.

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    Karmn = []
    RA = []
    DE = []
    d_pc = []
    ed_pc = []
    BT_mag = []
    B_mag = []
    mags = [[] for i in range(len(filters))]
    emags = [[] for i in range(len(filters))]
    for row in Mother:
        if row['Bool_delta'] == 'false' and row['Bool_dphot'] == 'false':
            if float(row['SpTnum']) <= 2.0:  # G1
                # if float(row['SpTnum']) > 2.0 and float(row['SpTnum']) <= 5.0:  # G2
                # if float(row['SpTnum']) > 5.0:  # G3
                Karmn.append(str(row['Karmn']))
                RA.append(float(row['RA_J2000']))
                DE.append(float(row['DE_J2000']))
                d_pc.append(str(row['d_pc']))
                ed_pc.append(str(row['ed_pc']))
                for i in range(len(filters)):
                    mags[i].append(str(row[filters[i]+'_mag']))
                    emags[i].append(str(row['e'+filters[i]+'_mag']))

# Write out
# Select the maximum size of the group.
# Write the group number for the output.

group_number = '1'

# To csv
filename = 'cif20.VOSA_G_'+group_number+'.txt'
with open(filename, mode='w') as mycsv:
    writer = csv.writer(mycsv, delimiter='\n')
    header = []
    rows = []
    for i in range(len(Karmn)):  # range(0, 1):
        for j in range(len(filters)):
            rows.append(("{} {} {} {}+-{} --- {} {} {} --- ---").format(
                Karmn[i], np.round(RA[i], 2), np.round(DE[i], 2), d_pc[i], ed_pc[i], VOSA_filters[j], mags[j][i], emags[j][i]))
    writer.writerow(rows)

# Reopen file and tune up.
# Copy in clipboard to paste directly in LaTeX table.

text = open(filename, 'r')
text = ''.join([i for i in text]).replace("   ---", " --- --- ---")
x = open(filename, 'w')
x.writelines(text)
# pyperclip.copy(text)  # copies into clipboard
x.close()
