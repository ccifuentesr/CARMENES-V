# %%

# GENERATOR: AVERAGE COLOURS
# Cifuentes et al. 2020 (Table A.2)
# Average colours for K5V to L8 sources.

# Read https://stackoverflow.com/questions/34050491/standard-deviation-in-numpy
# for 'ddof=1' in numpy std.
# Output ready for pasting into LaTeX tabular environment.

import numpy as np
import csv
import pyperclip

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
booleans = ['Bool_delta', 'Bool_young']

# Choose filter set and spectral type (i.e., a row of the table)

# filters = ['FUV', 'NUV', 'u', 'BT', 'B', 'g', 'BP', 'VT', 'V', 'GG', 'r']
filters = ['r', 'i', 'RP', 'J', 'H', 'Ks', 'W1', 'W2', 'W3', 'W4']

# SpT = [0]
SpT_range_pre = np.arange(0, 18.5, 0.5)  # M0.0 to L8.0
SpT_range = [-2, -1]  # K5 and K7
for i in range(len(SpT_range_pre)):
    SpT_range.append(np.ndarray.tolist(
        np.arange(0, 18.5, 0.5))[i])  # all range

SpT_name = ['K5\,V', 'K7\,V', 'M0.0\,V', 'M0.5\,V', 'M1.0\,V', 'M1.5\,V', 'M2.0\,V', 'M2.5\,V', 'M3.0\,V', 'M3.5\,V', 'M4.0\,V', 'M4.5\,V', 'M5.0\,V', 'M5.5\,V', 'M6.0\,V', 'M6.5\,V', 'M7.0\,V',
            'M7.5\,V', 'M8.0\,V', 'M8.5\,V', 'M9.0\,V', 'M9.5\,V', 'L0', 'L0.5', 'L1', 'L1.5', 'L2', 'L2.5', 'L3', 'L3.5', 'L4', 'L4.5', 'L5', 'L5.5', 'L6', 'L6.5', 'L7', 'L7.5', 'L8']

# Variables

colours = [[[] for i in range(len(filters))] for j in range(len(SpT_range))]
# mags = [[] for _ in range(len(filters))]  # a column for each filter
mean = [[] for i in range(len(SpT_range))]
stdev = [[] for i in range(len(SpT_range))]
items = [[] for i in range(len(SpT_range))]

# Mother: mean values

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    Mother = csv.DictReader(mycsv)
    for row in Mother:
        for m in range(len(SpT_range)-1):
            for n in range(len(filters)-1):
                if float(row['SpTnum']) == SpT_range[m]:  # Choose Spectral type
                    # Only clean stars
                    if row[filters[n]+'_mag'] != '' and row[filters[n+1]+'_mag'] != '':
                        if row[booleans[0]] == row[booleans[1]] == 'false':
                            if row['Qf_'+filters[n]] == 'false' and row['Qf_'+filters[n+1]] == 'false':
                                colours[m][n].append(float(row[filters[n]+'_mag']) -
                                                     float(row[filters[n+1]+'_mag']))


for m in range(len(SpT_range)):
    for n in range(len(filters)-1):
        mean[m].append(round(np.mean(colours[m][n]), 2))
        stdev[m].append(round(np.std(colours[m][n]), 2))
        items[m].append(len(colours[m][n]))

# WRITE OUT

# To csv
filename = 'cif01.average_colours.csv'
with open(filename, mode='w') as mycsv:
    writer = csv.writer(mycsv, delimiter=' ')
    header = []
    rows = []
    for m in range(len(SpT_range)):
        rows.append(SpT_name[m])
        for n in range(len(mean[m])):
            rows.append(("& {:.2f} $\pm$ {:.2f} ({})").format(
                mean[m][n], stdev[m][n], items[m][n]))
        rows.append("\\\\")
    writer.writerow(rows)

# Reopen file and make-up / Copy into clipboard / Ready to paste in LaTeX table.

text = open('cif01.average_colours.csv', 'r')
text = ''.join([i for i in text]).replace("&", " & ")
text = ''.join([i for i in text]).replace("$\pm$ 0.00", " ")
text = ''.join([i for i in text]).replace("nan $\pm$ nan (0)", "\ldots")
x = open('cif01.average_colours.csv', 'w')
x.writelines(text)
pyperclip.copy(text)
x.close()
