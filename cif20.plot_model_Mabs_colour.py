# %%

# PLOT: MODEL > ABSOLUTE MAGNITUDE - COLOUR
# Cifuentes et al. 2020
# MG vs. G - J
# Mr vs. r - J

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
import scipy.stats as stats
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.
# Uncomment to choose MG or Mr (#)

Mother_version = '01'

# filename = 'cif20_plot_MG_GJ' #
# x_axis = ['GG', 'J']
# y_axis = 'GG'

filename = 'cif20_plot_Mr_rJ'
x_axis = ['r', 'J']
y_axis = 'r'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Karmn_ = []
    SpT = []
    SpT_ = []
    SpT_GTO = []
    d_pc = []
    d_pc_ = []
    d_pc_GTO = []
    mag0 = []
    mag0_ = []
    mag0_GTO = []
    mag1 = []
    mag1_ = []
    mag1_GTO = []
    mag2 = []
    mag2_ = []
    mag2_GTO = []
    for row in reader:
        if row[y_axis+'_mag'] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
            if row['d_pc'] != '':
                Karmn_.append(str(row['Karmn']))
                SpT_.append(float(row['SpTnum']))
                d_pc_.append(float(row['d_pc']))
                mag0_.append(float(row[x_axis[0]+'_mag']))
                mag1_.append(float(row[x_axis[1]+'_mag']))
                mag2_.append(float(row[y_axis+'_mag']))
                if row['Qf_'+y_axis] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        SpT.append(float(row['SpTnum']))
                        d_pc.append(float(row['d_pc']))
                        mag0.append(float(row[x_axis[0]+'_mag']))
                        mag1.append(float(row[x_axis[1]+'_mag']))
                        mag2.append(float(row[y_axis+'_mag']))
                        if row['Bool_GTO'] == 'true' and row['Young_YMG'] != 'Active' and row['Young_YMG'] != 'Active?':
                            SpT_GTO.append(float(row['SpTnum']))
                            d_pc_GTO.append(float(row['d_pc']))
                            mag0_GTO.append(float(row[x_axis[0]+'_mag']))
                            mag1_GTO.append(float(row[x_axis[1]+'_mag']))
                            mag2_GTO.append(float(row[y_axis+'_mag']))

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
Mabs = [mag2[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag2))]
colour_x_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]
Mabs_ = [mag2_[i] - 5*np.log10(d_pc_[i]) + 5 for i in range(len(mag2_))]
colour_x_GTO = [mag0_GTO[i] - mag1_GTO[i] for i in range(len(mag0_GTO))]
Mabs_GTO = [mag2_GTO[i] - 5 *
            np.log10(d_pc_GTO[i]) + 5 for i in range(len(mag2_GTO))]

# Variables

x_cif_ = colour_x_
y_cif_ = Mabs_

x_cif = colour_x
y_cif = Mabs
z_cif = SpT

x_cif_GTO = colour_x_GTO
y_cif_GTO = Mabs_GTO
z_cif_GTO = SpT_GTO

# Statistics
# Choose degree k (#) for the polynomial fit.

k = 3

# Fit: all-range
# Uncomment to fit all range.

# xfit = np.asarray(x_cif)
# yfit = np.asarray(y_cif)
# zfit = np.asarray(z_cif)

# Fit: truncated-range
# Uncomment to fit a limited range for MG or Mr.

# fit_range = [min(x_cif_GTO), 4.0]  # MG
fit_range = [min(x_cif_GTO), 5.2]  # Mr

xfit_bis = [i for i in x_cif_GTO if i >= fit_range[0] and i <= fit_range[1]]

indices = []
for i in range(len(xfit_bis)):
    indices.append(x_cif_GTO.index(xfit_bis[i]))

yfit_bis = []
zfit_bis = []
for i in indices:
    yfit_bis.append(y_cif_GTO[i])
    zfit_bis.append(z_cif_GTO[i])

xfit = np.asarray(xfit_bis)
yfit = np.asarray(yfit_bis)
zfit = np.asarray(zfit_bis)

# Fit: Coefficients

# Fit data to polynomial and define a model based on it:
p, cov = np.polyfit(xfit, yfit, k, cov=True)
model = np.poly1d(p)
# Define custom data to draw the polynomial fit:
xfit_ = np.sort(xfit)
yfit_ = model(xfit_)
yfit_zero = [0 for i in range(len(xfit_))]
# Standard-deviation estimates for each coefficient:
perr = np.sqrt(np.diag(cov))
# Coefficient of correlation between x and y:
R2 = np.corrcoef(yfit, model(xfit))[0, 1]**2
# Residuals:
resid = yfit - model(xfit)
# Degrees of freedom:
n = yfit.size
m = p.size
dof = n - m
# Chi-squared, reduced chi-squared, standard deviation of the error:
chi2 = np.sum((resid/model(xfit))**2)
chi2red = chi2/(dof)
s_err = np.sqrt(np.sum(resid**2)/(dof))
# Print out:
print('Polynomial fitting:\n')
print('degree =', k, '\ncoeffs a-c =', p, '\nerr_coeffs a-c =',
      perr, '\nR2 =', np.round(R2, 4), '\nchi2 =', np.round(chi2, 4), '\nchi2_red =', np.round(chi2red, 4))

# PLOTTING

# Labels
# Uncomment to choose MG or Mr (#)

# xlabel = r'$G-J$ [mag]' #
# ylabel = r'M$_G$ [mag]'
xlabel = r'$r-J$ [mag]'
ylabel = r'M$_r$ [mag]'

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

# Plot: distribution

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap, zorder=1)
sc_GTO = ax.scatter(xfit, yfit, s=pointsize, facecolors='',
                    edgecolors='dimgrey', zorder=2)

# Plot: fitting

ax.plot(xfit_, yfit_, 'b--', lw=2)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
# ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f')) # MG
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))  # Mr
ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

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

# plt.xlim(1.4, 5.8) # MG
# plt.ylim(5.4, 21) # MG
plt.xlim(1.2, 8.0)  # Mr
plt.ylim(5, 23)  # Mr

ax.invert_yaxis()

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)  # Colorbar
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
sc.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()


# %%
# PLOTTING (W/RESIDUALS)

# Canvas & Colours
fig, ax = plt.subplots(2, 1, sharex='col', gridspec_kw={
                       'hspace': 0.1, 'wspace': 0.4, 'height_ratios': [10, 2], 'hspace': 0.03}, figsize=figsize)

# Plot: distribution
ax[0].scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
              edgecolors='gainsboro', zorder=0)
sc0 = ax[0].scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap, zorder=1)
sc1 = ax[0].scatter(xfit, yfit, facecolors='', edgecolors='dimgrey', zorder=2)
sc2 = ax[1].scatter(xfit, resid, s=pointsize,
                    facecolors='', edgecolors='dimgrey')

# Plot: fitting
ax[0].plot(xfit_, yfit_, '--', color='b', lw=2, zorder=3)
ax[1].plot(xfit_, yfit_zero, '--', color='b', lw=2, zorder=3)

# Aesthetics

# Axes: labels & legend

ax[0].set_ylabel(ylabel, size=labelsize)
ax[1].set_ylabel('O-C\n[mag]', size=labelsize*0.8)
ax[1].set_xlabel(xlabel, size=labelsize)
# ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.1f')) # MG
ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))  # Mr

# Axes: ticks

ax[0].tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both', labelbottom=False)
ax[0].tick_params(axis='y', labelsize=tickssize, direction='in',
                  right=True, labelright=False, which='both')
ax[1].tick_params(axis='x', labelsize=tickssize, direction='in',
                  top=True, labeltop=False, which='both')
ax[1].tick_params(axis='y', labelsize=tickssize*0.8, direction='in',
                  right=True, labelright=False, which='both')
ax[0].tick_params('both', length=10, width=1, which='major')
ax[0].tick_params('both', length=5, width=1, which='minor')
ax[1].tick_params('both', length=10, width=1, which='major')
ax[1].tick_params('both', length=5, width=1, which='minor')
ax[0].xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax[0].minorticks_on()
ax[1].xaxis.set_tick_params(which='minor', bottom=True, top=True)

# Axes: range & scale

ax[1].set_xlim(1, 8)
ax[1].set_ylim(-1.5, 1.5)
ax[0].invert_yaxis()
ax[1].invert_yaxis()

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

cbar = fig.colorbar(sc0, ax=ax.ravel().tolist(), pad=0.01, aspect=50)
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=45)
sc0.set_clim(vmin=-2, vmax=18)
# sc1.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
# cbar.ax.invert_yaxis()
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()

# %%

# MINITOOL: PARSECATOR


def Distance(mag1, mag2):
    """Distance in parsec based on the absolute magnitude against colour model fitted above.

    Args:
        mag1 (float): First magnitude.
        mag2 (float): Second magnitude.

    Returns:
        float: Distance in parsec.

    The range of validity is defined by the range of fitted data.
    """
    d_pc = 0
    colour = mag1 - mag2
    if min(xfit) < colour < max(xfit):
        Mabs = model(colour)
        d_pc = 10**((mag1 + 5 - Mabs)/5)
    else:
        print('Error: colour out of the valid range.')
    return(round(d_pc, 1))


print('PARSECATOR\n')
print('-----------\n')
print('INPUT:', x_axis[0], 'and', x_axis[1], '\n')
print('OUTPUT: Distance in pc\n')
print('VALID:', x_axis[0], '-', x_axis[1], 'from', round(min(xfit), 1),
      'to', round(max(xfit), 1), 'mag\n')

# %%
# MINITOOL: PARSECATOR
# Run for automatic calculation for all stars in the sample.

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    Karmn_pre = []
    mag0_pre = []
    mag1_pre = []
    d_pc_pre = []
    for row in reader:
        if row[x_axis[0]+'_mag'] != '':
            if row[booleans[2]] == 'true':  # Only photometric distances
                Karmn_pre.append(str(row['Karmn']))
                d_pc_pre.append(str(row['d_pc']))
                mag0_pre.append(float(row[x_axis[0]+'_mag']))
                mag1_pre.append(float(row[x_axis[1]+'_mag']))
                colour_x_pre = [mag0_pre[i] - mag1_pre[i]
                                for i in range(len(mag0_pre))]

Karmn_est = []
colour_x_est = []
Mabs_est = []
mag0_est = []
d_pc_est = []

# Variables

for i in range(len(colour_x_pre)):
    if colour_x_pre[i] >= min(xfit) and colour_x_pre[i] <= max(xfit):  # valid range
        Karmn_est.append(Karmn_pre[i])
        mag0_est.append(mag0_pre[i])
        colour_x_est.append(colour_x_pre[i])
        Mabs_est.append(np.polyval(p, colour_x_pre[i]))

for i in range(len(Mabs_est)):
    d_pc_est.append(10**((5 - Mabs_est[i] + mag0_est[i])/5))

# Output to csv

with open('cif20.distance_estimated.csv', mode='w') as mycsv:
    writer = csv.writer(mycsv, delimiter=',')
    writer.writerow(['Karmn', 'd_pc', 'd_cif', 'r_J'])
    for i in range(len(Karmn_est)):
        writer.writerow([Karmn_est[i], d_pc_pre[i],
                         d_pc_est[i], colour_x_est[i]])
