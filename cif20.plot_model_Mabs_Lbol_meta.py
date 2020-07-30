# %%

# PLOT: MODEL > ABSOLUTE MAGNITUDE - LUMINOSITY (W/METALLICITY)
# Cifuentes et al. 2020
# MG vs. G - J

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import csv
from scipy.optimize import curve_fit
from scipy.stats import pearsonr, spearmanr
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
filename = 'cif20_plot_MJ_Lbol_meta'

y_axis = 'Lbol'
x_axis = 'J'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    d_pc = []
    d_pc_ = []
    mag0 = []
    mag0_ = []
    Lbol = []
    Lbol_ = []
    eLbol = []
    eLbol_ = []
    Teff = []
    Teff_ = []
    Meta_lit = []
    Meta_lit_ = []
    eMeta_lit = []
    eMeta_lit_ = []
    for row in reader:
        if row[x_axis+'_mag'] != '' and row['Qf_'+x_axis] == 'false':
            if row['d_pc'] != '' and row['Teff'] != '':
                if row[y_axis] != '':
                    SpT_.append(float(row['SpTnum']))
                    d_pc_.append(float(row['d_pc']))
                    Lbol_.append(float(row[y_axis]))
                    eLbol_.append(float(row['Lberr']))
                    Teff_.append(float(row['Teff']))
                    mag0_.append(float(row[x_axis+'_mag']))
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        if row['FeH_lit'] != '':
                            SpT.append(float(row['SpTnum']))
                            d_pc.append(float(row['d_pc']))
                            Lbol.append(float(row[y_axis]))
                            eLbol.append(float(row['Lberr']))
                            Teff.append(float(row['Teff']))
                            mag0.append(float(row[x_axis+'_mag']))
                            Meta_lit.append(float(row['FeH_lit']))
                            eMeta_lit.append(float(row['eFeH_lit']))

Lbol_log = [np.log10(Lbol[i]) for i in range(len(Lbol))]
Lbol_log_ = [np.log10(Lbol_[i]) for i in range(len(Lbol_))]
Lbol = [Lbol[i] for i in range(len(Lbol))]
Lbol_ = [Lbol_[i] for i in range(len(Lbol_))]
Mabs = [mag0[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag0))]
Mabs_ = [mag0_[i] - 5*np.log10(d_pc_[i]) + 5 for i in range(len(mag0_))]

# Variables

x_cif = Mabs
y_cif = Lbol
z_cif = Meta_lit
ey_cif = eLbol

x_cif_ = Mabs_
y_cif_ = Lbol_
ey_cif_ = eLbol_

xfit = x_cif
yfit = y_cif
zfit = z_cif

# Statistics
# Choose degree k (#) for the polynomial fit.

k = 3

# Fit data to polynomial and define a model based on it:
p, cov = np.polyfit(xfit, np.log10(
    yfit), k, cov=True)
model = np.poly1d(p)
# Define custom data to draw the polynomial fit:
xfit_p = np.sort(xfit)
yfit_p = 10**model(xfit_p)
yfit_p_zero = [0 for i in range(len(xfit_p))]
# Standard-deviation estimates for each coefficient:
perr = np.sqrt(np.diag(cov))
# Coefficient of correlation between x and y:
R2 = np.corrcoef(yfit, 10**model(xfit))[0, 1]**2
# Residuals:
resid = yfit - 10**model(xfit)
resid_p = (yfit - 10**model(xfit))/yfit
# Degrees of freedom:
n = len(yfit)
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

xlabel = r'$M_J$ [mag]'
ylabel = r'$L [L_{\rm sol}]$'
cbarlabel = r'[Fe/H]'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (12, 12)
pointsize = 40
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('copper')

# Plot: distribution

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='none',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

# Plot: fitting

ax.semilogy(xfit_p, yfit_p, 'b--', lw=2)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
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

ax.set_xlim(5, 11.1)
ax.set_ylim(3E-4, 2E-1)
ax.set_yscale('log')

# Colorbar

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)
sc.set_clim(vmin=-0.3, vmax=0.5)
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=30)
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

# Plot: global

resid = np.array(resid).tolist()

# Plot: distribution

sc0 = ax[0].scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)
sc1 = ax[1].scatter(xfit, resid_p, s=pointsize, c=zfit, cmap=cmap, alpha=0.6)

# Plot: fitting

ax[0].semilogy(xfit_p, yfit_p, 'b--', lw=2)
ax[1].plot(xfit_p, yfit_p_zero, 'b--', lw=2)

# Aesthetics

# Axes: labels & legend

ax[0].set_ylabel(ylabel, size=labelsize)
ax[1].set_ylabel('(O-C)/C', size=labelsize*0.8)
ax[1].set_xlabel(xlabel, size=labelsize)
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

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
ax[1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# Axes: range & scale

ax[0].set_xlim(5, 11.1)
ax[0].set_ylim(3E-4, 2E-1)

ax[1].set_ylim([-.2, .2])  # percentage
ax[0].set_yscale('log')

# Colorbar
# Only keep (#) to hide colorbar with consistent colouring.

# cbar = fig.colorbar(sc0, ax=ax.ravel().tolist(), pad=0.01, aspect=50)
# cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=45)
sc0.set_clim(vmin=-1.0, vmax=0.6)
sc1.set_clim(vmin=-1.0, vmax=0.6)
# cbar.ax.tick_params(labelsize=cblabsize)
# cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
