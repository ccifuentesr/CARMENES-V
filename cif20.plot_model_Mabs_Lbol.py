# %%

# PLOT: MODEL > ABSOLUTE MAGNITUDE - LUMINOSITY
# Cifuentes et al. 2020
# MG vs. Lbol
# MJ vs. Lbol

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
# Uncomment to choose MG or MJ (#)

Mother_version = '01'

y_axis = 'Lbol'

# x_axis = 'GG'
# filename = 'cif20_plot_MG_Lbol'

x_axis = 'J'
filename = 'cif20_plot_MJ_Lbol'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    d_pc = []
    mag0 = []
    Lbol = []
    eLbol = []
    Teff = []
    for row in reader:
        if row[x_axis+'_mag'] != '' and row['Qf_'+x_axis] == 'false':
            if row['d_pc'] != '' and row['Teff'] != '' and row[y_axis] != '':
                if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                    SpT.append(float(row['SpTnum']))
                    d_pc.append(float(row['d_pc']))
                    Lbol.append(float(row[y_axis]))
                    eLbol.append(float(row['Lberr']))
                    Teff.append(float(row['Teff']))
                    mag0.append(float(row[x_axis+'_mag']))

Lbol_log = [np.log10(Lbol[i]) for i in range(len(Lbol))]
Lbol = [Lbol[i] for i in range(len(Lbol))]
Mabs = [mag0[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag0))]

# Vertical line

yp_v = np.linspace(7.5E-3-1, 7.5E-3+1, 100)
yp_vr = np.linspace(0-0.25, 0+0.25, 100)

xp_v = np.linspace(11.3, 11.3, 100)  # for J
# xp_v = np.linspace(14.0, 14.0, 100)  # for G

# Variables

x_cif = Mabs
y_cif = Lbol
z_cif = SpT

ey_cif = eLbol

# Fit: truncated-range 1
# Uncomment to fit a limited range for MG or MJ.

# xfit_1 = [i for i in x_cif if i <= 14.0 and i >= 0]  # MG
xfit_1 = [i for i in x_cif if i <= 11.3 and i >= 0]  # MJ
indices_1 = []
for i in range(len(xfit_1)):
    indices_1.append(x_cif.index(xfit_1[i]))
yfit_1 = []
for i in indices_1:
    yfit_1.append(y_cif[i])
zfit_1 = []
for i in indices_1:
    zfit_1.append(z_cif[i])
eyfit_1 = []
for i in indices_1:
    eyfit_1.append(ey_cif[i])

# Fit: truncated-range 2
# Uncomment to fit a limited range for MG or MJ.

# xfit_2 = [i for i in x_cif if i <= 21 and i > 14.1]  # MG
xfit_2 = [i for i in x_cif if i <= 20 and i > 11.3]  # MJ
indices_2 = []
for i in range(len(xfit_2)):
    indices_2.append(x_cif.index(xfit_2[i]))
yfit_2 = []
for i in indices_2:
    yfit_2.append(y_cif[i])
zfit_2 = []
for i in indices_2:
    zfit_2.append(z_cif[i])

# Statistics
# Two polynomials are fitted to the data in ranges 1 and 2.
# Choose degrees k_1 and k_2 for the two polynomial fits.

k_1 = 3
k_2 = 2

# Fit data to polynomial and define a model based on it:
p_1, cov_1 = np.polyfit(xfit_1, np.log10(yfit_1), k_1, cov=True)
p_2, cov_2 = np.polyfit(xfit_2, np.log10(yfit_2), k_2, cov=True)
model_1 = np.poly1d(p_1)
model_2 = np.poly1d(p_2)
# Define custom data to draw the polynomial fit:
xfit_1p = np.sort(xfit_1)
yfit_1p = 10**model_1(xfit_1p)  # Equivalent to np.polyval(p, xfit_1_)
yfit_1p_zero = [0 for i in range(len(xfit_1p))]
xfit_2p = np.sort(xfit_2)
yfit_2p = 10**model_2(xfit_2p)
yfit_2p_zero = [0 for i in range(len(xfit_2p))]
# Standard-deviation estimates for each coefficient:
perr_1 = np.sqrt(np.diag(cov_1))
perr_2 = np.sqrt(np.diag(cov_2))
# Coefficient of correlation between x and y:
R2_1 = np.corrcoef(yfit_1, 10**model_1(xfit_1))[0, 1]**2
R2_2 = np.corrcoef(yfit_2, 10**model_2(xfit_2))[0, 1]**2
# Residuals:
resid_1 = yfit_1 - 10**model_1(xfit_1)
resid_1p = (yfit_1 - 10**model_1(xfit_1))/yfit_1
resid_2 = yfit_2 - 10**model_2(xfit_2)
resid_2p = (yfit_2 - 10**model_2(xfit_2))/yfit_2
# Degrees of freedom:
n_1 = len(yfit_1)
m_1 = p_1.size
dof_1 = n_1 - m_1  # degrees of freedom
n_2 = len(yfit_2)
m_2 = p_2.size
dof_2 = n_2 - m_2
# Chi-squared, reduced chi-squared, standard deviation of the error:
chi2_1 = np.sum((resid_1/model_1(xfit_1))**2)
chi2red_1 = chi2_1/(dof_1)
s_err_1 = np.sqrt(np.sum(resid_1**2)/(dof_1))
chi2_2 = np.sum((resid_2/model_2(xfit_2))**2)
chi2red_2 = chi2_2/(dof_2)
s_err_2 = np.sqrt(np.sum(resid_2**2)/(dof_2))
# Print out:
print('First polynomial fitting:\n')
print('degree =', k_1, '\nR2 =', R2_1, '\nchi2 =', chi2_1, '\nchi2_red =',
      chi2red_1, '\ncoeffs a-c =', p_1, '\nerr_coeffs a-c =', perr_1)
print('\nSecond polynomial fitting:\n')
print('degree =', k_2, '\nR2 =', R2_2, '\nchi2 =', chi2_2, '\nchi2_red =',
      chi2red_2, '\ncoeffs a-c =', p_2, '\nerr_coeffs a-c =', perr_2)

# PLOTTING

# Labels
# Uncomment to choose MG or MJ (#)

# xlabel = r'$M_G$ [mag]' #
xlabel = r'$M_J$ [mag]'
ylabel = r'$L [L_{\rm sol}]$'
cbarlabel = r'Spectral type'
cbarlabel2 = r'T$_{\rm eff}$ [K]'
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

sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

# Plot: fitting

ax.semilogy(xfit_1p, yfit_1p, 'b--', lw=2)
ax.semilogy(xfit_2p, yfit_2p, 'b--', lw=2)

# Vertical line
ax.plot(xp_v, yp_v, '--', c='grey', lw=2)

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

plt.ylim(9E-6, 1)
ax.set_yscale('log')

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

# Plot: global

xfit_all = xfit_1 + xfit_2
zfit_all = zfit_1 + zfit_2
resid_all = np.array(resid_1).tolist() + np.array(resid_2).tolist()
resid_allp = np.array(resid_1p).tolist() + \
    np.array(resid_2p).tolist()  # p for percentage


# Plot: distribution

sc0 = ax[0].scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)
sc1 = ax[1].scatter(xfit_all, resid_allp, s=pointsize, c=zfit_all, cmap=cmap)

# Plot: fitting

ax[0].plot(xfit_1p, yfit_1p, 'b--', lw=2)
ax[0].plot(xfit_2p, yfit_2p, 'b--', lw=2)
ax[1].plot(xfit_1p, yfit_1p_zero, 'b--', lw=2)
ax[1].plot(xfit_2p, yfit_2p_zero, 'b--', lw=2)

# Vertical line

ax[0].plot(xp_v, yp_v, '--', c='grey', lw=2)
ax[1].plot(xp_v, yp_vr, '--', c='grey', lw=2)

# Aesthetics

# Plot: empty-faced

sc1.set_facecolor('none')

# Axes: labels & legend

ax[0].set_ylabel(ylabel, size=labelsize)
ax[1].set_ylabel('(O-C)/C', size=labelsize*0.8)
# ax[1].set_ylabel('O-C', size=labelsize*0.8)
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


ax[0].set_ylim([7E-6, 1])
# ax[1].set_ylim([0-1.5E-2, 0+1.5E-2])
ax[1].set_ylim([-.2, .2])  # percentage
ax[0].set_yscale('log')
# ax[1].set_yscale('log')
# ax.invert_xaxis()
# ax.invert_yaxis()

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

cbar = fig.colorbar(sc0, ax=ax.ravel().tolist(), pad=0.01, aspect=50)
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=45)
sc1.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
# cbar.ax.invert_yaxis()
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()

# %%

# MODEL TESTING
# Compares model predictions with VOSA values.

# Data
# Defined previously.


ylabel = r'L$_{\rm VOSA}$/L$_J$ [L$_\odot]$'
xlabel = r'L$_{\rm VOSA}$ [L$_\odot]$'

filename = 'cif20_plot_L_model'

# Data: Model


def Luminosity(J_mag, d_pc):
    Lbol = 1
    MJ = J_mag - 5*np.log10(d_pc) + 5
    if min(xfit_1) < MJ < max(xfit_1):
        Lbol = 10**model_1(MJ)
    if min(xfit_2) < MJ < max(xfit_2):
        Lbol = 10**model_2(MJ)
    return(Lbol)

# Data: Variables


Lbol_estimated = []
for i in range(len(Lbol)):
    Lbol_estimated.append(Luminosity(mag0[i], d_pc[i]))

Lbol_log = [np.log10(Lbol[i]) for i in range(len(Lbol))]
Lbol_ratio = [Lbol[i]/Lbol_estimated[i] for i in range(len(Lbol))]

Delta_Lbol = [abs((Lbol_estimated[i] - Lbol[i])/Lbol[i])
              for i in range(len(Lbol))]

x = Lbol
y = Lbol_ratio
z = SpT
xp = np.logspace(-6, 0, 1000)
yp = xp

# PLOTTING

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('magma_r')

# Plot: distributions

sc0 = ax.scatter(x, y, c=z, cmap=cmap, s=pointsize,
                 marker='o', label='This work')

# Plot: lines

ax.axhline(y=np.mean(Lbol_ratio), color='grey', linestyle='--', lw=2, zorder=0)
ax.axhline(y=np.mean(Lbol_ratio)+np.std(Lbol_ratio),
           color='grey', linestyle='--', lw=2, zorder=0)
ax.axhline(y=np.mean(Lbol_ratio)-np.std(Lbol_ratio),
           color='grey', linestyle='--', lw=2, zorder=0)

# Aesthetics

# Plot: empty-faced

sc0.set_facecolor('none')
sc1.set_facecolor('none')

# Axes: range & scale

ax.set_xlim(1e-5, 1)
ax.set_ylim(1-.3, 1+.3)
ax.set_xscale('log')

# Axes: labels & legend

ax.set_ylabel(ylabel, size=labelsize)
ax.set_xlabel(xlabel, size=labelsize)

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both', labelbottom=True)
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=True, labelright=False, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax.xaxis.set_tick_params(which='minor', bottom=True, top=True)
ax.minorticks_on()

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
# cbar.ax.invert_yaxis()
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()

# %%
# MINITOOL: ILLUMINATOR


def Luminosity(mag, d_pc):
    """Bolometric luminosity in solar units from absolute magnitude in J (preferred) or G passbands.
    Args:
        mag (float): J or G magnitude.
        d_pc (float): Distance in parsec.

    Returns:
        float: Luminosity in solar units.

    The range of validity is defined by the range of fitted data.
    """
    Lbol = 0
    Mabs = mag - 5*np.log10(d_pc) + 5
    if min(xfit_1) < Mabs < max(xfit_1):
        Lbol = 10**model_1(Mabs)
    elif min(xfit_2) < Mabs < max(xfit_2):
        Lbol = 10**model_2(Mabs)
    else:
        print('Error: Absolute magnitude out of range.')
    return(round(Lbol, 6))


print('ILLUMINATOR\n')
print('-----------\n')
print('INPUT:', x_axis[0], 'and d_pc\n')
print('OUTPUT: Luminosity in Lsol\n')
print('VALID: M', x_axis[0], 'from', round(min(xfit_1), 1),
      'to', round(max(xfit_2), 1), 'mag\n')
