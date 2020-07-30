# %%

# PLOT: MODEL > BOLOMETRIC CORRECTION - COLOUR
# Cifuentes et al. 2020
# BCG vs. G - J
# BCr vs. r - J
# BCJ vs. G - J
# BCW3 vs. G - J

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from scipy import stats
from scipy import optimize
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data
# Booleans allow to select only clean data.
# Uncomment to choose BCG, BCr, BCJ or BCW3 (#)

Mother_version = '01'
booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# y_axis = 'GG'
# x_axis = ['GG', 'J']
# filename = 'cif20_plot_BCG_GJ'

y_axis = 'r'
x_axis = ['r', 'J']
filename = 'cif20_plot_BCr_rJ'

# y_axis = 'J'
# x_axis = ['GG', 'J']
# filename = 'cif20_plot_BCJ_GJ'

# y_axis = 'W3'
# x_axis = ['GG', 'J']
# filename = 'cif20_plot_BCW3_GJ'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    d_pc = []
    d_pc_ = []
    Lbol = []
    Lbol_ = []
    mag0 = []
    mag0_ = []
    mag1 = []
    mag1_ = []
    mag2 = []
    mag2_ = []
    Teff = []
    Teff_ = []
    for row in reader:
        if row[y_axis+'_mag'] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
            if row['d_pc'] != '' and row['Lbol'] != '':
                SpT_.append(float(row['SpTnum']))
                d_pc_.append(float(row['d_pc']))
                Lbol_.append(float(row['Lbol']))
                mag0_.append(float(row[x_axis[0]+'_mag']))
                mag1_.append(float(row[x_axis[1]+'_mag']))
                mag2_.append(float(row[y_axis+'_mag']))
                Teff_.append(float(row['Teff']))
                if row['Qf_'+y_axis] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                    if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                        SpT.append(float(row['SpTnum']))
                        d_pc.append(float(row['d_pc']))
                        Lbol.append(float(row['Lbol']))
                        mag0.append(float(row[x_axis[0]+'_mag']))
                        mag1.append(float(row[x_axis[1]+'_mag']))
                        mag2.append(float(row[y_axis+'_mag']))
                        Teff.append(float(row['Teff']))


colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
colour_x_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]

Lsun = 3.828E26  # in Watts (IAU B2 Resolution)
Mbol = [71.197425-2.5*np.log10(Lbol[i]*Lsun) for i in range(len(Lbol))]
Mbol_ = [71.197425-2.5*np.log10(Lbol_[i]*Lsun) for i in range(len(Lbol_))]

Mabs = [mag2[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag2))]
Mabs_ = [mag2_[i] - 5*np.log10(d_pc_[i]) + 5 for i in range(len(mag2_))]

BC = [-Mabs[i] + Mbol[i] for i in range(len(Mabs))]
BC_ = [-Mabs_[i] + Mbol_[i] for i in range(len(Mabs_))]

# Variables

x_cif_ = colour_x_
y_cif_ = BC_

x_cif = colour_x
y_cif = BC
z_cif = SpT

# Statistics
# Choose degree k (#) for the polynomial fit.

k = 4

# Fit: all-range
# Uncomment to fit all range.
# Use for BCG, BCr.

xfit = x_cif
yfit = y_cif
zfit = z_cif

# Fit: truncated-range
# Uncomment to fit a limited range.
# Use for BCJ, BCW3.

# xfit_sortedbis = [i for i in x_cif if i <= 4]  # BCJ or BCW3
# indices = []
# for i in range(len(xfit_sortedbis)):
#     indices.append(x_cif.index(xfit_sortedbis[i]))

# yfit_sortedbis = []
# zfit_bis = []
# for i in indices:
#     yfit_sortedbis.append(y_cif[i])
#     zfit_bis.append(z_cif[i])

# xfit = np.asarray(xfit_sortedbis)
# yfit = np.asarray(yfit_sortedbis)
# zfit = np.asarray(zfit_bis)

# Fit: Coefficients

# Fit data to polynomial and define a model based on it:
p, cov = np.polyfit(xfit, yfit, k, cov=True)
model = np.poly1d(p)
# Define custom data to draw the polynomial fit:
# Uncomment (#) to use.
xr = np.linspace(np.min(xfit), np.max(xfit))
yr = model(xr)
yzero = [0 for i in range(len(xr))]
# Standard-deviation estimates for each coefficient:
perr = np.sqrt(np.diag(cov))
# Coefficient of correlation between x and y:
R2 = np.corrcoef(yfit, model(xfit))[0, 1]**2
# Residuals:
resid = yfit - model(xfit)
residp = (yfit - model(xfit))/yfit
# Degrees of freedom:
n = len(yfit)
m = p.size
dof = n - m  # degrees of freedom
# Chi-squared, reduced chi-squared, standard deviation of the error:
chi2 = np.sum((resid/yfit)**2)
chi2red = chi2/(dof)
s_err = np.sqrt(np.sum(resid**2)/(dof))
# Print out:
print('Polynomial fitting:\n')
print('degree =', k, '\ncoeffs a-c =', p, '\nerr_coeffs a-c =',
      perr, '\nR2 =', np.round(R2, 4), '\nchi2 =', np.round(chi2, 4), '\nchi2_red =', np.round(chi2red, 4))

# PLOTTING

# Labels
# Uncomment to choose BCG, BCr, BCJ or BCW3 (#)

# ylabel = r'BC$_G$ [mag]'
# xlabel = r'$G-J$ [mag]'

ylabel = r'BC$_r$ [mag]'
xlabel = r'$r-J$ [mag]'

# ylabel = r'BC$_J$ [mag]'
# xlabel = r'$G-J$ [mag]'

# ylabel = r'BC$_{W3}$ [mag]'
# xlabel = r'$G-J$ [mag]'

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
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

# Plot: fitting

ax.plot(xr, yr, 'b--', lw=2)

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
# Uncomment to choose BCG, BCr, BCJ or BCW3 (#)

# plt.xlim(1.3, 5.6)  # BCG
# plt.ylim(-4.1, 0.2)  # BCG

plt.xlim(1.0, 8.0)  # BCr
plt.ylim(-6.2, 0.2)  # BCr

# plt.xlim(1.3, 5.6)  # BCJ
# plt.ylim(1.35, 2.25)  # BCJ

# plt.xlim(1.2, 5.8)  # BCW3
# plt.ylim(1.9, 5.4)  # BCW3

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
sc0 = ax[0].scatter(xfit, yfit, s=pointsize, c=zfit, cmap=cmap)
sc1 = ax[1].scatter(xfit, resid, s=pointsize, c=zfit, cmap=cmap)
ax[0].scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
              edgecolors='gainsboro', zorder=0)

# Plot: fitting

ax[0].plot(xr, yr, 'b--', lw=2)
ax[1].plot(xr, yzero, 'b--', lw=2)
# ax[1].axhline(y=0, color='grey', linestyle='--', lw=2, zorder=0)
# ax[1].axhline(y=0.05, color='grey', linestyle='--', lw=2, zorder=0)
# ax[1].axhline(y=-0.05, color='grey', linestyle='--', lw=2, zorder=0)

# Aesthetics

# Plot: empty-faced
sc1.set_facecolor('none')

# Axes: labels & legend

ax[0].set_ylabel(ylabel, size=labelsize)
# ax[1].set_ylabel('(O-C)/O', size=labelsize*0.8)
ax[1].set_ylabel('O-C', size=labelsize*0.8)
ax[1].set_xlabel(xlabel, size=labelsize)
ax[0].xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

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

# ax[0].set_xlim(1.3, 5.6)  # BCG
# ax[0].set_ylim(-4.1, 0.2)  # BCG

ax[0].set_xlim(1.0, 8.0)   # BCr
ax[0].set_ylim(-6.2, 0.2)  # BCr

# ax[0].set_xlim(1.3, 5.6)  # BCG
# ax[0].set_ylim(1.35, 2.25)  # BCG

# ax[0].set_xlim(1.2, 5.8)  # BCW3
# ax[0].set_ylim(1.9, 5.4)  # BCW3

ax[1].set_ylim(0-0.25, 0+0.25)

# Colorbar

# Use 1/2 to show all/half SpT
# Use SpT_name/SpT_half to show all/half SpT

cbar = fig.colorbar(sc0, ax=ax.ravel().tolist(), pad=0.01, aspect=50)
cbar.set_label(cbarlabel, rotation=270, fontsize=labelsize, labelpad=45)
sc0.set_clim(vmin=-2, vmax=18)
sc1.set_clim(vmin=-2, vmax=18)
cbar.set_ticks(np.arange(-2, 19, 2))
cbar.ax.set_yticklabels(SpT_half)
cbar.ax.tick_params(labelsize=cblabsize)
cbar.outline.set_visible(False)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()

# %%

# MINITOOL: ILLUMINATOR (BC VERSION)


def Luminosity(mag1, mag2, d_pc):
    """Distance in parsec based on the bolometric correction against colour model fitted above.

    Args:
        mag1 (float): First magnitude.
        mag2 (float): Second magnitude.
        d_pc (float): Distance in parsec.

    Returns:
        float: Luminosity in solar units.

    The range of validity is defined by the range of fitted data.
    """
    Lbol = 0
    Lsol = 3.828E26  # in Watts (IAU B3)
    col = mag1 - mag2
    Mabs = mag1 - 5*np.log10(d_pc) + 5
    BC = model(col)
    Mbol = Mabs - BC
    if min(xfit) < col < max(xfit):
        Lbol = 10**((Mbol - 71.197425)/(-2.5))  # in Watts (IAU B3)
    else:
        print('Error: Bolometric correction out of range.')
    return(round(Lbol/Lsol, 6))
