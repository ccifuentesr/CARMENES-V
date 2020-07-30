
# %%

# PLOT: MODEL > ABSOLUTE MAGNITUDE - COLOUR (W/METALLICITY)
# Cifuentes et al. 2020
# MG vs. G - J

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

# Data
# Booleans allow to select only clean data.

Mother_version = '01'
filename = 'cif20_plot_MG_GJ_meta'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

x_axis = ['GG', 'J']
y_axis = 'GG'

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    SpT = []
    SpT_ = []
    d_pc = []
    d_pc_ = []
    mag0 = []
    mag0_ = []
    mag1 = []
    mag1_ = []
    mag2 = []
    mag2_ = []
    FeH = []
    for row in reader:
        if row[y_axis+'_mag'] != '' and row[x_axis[0]+'_mag'] != '' and row[x_axis[1]+'_mag'] != '':
            if row['d_pc'] != '':
                SpT_.append(float(row['SpTnum']))
                d_pc_.append(float(row['d_pc']))
                mag0_.append(float(row[x_axis[0]+'_mag']))
                mag1_.append(float(row[x_axis[1]+'_mag']))
                mag2_.append(float(row[y_axis+'_mag']))
                if row['Qf_'+y_axis] == row['Qf_'+x_axis[0]] == row['Qf_'+x_axis[1]] == 'false':
                    if row['FeH_lit'] != '':
                        if row[booleans[0]] == row[booleans[1]] == row[booleans[2]] == row[booleans[3]] == 'false':
                            SpT.append(float(row['SpTnum']))
                            d_pc.append(float(row['d_pc']))
                            mag0.append(float(row[x_axis[0]+'_mag']))
                            mag1.append(float(row[x_axis[1]+'_mag']))
                            mag2.append(float(row[y_axis+'_mag']))
                            FeH.append(float(row['FeH_lit']))

colour_x = [mag0[i] - mag1[i] for i in range(len(mag0))]
absmag = [mag2[i] - 5*np.log10(d_pc[i]) + 5 for i in range(len(mag2))]
colour_x_ = [mag0_[i] - mag1_[i] for i in range(len(mag0_))]
absmag_ = [mag2_[i] - 5*np.log10(d_pc_[i]) + 5 for i in range(len(mag2_))]

# Variables

x_cif_ = colour_x_
y_cif_ = absmag_

x_cif = colour_x
y_cif = absmag
z_cif = FeH

xfit = x_cif
yfit = y_cif
zfit = z_cif


# Statistics
# Choose degree k (#) for the polynomial fit.

k = 3

# Fit: Coefficients

# Fit data to polynomial and define a model based on it:
p, cov = np.polyfit(xfit, yfit, k, cov=True)
model = np.poly1d(p)
# Define custom data to draw the polynomial fit:
xfit_p = np.linspace(min(xfit), max(xfit), 1000)
yfit_p = model(xfit_p)
yfit_p_zero = [0 for i in range(len(xfit_p))]
# Standard-deviation estimates for each coefficient:
perr = np.sqrt(np.diag(cov))
# Coefficient of correlation between x and y:
R2 = np.corrcoef(yfit, model(xfit))[0, 1]**2
# Residuals:
resid = yfit - model(xfit)
resid_p = (yfit - model(xfit))/yfit
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

xlabel = r'$G-J$ [mag]'
ylabel = r'M$_G$ [mag]'
cbarlabel = r'[Fe/H]'

# Sizes

figsize = (12, 12)
pointsize = 30
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('copper')

# Plots: distribution

sc_ = ax.scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
                 edgecolors='gainsboro', zorder=0)
sc = ax.scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)

# Plot: fitting

ax.plot(xfit_p, yfit_p, 'b--', lw=2)

# Aesthetics

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
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

plt.xlim(1.3, 5.8)
ax.invert_yaxis()

# Colorbar

divider = make_axes_locatable(plt.gca())
cax = divider.append_axes("right", "2%", pad="1%")
cbar = plt.colorbar(sc, cax=cax)
sc.set_clim(vmin=-1.0, vmax=0.6)
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

# Plot: distribution

resid = np.array(resid).tolist()
sc0 = ax[0].scatter(x_cif, y_cif, s=pointsize, c=z_cif, cmap=cmap)
sc1 = ax[1].scatter(xfit, resid, s=pointsize, c=zfit, cmap=cmap, alpha=0.6)

# Plot: fitting

ax[0].scatter(x_cif_, y_cif_, s=pointsize, facecolors='',
              edgecolors='gainsboro', zorder=0)
ax[0].plot(xfit_p, yfit_p, 'b--', lw=2)
ax[1].plot(xfit_p, yfit_p_zero, '--', color='b', lw=2, zorder=3)

# Aesthetics

# Axes: labels & legend

ax[0].set_ylabel(ylabel, size=labelsize)
ax[1].set_ylabel('O-C', size=labelsize*0.8)
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

ax[0].set_xlim(1.3, 5.8)
ax[0].set_ylim(5.8, 21)
ax[1].set_ylim(-1.5, 1.5)
ax[0].invert_yaxis()

# Colorbar

sc0.set_clim(vmin=-1.0, vmax=0.6)
sc1.set_clim(vmin=-1.0, vmax=0.6)

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
