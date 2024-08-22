# %%

# PLOT: SED
# Cifuentes et al. 2020+


from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import FormatStrFormatter
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import csv

rc('font', **{'family': 'serif', 'serif': ['DejaVu Serif']})
rc('text', usetex=False)
plt.rcParams['mathtext.fontset'] = 'dejavuserif'

# Data

starname = 'J10384+485_3300_55.csv'
spectraname = 'lte033.0-5.5-0.0a+0.0.BT-Settl.spec.7.dat'
filtername = 'Filters_EW_BT.csv'
filename = 'cif20_plot_SED'

# Filters

with open('Filters/'+filtername, 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    filter_wavelength = []
    filter_weff = []
    for row in reader:
        filter_wavelength.append(float(row['Wavelength']))
        filter_weff.append(float(row['W_eff']))

# Star (VOSA output)

with open('Stars/'+starname, 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    wavelength = []
    w_eff = []
    obs_flux = []
    obs_error = []
    model_flux = []
    for row in reader:
        wavelength.append(float(row['Wavelength']))
        obs_flux.append(float(row['Obs.Flux']))
        obs_error.append(float(row['Obs.Error']))
        model_flux.append(float(row['FluxMod']))

lambda_obs_flux = np.array(wavelength)*np.array(obs_flux)*1E-7*1e4  # J/s/m2
lambda_obs_error = np.array(wavelength)*np.array(obs_error)*1E-7*1e4  # J/s/m2
lambda_model_flux = np.array(
    wavelength)*np.array(model_flux)*1E-7*1e4  # J/s/m2

# Spectra

with open('Spectra/'+spectraname, 'r') as mycsv:
    dat_reader = csv.reader(mycsv)
    dat_ref = list(dat_reader)
    lref = []  # lambdas
    fref = []  # fluxes
    for i in range(0, len(dat_ref)-1):
        lref.append(float(dat_ref[i][0].split()[0]))  # Angstrom
        fref.append(float(dat_ref[i][0].split()[1]))  # erg/cm2/s/A

l_fref = list(np.array(lref)*np.array(fref))  # lambda x flux model
l_fref_norm = np.array(l_fref)/max(l_fref) * \
    max(lambda_model_flux)  # normalise to J/s/cm2

# Add UV passbands (FUV, NUV, u)

uv_mags = [21.832, 19.88, 16.787]
uv_lambda = [1549.0, 2304.7, 3594.9]
uv_zeropoints = [6.506e-9*1E-7*1e4, 4.450e-9 *
                 1E-7*1e4, 3.639e-9*1E-7*1e4]  # J m2 s A
uv_fluxes = [uv_zeropoints[n]*10**(-uv_mags[n]/2.5)
             for n in range(len(uv_mags))]
uv_lfluxes = [uv_lambda[n]*uv_fluxes[n] for n in range(len(uv_mags))]

# PLOTTING

# Labels

xlabel = r'$\lambda$ [$\AA$]'
ylabel = r'$\lambda F_\lambda \,[{\rm W/m}^{2}]$'

# Sizes

figsize = (16, 10)
pointsize = 60
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
cmap = plt.get_cmap('rainbow')

# Plots: theoretical spectrum

ax.plot(lref, l_fref_norm, c='gainsboro', linewidth=1.2, zorder=0)

# Plots: SED

for n in range(len(wavelength)):
    ax.errorbar(wavelength[n], lambda_obs_flux[n], yerr=lambda_obs_error[n], xerr=filter_weff[n]/2, ms=10, fmt='o', color=cmap(
        (n+3)/len(wavelength)), markerfacecolor='none', ecolor=cmap((n+3)/len(wavelength)), capthick=2, zorder=2)
    ax.scatter(wavelength[n], lambda_model_flux[n],
               s=pointsize, marker='o', edgecolors='black', color='none', zorder=1)

# Plots: UV passbands

for n in range(len(uv_mags)):
    ax.scatter(uv_lambda[n], uv_lfluxes[n], s=pointsize,
               marker='x', color=cmap((n)/len(wavelength)), zorder=1)

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

plt.xlim(1.3e3, 4e5)
plt.ylim(min(lambda_obs_flux)*0.01, max(lambda_obs_flux)*2.0)
ax.set_xscale('log')
ax.set_yscale('log')

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
