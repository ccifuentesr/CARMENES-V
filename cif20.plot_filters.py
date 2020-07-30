# %%

# PLOT: FILTERS
# Cifuentes et al. 2020
# Filters information from VOSA's Filter Profile Service.
# http://svo2.cab.inta-csic.es/theory/fps/

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
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
filename = 'cif20_plot_filters'
starname = 'J11026+219_3700_50.csv'

booleans = ['Bool_delta', 'Bool_young', 'Bool_dphot', 'Bool_RUWE']

# Star

filters = ['FUV', 'NUV', 'u', 'BT', 'B', 'g', 'BP', 'VT', 'V',
           'GG', 'r', 'i', 'RP', 'J', 'H', 'Ks', 'W1', 'W2', 'W3', 'W4']

mags = [[] for x in range(len(filters))]
emags = [[] for x in range(len(filters))]

# Mother

with open('cif20.Mother.v'+Mother_version+'.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    for row in reader:
        if row['Karmn'] == 'J11026+219':
            for i in range(len(filters)):
                mags[i].append(float(row[filters[i]+'_mag']))
                emags[i].append(float(row['e'+filters[i]+'_mag']))

# Filters
# Normalisation

filters_ = [[] for x in range(len(filters))]
filters_lambda = [[] for x in range(len(filters))]
filters_transm = [[] for x in range(len(filters))]
filters_transm_norm = [[] for x in range(len(filters))]

for n in range(len(filters)):
    with open('Filters/Passband_'+filters[n]+'.dat', 'r') as fobj:
        for line in fobj:
            for num in line.split():
                filters_[n].append(float(num))

for n in range(len(filters)):
    for k in range(0, int(len(filters_[n])/2)):
        filters_lambda[n].append(filters_[n][2*k])
        filters_transm[n].append(filters_[n][abs(2*k+1)])

for n in range(len(filters)):
    filters_transm_norm[n].append(
        np.array(filters_transm[n])/max(filters_transm[n]))

# Filters information (VOSA)
# columns: Filter,Phot_Sys,Lambda_eff,VEGA_erg,VEGA_J

with open('Filters/VOSA_filters.csv', 'r') as mycsv:
    reader = csv.DictReader(mycsv)
    l = []  # lambdas
    z = []  # zero fluxes
    for row in reader:
        l.append(float(row['lambda_eff']))
        z.append(float(row['VEGA_erg']))

# erg/s/cm2/A to W/m2/A
fluxes = [z[i]*10**(mags[i][0]/-2.5)*1E-3 for i in range(len(z))]
lfluxes = [fluxes[i]*l[i] for i in range(len(l))]
fluxes_norm = [(fluxes[i]/max(fluxes)) for i in range(len(fluxes))]
lfluxes = np.array(l)*np.array(fluxes)  # W/m2
lfluxes_norm = [(lfluxes[i]/max(lfluxes)) for i in range(len(lfluxes))]

# Star (VOSA output)
# From BT to W4

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
lambda_obs_flux_norm = [lambda_obs_flux[i] /
                        max(lambda_obs_flux) for i in range(len(lambda_obs_flux))]

# PLOTTING

# Labels

xlabel = r'$\lambda$ ($\AA$)'
ylabel = r'Transmitance (normalised)'
SpT_name = ['K5', 'K7', 'M0', 'M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7',
            'M8', 'M9', 'L0', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8']
SpT_half = [SpT_name[2*i] for i in range(0, int(round(len(SpT_name)/2, 0))+1)]

# Sizes

figsize = (16, 10)
pointsize = 40
tickssize = 22
labelsize = 22
legendsize = 18
cblabsize = 18

# Canvas & Colours

fig, ax = plt.subplots(figsize=figsize)
ax_ = ax.twinx()

cmap = plt.get_cmap('rainbow')

# Plot: filter curves

for n in range(0, len(filters)):
    ax.plot(filters_lambda[n], filters_transm_norm[n][0], c=cmap(
        n/len(filters)), linewidth=1.2, label=filters[n], zorder=1)

# Plot: SED

for n in range(len(l)):  # Bluer than BT
    ax_.scatter(l[n], lfluxes[n], s=pointsize*2, c=cmap(n/len(l)))

# Aesthetics

# Axes: range & scale

ax.set_xlim(1e3, 3e5)
ax.set_ylim(-0.01, 1.05)
ax_.set_ylim(min(lfluxes)*0.85, max(lfluxes)*1.5)
ax.set_xscale('log')
ax_.set_yscale('log')

# Axes: labels & legend

ax.set_xlabel(xlabel, size=labelsize)
ax.set_ylabel(ylabel, size=labelsize)
ax_.set_ylabel(r'$\lambda F_\lambda$ (W/m$^2$)', size=labelsize)

# Axes: ticks

ax.tick_params(axis='x', labelsize=tickssize, direction='in',
               top=True, labeltop=False, which='both')
ax.tick_params(axis='y', labelsize=tickssize, direction='in',
               right=False, labelright=False, which='both')
ax_.tick_params(axis='y', labelsize=tickssize, direction='in',
                right=True, labelright=True, which='both')
ax.tick_params('both', length=10, width=1, which='major')
ax.tick_params('both', length=5, width=1, which='minor')
ax_.tick_params('both', length=10, width=1, which='major')
ax_.tick_params('both', length=5, width=1, which='minor')
ax.minorticks_on()

# Show & Save

plt.savefig(filename+'.eps', bbox_inches='tight')
plt.show()
