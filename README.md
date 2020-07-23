# CARMENES input catalogue of M dwarfs. 
## V. Luminosities, colours, and spectral energy distributions.
  
> This repository contains the pieces of code necessary to produce all figures, tables and models in Cifuentes et al. 2020. 

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Publication](https://img.shields.io/badge/Published%3F-waiting-orange.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://GitHub.com/ccifuentesr)

## Table of Contents

- [Installation](#installation)
- [Structure](#structure)
- [Support](#support)
- [License](#license)
- [Suggested Resources](#resources)

---

## Installation

> The files are self-contained, self-consistent, homogenoeusly formatted, fairly self-explanatory.

[![Usage](https://github.com/ccifuentesr/CARMENES-V/blob/master/code.png)]()

- The code is provided as `*.py` files meant to be run individually.
- They can be run as Python Notebooks. The symbol `# %%` starts a cell that can be run separately.
- Some files requires of additional data contained in the folders stored in the repository.
- The installation of basic libraries such as `numpy` is a prerequisite. 

### Clone

- Clone this repo to your local machine using `git clone https://github.com/ccifuentesr/CARMENES-V`.

## Structure

### Directories

- Directory ./: Includes all the code files classified as detailed below and the master table (`Mother.v01.csv`).
- Directory ./Baraffe: (see https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/).
- Directory ./Filters: Includes the transmission curves for 20 filters (FUV to W4). 
- Directory ./Literature: Includes the individual tables with data from the literature.
- Directory ./Spectra: Includes PHOENIX synthetic spectra.
- Directory ./Stars: Includes data for individual stars to produce custom SEDs.

### Files

All files are named `cif20.xxx_yyy_zzz.py`, where `xxx` define the kind of output that it produces, `yyy` gives additional information about the output, and `zzz` enumerates the main variables involved. For example, the script `cif20.plot_literature_Mabs_SpT.py` produces an absolute magnitude vs. spectral type plot, and compares the values with those of the literature. The complete list of files and their description can be found below.

| File | Description | Input<sup id="a1">[1](#f1)</sup>| Output | 
| --- | --- | --- | --- | 
| cif20.calculator_averagecolors.py | Average colours and standard deviations for each spectral type from K5V to L8. | Magnitudes from FUV to W4 | Table A.2 |
| cif20.calculator_Karmn.py | A generator of Carmencita identification name (Karmn JHHMMm+DDdAAA). | Simbad valid name | Karmn ID |
|	cif20.calculator_MR.py	|	Calculator of mass and radius in solar units.	|	Lbol, Teff	|	... |	
|	cif20.calculator_completeness.py |	Estimator of the completeness in volume of the sample.	|	d	|	... |
|	cif20.histogram_completeness.py	|	A graphical representation of the completeness of the photometric sample in relation to 2MASS J.	|	Magnitudes from FUV to W4 |	Figure X	|	
|	cif20.histogram_distance.py	|	Produces a histogram of distances.	|	d	|	Figure	X |
|	cif20.histogram_logg.py	|	Produces a histogram of log g.	|	logg	|	Figure	X |
|	cif20.histogram_luminosities .py	|	Produces a histogram of luminosities.	|	Lbol	|	Figure	X |
|	cif20.histogram_magnitudes.py	|	Produces a histogram of magnitudes for each passband.	|	Magnitudes from FUV to W4	|	Figure	X |
|	cif20.histogram_RUWE.py	|	Produces a histogram of RUWE.	|	RUWE<sup id="a2">[2](#f2)</sup>	|	Figure	X |
|	cif20.histogram_SpT.py	|	Histogram of spectral types.	|	 SpTnum<sup id="a3">[3](#f3)</sup>	|	Figure	X |
|	cif20.histogram_Teff.py	|	Produces a histogram of effective temperatures.	|	Teff	|	Figure	X |
|	cif20.plot_appendix_colour_colour.py	|	Plots color vs. colour for Appendix.	|	Pairs of selected magnitudes.	|	Figure A.1. |
|	cif20.plot_appendix_colour_SpT.py	|	Plots color vs. SpT for Appendix.	|	Pairs of selected magnitudes.	|	Figure A.2.	|
|	cif20.plot_binaries_deltaG.py	|	Plots ∆G vs. ρ for binary stars resolved by *Gaia* DR2.	|	G mag, RA, DE	|	Figure	X |
|	cif20.plot_binaries_distances.py	|	Plots a comparison of parallaxes for primary and secondary components of binary systems (physical or not) resolved by *Gaia* DR2.	|	Plx	|	Figure	X |
|	cif20.plot_box_meta.py	|	Plots box & whiskers diagram comparing VOSA's BT-Settl against metallicity values from the literature.	|	[Fe/H] |	Figure	X |	
|	cif20.plot_literature_colour_colour_1.py	|	Plots J-H vs. g-i diagram plus literature values (Covey+07, Davenport+14).	|	JHgi mag	|	Figure X	|
|	cif20.plot_literature_colour_colour_2.py	|	Plots NUV-Ks vs. V-J diagram and the ‘basal NUV-Ks’ calculations for young stars from Ansdell+15	|	NUVVJKs magnitudes.	|	... |
|	cif20.plot_literature_colour_colour_3.py	|	Plots g-r vs. r-i diagram plus literature values (Bochanski+08, Davenport+14).	|	gri mag	|	Figure X	|
|	cif20.plot_literature_colour_SpT.py	|	Plots r-i vs. SpT diagram plus literature values (Hawley+02, Covey+07, Bochanski+08).	|	ri mag	|	Figure X	|	
|	cif20.plot_literature_LR_Teff.py	|	Plots Luminosity or Radius vs. Teff diagram plus Rabus+18 alleged discontinuity, alsong with Baraffe’s isochrones BCAH98 and DUSTY00.	|	Lbol, Teff. Baraffe’s isochrones.	|	-	|	
|	cif20.plot_literature_LTMR.py	|		Plots comparing Luminosity, Effective Temperature, Mass or Radius values with literature (many references, see main text in the article).	|	Lbol, Teff.	| Figure X |
|	cif20.plot_literature_LT.py	|	Plots Lbol vs. Teff diagram plus literature values (Pecaut+13, Newton+15, Faherty+16).	|	Lbol, Teff |	Figure X	|	
|	cif20.plot_literature_Mabs_colour.py	|	Plots MJ vs. J-Ks diagram plus literature values (Lepine+11,13, Knapp+04).	|	JKs mag, d	| Figure X |
|	cif20.plot_literature_Mabs_SpT.py	|	Plots MJ vs. SpT diagram plus literature values (Hawley+02, Kiman+19).	|	J mag, d, SpTnum	|	 Figure X	|	
|	cif20.plot_literature_Mabs_Teff.py	|	Plots MJ vs. Teff diagram plus literature values (Dahn+02, Lepine+14, Gaidos+14).	|	J mag, d, Teff	|	 Figure X	|
|	cif20.plot_literature_Mass.py	|	Plots a comparison of masses from this work vs. values derived from MKs-Mass literature models (Delfosse+00, Benedict+16, Mann+19).	|	Ks mag, d, Lbol, Teff	| Figure X |
|	cif20.plot_literature_MR_SpT.py	|	Plots Mass or Radius vs. SpT diagram plus literature values (Pecaut+13, Mann+15).	|	Lbol, Teff |	Figure X |
|	cif20.plot_literature_Teff.py	|	Plots Teff vs. Teff diagram plus literature values (Passegger+19).	|	Teff	|	Pas19	|	Figure X |
|	cif20.plot_literature_Teff_colour.py	|	Plots Teff vs. V-J diagram plus literature values (Casagrande+08, Pecaut+13).	|	VJ mag, Teff	|	Figure X	|	
|	cif20.plot_literature_Lbol_meta.py	|	Plots comparison of luminosities from BT-Settl and from BT-Settl CIFIST.	|	Lbol, [Fe/H]	|	Figure X	|
|	cif20.plot_literature_Meta_Meta.py	|	Plots comparison of metallicities from BT-Settl and from the literature (many references).	|	[Fe/H]	|	Figure	X |
|	cif20.plot_literature_Teff_meta.py	|	Plots comparison of effective temperatures from BT-Settl and from BT-Settl CIFIST.	|	Effective temperatures. Metallicities.	|	Figure	|	Plots comparison of effective temperatures from BT-Settl and BT-Settl CIFIST [Fe/H].
|	cif20.plot_model_BC_colour_meta.py	|	Plots BCG vs. G-J colour-coded by metallicity and fits a polynomial.	|	GJ mag, d, Lbol	|	Figure X, Table 5	|	
|	cif20.plot_model_BC_colour.py	|	Plots BCG, BCr, BCJ, BCW3 vs. G-J,  r-J, G-J, G-J, and fits two polynomial (2 ranges).	|	GrJW3 mag, d, Lbol	|	Figure X, Table 5,  Minitool: Illuminator BC™	|
|	cif20.plot_model_Mabs_colour_meta.py	|	Plots MG vs. G-J colour-coded by metallicity and fits a polynomial.	|	GJ mag, d	|	Figure X, Table 5. |
|	cif20.plot_model_Mabs_Lbol_meta.py	|	Plots Luminosity vs. MG or MJ colour-coded by metallicity and fits a polynomial.	|	J mag, d, Lbol |	Figure X, Table 5. |
|	cif20.plot_model_Mabs_colour.py	|	Plots MG vs. G-J or Mr vs. r-J and fits a polynomial.	|	GJr mag, d |	Figure X, Table 5, Minitool: Parsecator™	|
|	cif20.plot_model_Mabs_Lbol.py	|	Plots Luminosity vs. MG or MJ and fits two polynomials (2 ranges).|	GJ mag, d, Lbol	|	Figure X, Table 5, Minitool: Parsecator™ |
|	cif20.plot_standalone_cbar.py	|	Produces a standalone colorbar with the selected colormap.	|	A <a href="https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html" target="_blank">colormap</a>  |	Figures X and Y	|	
|	cif20.plot_BC_colour.py	|	Plots BC vs. G-J for every passband (FUV...W4) in a single canvas. | FUV to W4 mag, d, Lbol |	Figure X	|
|	cif20.plot_colour_Teff.py	|	Plots G-RP vs. Teff diagram with average values.	|	GRP mag, Teff	|	Figure X	|
|	cif20.plot_colour_SpT_meta.py	|	Plots G-J vs. SpT, colour-coded by metallicity.	|	GJ mag, [Fe/H]	|	Figure	X |	
|	cif20.plot_colour_colour.py	|	Plots BP-RP vs. G-J or J-H vs. H-Ks.	|	2MASS or *Gaia* DR2 mag, SpTnum	|	Figure	X |	
|	cif20.plot_colour_SpT.py	|	Plots G-J vs. SpT with mean values for each SpT.	|	GJ mag, Teff	|	Figure X	|
|	cif20.plot_Mabs_colour.py	|	Plots MG vs. G-J.	|	GJ mag, d, SpTnum	|	Figure	X |	
|	cif20.plot_RA_DE.py	|	Plots equatorial coordinates (Right Ascension and Declination) on a cartesian plane.	|	RA, DE	|	Figure	X |
|	cif20.plot_SED.py	|	Plots observed and theoretical SED on top of a PHOENIX synthetic spectrum.	|	Spectrum, transmission curves for the filters, sample stellar SED	|	Figure X	|	
|	cif20.plot_skymap.py	|	Plots with the equatorial or galactic positions of the stars in the J2000 equinox.	|	RA, DE	|	Figure X	|	
|	cif20.plot_filters.py	|	Plots a sample SED on top of the transmission curves of NUV-W4 passbands.	|	Transmission curves for the filters, sample stellar SED	|	Figure	X |	
|	cif20.plot_Lbol_Teff.py	|	Plots luminosity vs. Teff. with young stars shown in different colour. |	Lbol, Teff	|	Figure	X |	
|	cif20.plot_colour_excess.py	|	Plots *Gaia* DR2 excess factor vs. BP - RP with the colour excess limits from Evans+18 plotted.	|	BP RP magnitudes, *Gaia* DR2 `phot_bp_rp_excess_factor`	|	- |
|	cif20.table_generator.py	|	Produces Tables 6, 7 and 8. Plots mean values of Teff, Lbol, Mass, Radius.	|	All parameters in Tabes 6-8.	|	Tables 6, 7 and 8 in TeX format, ready to paste. Figures X and Y |
|	cif20.VOSA_input_generator.py	|	Generates an ASCII file from a set of photometric data, compatible with the VOSA input format.	|	Mags, d |	A VOSA compatible .txt with all stars. |

1. <small id="f1"> Lbol = bolometric luminosity (Lsol); Teff = effective temperature (K); M = mass (Msol); R = radius (Rsol); logg = surface gravity (dex); 	[Fe/H] = metallicity (dex); d = distance (pc); Plx = parallax (mas); RA, DE = equatorial coordinates in the J2000 equinox; 'mag' = magnitudes. Uncertainties are almost always used as an input, but they are omitted here for simplicity. Also, spectral types are almost always needed to colour-code the scatter plots. </small> [↩](#a1) 
2. <small id="f2"> Now available in the <a href="https://gea.esac.esa.int/archive/" target="_blank">*Gaia* Archive</a>. </small> [↩](#a2)
3. <small id="f3"> Spectral types in numerical form: K5V = -2, K7V = -1, M0.0 = 0.0, and so on. </small> [↩](#a3)

### Mother: the master table

`Mother.v01.csv` is the full version of the table associated to Cifuentes et al. 2020 and stored in [link to CDS]. It contains 2486 rows and 183 columns.

- Versions
  - First: v01 (July 2020)
  - Current: v01 (July 2020)

(Row-by-row description of the table available here).

---

## Support

Reach out to me at <a href="mailto:ccifuentes@cab.inta-csic.es">`ccifuentes@cab.inta-csic.es`</a>.
<Website at <a href="http://aplaceformyhead.es" target="_blank">`aplaceformyhead.es`</a>>

---

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**

---

## Suggested Resources

- <a href="https://www.python.org/dev/peps/pep-0008/" target="_blank">Style Guide for Python Code (PEP 8)</a>
