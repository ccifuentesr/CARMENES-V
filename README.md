# CARMENES input catalogue of M dwarfs V. Luminosities, colours, and spectral energy distributions

The repository contains:

- Directory ./: Includes all the code files named as `cif03.xxx_***.py`, where `xxx` defines the usage for the particular code ('plot', 'utilities', 'calculator', 'model'). 
- Directory ./Data: Includes all necessary auxiliary files.
- The **main table** (`cif03.full_table.csv`) [2634 rows, 132 columns] contains astrometric and photometric data, fundamental parameters, and multiplicity information of all the stars in the sample and their physically bound companions.

Reach out to me at <a href="mailto:ccifuentes@cab.inta-csic.es">`ccifuentes@cab.inta-csic.es`</a>.

The *CARMENES input catalogue of M dwarfs* series:

- <a href="https://ui.adsabs.harvard.edu/abs/2024arXiv241212264C/abstract" target="_blank">IX. Multiplicity from close spectroscopic binaries to ultrawide systems</a>  Cifuentes et al 2024.
- <a href="https://ui.adsabs.harvard.edu/abs/2024A%26A...692A.206C/abstract" target="_blank">VIII. Kinematics in the solar neighbourhood</a>  Cortés-Contreras et al. 2024.
- <a href="https://ui.adsabs.harvard.edu/abs/2024A%26A...684A...9S/abstract" target="_blank">VII. New rotation periods for the survey stars and their correlations with stellar activity</a>  Shan et al. 2024
- <a href="https://ui.adsabs.harvard.edu/abs/2021A%26A...652A.116P/abstract" target="_blank">VI. A time-resolved Ca II H&K catalog from archival data</a>  Perdelwitz et al 2021.
- <a href="https://ui.adsabs.harvard.edu/abs/2020A%26A...642A.115C/abstract" target="_blank">**V. Luminosities, colours, and spectral energy distributions**</a>  **Cifuentes et al 2020 (this work).**
- <a href="https://ui.adsabs.harvard.edu/abs/2019A%26A...621A.126D/abstract" target="_blank">IV. New rotation periods from photometric time series </a> Díez Alonso et al. 2019.
- <a href="https://ui.adsabs.harvard.edu/abs/2018A%26A...614A..76J/abstract" target="_blank">III. Rotation and activity from high-resolution spectroscopic observations </a> Jeffers et al. 2018.
- <a href="https://ui.adsabs.harvard.edu/abs/2017A%26A...597A..47C/abstract" target="_blank">II. High-resolution imaging with FastCam</a> Cortés-Contreras et al. 2017.
- <a href="https://ui.adsabs.harvard.edu/abs/2015A%26A...577A.128A/abstract" target="_blank">I. Low-resolution spectroscopy with CAFOS</a> Alonso-Floriano et al. 2015.

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Publication](https://img.shields.io/badge/Published%3F-yes-brightgreen.svg)](https://www.aanda.org/articles/aa/abs/2020/10/aa38295-20/aa38295-20.html)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://GitHub.com/ccifuentesr)

See the <a href="https://carmenes.caha.es" target="_blank">CARMENES Website</a>.

## Mother: the master table

`cif20.Mother.v01.csv` (2483 rows and 175 columns) is the full version of Table A.3 (summary table) in <a href="https://arxiv.org/abs/2007.15077" target="_blank">Cifuentes et al. 2020</a>, available at the CDS via anonymous ftp to cdsarc.u-strasbg.fr(130.79.128.5).

## Directories

- Directory ./: Includes all the code files classified as detailed below and the master table (`cif20.Mother.v01.csv`). [57 files, 3.1 MB]
- Directory ./Baraffe: (see <a href="https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/" target="_blank">Baraffe et al. 2015</a>). [7 items, 33 KB]
- Directory ./Filters: Includes the transmission curves for 20 filters (FUV to W4). [24 items 138 KB]
- Directory ./Literature: Includes the individual tables with data from the literature. [33 items, 3.9 MB]
- Directory ./Spectra: Includes PHOENIX synthetic spectra. [1 item, 22 MB]
- Directory ./Stars: Includes data for individual stars to produce custom SEDs. [4 items, 20 KB]

The total size is 29 MB. If the `./Spectra` directory is not needed (only used in `cif20.plot_SED.py`, see below for more details), the total size shrinks to less than 7 MB.

## Files

All filenames are formatted as `cif20.xxx_yyy_zzz.py`, where `xxx` defines the kind of output that it produces, `yyy` gives additional information about the output, and `zzz` enumerates the main variables involved.

| File | Description | Input<sup id="a1">[1](#f1)</sup>| Output<sup id="a2">[2](#f2)</sup> | 
| --- | --- | --- | --- | 
| cif20.calculator_averagecolors.py | Average colours and standard deviations for each spectral type from K5V to L8. | Magnitudes from FUV to W4 | Table A.2 |
  |	cif20.calculator_completeness.py |	Estimator of the completeness in volume of the sample.	|	d	|	Main text. |
| cif20.calculator_Karmn.py | A generator of Carmencita identification name (Karmn JHHMMm+DDdAAA). | Simbad valid name | Karmn ID |
|	cif20.calculator_MR.py	|	Calculator of mass and radius in solar units.	|	Lbol, Teff	|	- |	
|	cif20.histogram_completeness.py	|	A graphical representation of the completeness of the photometric sample in relation to 2MASS J.	|	Magnitudes from FUV to W4 |	Figure 4 |	
|	cif20.histogram_distance.py	|	Produces a histogram of distances.	|	d	|	Figure	7 |
|	cif20.histogram_logg.py	|	Produces a histogram of log g.	|	logg	|	Figure	11 |
|	cif20.histogram_luminosities .py	|	Produces a histogram of luminosities.	|	Lbol	|	Figure	11 |
|	cif20.histogram_magnitudes.py	|	Produces a histogram of magnitudes for each passband.	|	Magnitudes from FUV to W4	|	Figure 5 |
|	cif20.histogram_RUWE.py	|	Produces a histogram of RUWE.	|	RUWE<sup id="a3">[3](#f3)</sup>	|	Figure	7 |
|	cif20.histogram_SpT.py	|	Histogram of spectral types.	|	 SpTnum<sup id="a4">[4](#f4)</sup>	|	Figure	1 |
|	cif20.histogram_Teff.py	|	Produces a histogram of effective temperatures.	|	Teff	|	Figure	11 |
|	cif20.plot_appendix_colour_colour.py	|	Plots color vs. colour for Appendix.	|	Pairs of selected magnitudes.	|	Figure A.2. |
|	cif20.plot_appendix_colour_SpT.py	|	Plots color vs. SpT for Appendix.	|	Pairs of selected magnitudes.	|	Figure A.1.	|
|	cif20.plot_BC_colour.py	|	Plots BC vs. G-J for every passband (FUV...W4) in a single canvas. | FUV to W4 mag, d, Lbol |	Figure 17	|
|	cif20.plot_binaries_deltaG.py	|	Plots ∆G vs. ρ for binary stars resolved by *Gaia* DR2.	|	G mag, RA, DE	|	Figure	9 |
|	cif20.plot_binaries_distances.py	|	Plots a comparison of parallaxes for primary and secondary components of binary systems (physical or not) resolved by *Gaia* DR2.	|	Plx	|	Figure	8 |
|	cif20.plot_box_meta.py	|	Plots box & whiskers diagram comparing VOSA's BT-Settl against metallicity values from the literature.	|	[Fe/H] |	Figure	22 |
|	cif20.plot_colour_colour.py	|	Plots BP-RP vs. G-J or J-H vs. H-Ks.	|	2MASS or *Gaia* DR2 mag, SpTnum	|	Figure	14 |	
|	cif20.plot_colour_excess.py	|	Plots *Gaia* DR2 excess factor vs. BP - RP with the colour excess limits from Evans+18 plotted.	|	BP RP magnitudes, *Gaia* DR2 `phot_bp_rp_excess_factor`	|	- |
|	cif20.plot_colour_SpT.py	|	Plots G-J vs. SpT with mean values for each SpT.	|	GJ mag, Teff	|	Figure 13	|
|	cif20.plot_colour_SpT_meta.py	|	Plots G-J vs. SpT, colour-coded by metallicity.	|	GJ mag, [Fe/H]	|	Figure	21 |
|	cif20.plot_colour_Teff.py	|	Plots G-RP vs. Teff diagram with average values.	|	GRP mag, Teff	|	-	|
|	cif20.plot_filters.py	|	Plots a sample SED on top of the transmission curves of NUV-W4 passbands.	|	Transmission curves for the filters, sample stellar SED	|	Figure	3 |
|	cif20.plot_Lbol_Teff.py	|	Plots luminosity vs. Teff. with young stars shown in different colour. |	Lbol, Teff	|	Figure	12 |
|	cif20.plot_Mabs_colour.py	|	Plots MG vs. G-J.	|	GJ mag, d, SpTnum	|	Figure	12 |	
|	cif20.plot_RA_DE.py	|	Plots equatorial coordinates (Right Ascension and Declination) on a cartesian plane.	|	RA, DE	|	- |
|	cif20.plot_SED.py	|	Plots observed and theoretical SED on top of a PHOENIX synthetic spectrum.	|	Spectrum, transmission curves for the filters, sample stellar SED	|	Figure 10	|	
|	cif20.plot_skymap.py	|	Plots with the equatorial or galactic positions of the stars in the J2000 equinox.	|	RA, DE	|	Figure 2	|	
|	cif20.plot_standalone_cbar.py	|	Produces a standalone colorbar with the selected colormap.	|	A <a href="https://matplotlib.org/3.1.1/gallery/color/colormap_reference.html" target="_blank">colormap</a>  |	Figures 2, 19, 20, 21, 23, 24, 26, A.1 and A.2	|	
|	cif20.plot_literature_colour_colour_1.py	|	Plots J-H vs. g-i diagram plus literature values (Covey+07, Davenport+14).	|	JHgi mag	|	-	|
|	cif20.plot_literature_colour_colour_2.py	|	Plots NUV-Ks vs. V-J diagram and the ‘basal NUV-Ks’ calculations for young stars from Ansdell+15	|	NUVVJKs magnitudes.	|	- |
|	cif20.plot_literature_colour_colour_3.py	|	Plots g-r vs. r-i diagram plus literature values (Bochanski+08, Davenport+14).	|	gri mag	|	Figure 26	|
|	cif20.plot_literature_colour_SpT.py	|	Plots r-i vs. SpT diagram plus literature values (Hawley+02, Covey+07, Bochanski+08).	|	ri mag	|	Figure 26	|	
|	cif20.plot_literature_Lbol_meta.py	|	Plots comparison of luminosities from BT-Settl and from BT-Settl CIFIST.	|	Lbol, [Fe/H]	|	Figure 23	|
|	cif20.plot_literature_LTMR.py	|		Plots comparing Luminosity, Effective Temperature, Mass or Radius values with literature (many references, see main text in the article).	|	Lbol, Teff.	| Figures 19, 20, 24 and 25 |
|	cif20.plot_literature_LT.py	|	Plots Lbol vs. Teff diagram plus literature values (Pecaut+13, Newton+15, Faherty+16).	|	Lbol, Teff |	Figure 20	|	
|	cif20.plot_literature_LR_Teff.py	|	Plots Luminosity or Radius vs. Teff diagram plus Rabus+18 alleged discontinuity, alsong with Baraffe’s isochrones BCAH98 and DUSTY00.	|	Lbol, Teff. Baraffe’s isochrones.	|	Figure 24	|	
|	cif20.plot_literature_Mabs_colour.py	|	Plots MJ vs. J-Ks diagram plus literature values (Lepine+11,13, Knapp+04).	|	JKs mag, d	| Figure 26 |
|	cif20.plot_literature_Mabs_SpT.py	|	Plots MJ vs. SpT diagram plus literature values (Hawley+02, Kiman+19).	|	J mag, d, SpTnum	|	 -	|	
|	cif20.plot_literature_Mabs_Teff.py	|	Plots MJ vs. Teff diagram plus literature values (Dahn+02, Lepine+14, Gaidos+14).	|	J mag, d, Teff	|	 Figure 20	|
|	cif20.plot_literature_Mass.py	|	Plots a comparison of masses from this work vs. values derived from MKs-Mass literature models (Delfosse+00, Benedict+16, Mann+19).	|	Ks mag, d, Lbol, Teff	| Figure 25 |
|	cif20.plot_literature_Meta_Meta.py	|	Plots comparison of metallicities from BT-Settl and from the literature (many references).	|	[Fe/H]	|	Figure	22 |
|	cif20.plot_literature_MR_SpT.py	|	Plots Mass or Radius vs. SpT diagram plus literature values (Pecaut+13, Mann+15).	|	Lbol, Teff |	Figure 24 |
|	cif20.plot_literature_Teff.py	|	Plots Teff vs. Teff diagram plus literature values (Passegger+19).	|	Teff	|	Pas19	|	- |
|	cif20.plot_literature_Teff_colour.py	|	Plots Teff vs. V-J diagram plus literature values (Casagrande+08, Pecaut+13).	|	VJ mag, Teff	|	Figure 26	|	
|	cif20.plot_literature_Teff_meta.py	|	Plots comparison of effective temperatures from BT-Settl and from BT-Settl CIFIST.	|	Effective temperatures. Metallicities.	|	- |	
|	cif20.plot_model_BC_colour.py	|	Plots BCG, BCr, BCJ, BCW3 vs. G-J,  r-J, G-J, G-J, and fits two polynomial (2 ranges).	|	GrJW3 mag, d, Lbol	|	Figure 18, Table 5,  Minitool: Illuminator BC<sup id="a5">[5](#f5)</sup>	|
|	cif20.plot_model_BC_colour_meta.py	|	Plots BCG vs. G-J colour-coded by metallicity and fits a polynomial.	|	GJ mag, d, Lbol	|	Table 5	|
|	cif20.plot_model_Mabs_colour.py	|	Plots MG vs. G-J or Mr vs. r-J and fits a polynomial.	|	GJr mag, d |	Figure 15, Table 5, Minitool: Parsecator<sup id="a5">[5](#f5)</sup>	|
|	cif20.plot_model_Mabs_colour_meta.py	|	Plots MG vs. G-J colour-coded by metallicity and fits a polynomial.	|	GJ mag, d	|	Figure 21, Table 5. |
|	cif20.plot_model_Mabs_Lbol.py	|	Plots Luminosity vs. MG or MJ and fits two polynomials (2 ranges).|	GJ mag, d, Lbol	|	Figure 16, Table 5, Minitool: Illuminator<sup id="a5">[5](#f5)</sup> |
|	cif20.plot_model_Mabs_Lbol_meta.py	|	Plots Luminosity vs. MG or MJ colour-coded by metallicity and fits a polynomial.	|	J mag, d, Lbol |	Figure 21, Table 5. |
|	cif20.table_generator.py	|	Produces Tables 6, 7 and 8. Plots mean values of Teff, Lbol, Mass, Radius.	|	All parameters in Tabes 6-8	|	Tables 6, 7 and 8 in TeX format, ready to paste. Figures 19 and 20 |
|	cif20.VOSA_input_generator.py	|	Generates an ASCII file from a set of photometric data, compatible with the VOSA input format.	|	Mags, d |	A VOSA compatible .txt with all stars. |

1. <small id="f1"> Lbol = bolometric luminosity (Lsol); Teff = effective temperature (K); M = mass (Msol); R = radius (Rsol); logg = surface gravity (dex); 	[Fe/H] = metallicity (dex); d = distance (pc); Plx = parallax (mas); RA, DE = equatorial coordinates in the J2000 equinox; 'mag' = magnitudes. Uncertainties are almost always used as an input, but they are omitted here for simplicity. Also, spectral types are almost always needed to colour-code the scatter plots. </small> [↩](#a1) 
2. <small id="f2"> Specific element in Cifuentes et al. 2020 that each script produces. </small> [↩](#a2)
3. <small id="f3"> Now available in the <a href="https://gea.esac.esa.int/archive/" target="_blank">*Gaia* Archive</a>.</small> [↩](#a3)
4. <small id="f4"> Spectral types in numerical form: K5V = -2, K7V = -1, M0.0 = 0.0, and so on. </small> [↩](#a4)
5. <small id="f5"> The minitools *Parsecator* and *Illuminator* estimate distances and luminosities, respectively, based on the models produced in their respective files. They are provided as individual code cells meant to be run as standalone pieces of code, always after the models have been produced. </small> [↩](#a5)
