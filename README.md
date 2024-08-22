# CARMENES input catalogue of M dwarfs V. Luminosities, colours, and spectral energy distributions
  
> This repository contains the pieces of code necessary to produce all figures, tables and models in <a href="https://arxiv.org/abs/2007.15077" target="_blank">Cifuentes et al. 2020</a>.

The *CARMENES input catalogue of M dwarfs* series contains:

- <a href="https://www.aanda.org/articles/aa/abs/2020/10/aa38295-20/aa38295-20.html" target="_blank">**CARMENES input catalogue of M dwarfs V. Luminosities, colours, and spectral energy distributions** </a>  **Cifuentes et al 2020 (This work)**.
- <a href="https://ui.adsabs.harvard.edu/abs/2019A%26A...621A.126D/abstract" target="_blank">*CARMENES input catalogue of M dwarfs. IV. New rotation periods from photometric time series* </a> Díez Alonso et al. 2019.
- <a href="https://ui.adsabs.harvard.edu/abs/2018A%26A...614A..76J/abstract" target="_blank">*CARMENES input catalogue of M dwarfs. III. Rotation and activity from high-resolution spectroscopic observations* </a> Jeffers et al. 2018.
- <a href="https://ui.adsabs.harvard.edu/abs/2017A%26A...597A..47C/abstract" target="_blank">*CARMENES input catalogue of M dwarfs. II. High-resolution imaging with FastCam*</a> Cortés-Contreras et al. 2017.
- <a href="https://ui.adsabs.harvard.edu/abs/2015A%26A...577A.128A/abstract" target="_blank">*CARMENES input catalogue of M dwarfs. I. Low-resolution spectroscopy with CAFOS*</a> Alonso-Floriano et al. 2015.

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Publication](https://img.shields.io/badge/Published%3F-yes-brightgreen.svg)](https://www.aanda.org/articles/aa/abs/2020/10/aa38295-20/aa38295-20.html)
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
- They may be run as Python Notebooks. The symbol `# %%` starts a cell that can be run separately.
- Most of the `*.py` files require of additional data contained in the folders stored in the repository.
- Cloning or downloading the complete repository is strongly recommended (see below).
- The installation of some basic libraries is a prerequisite: `numpy`, `scipy`, `astropy`, `matplotlib` and `pyperclip`. Other modules are included in the Python distribution and do not need additional installation (e.g., `csv`).

### Clone

- Clone this repo to your local machine using `git clone https://github.com/ccifuentesr/CARMENES-V`, or
- Download this repo as a .zip and run the scripts in your local machine.

## Structure

### Directories

- Directory ./: Includes all the code files classified as detailed below and the master table (`cif20.Mother.v01.csv`). [57 files, 3.1 MB]
- Directory ./Baraffe: (see <a href="https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/" target="_blank">Baraffe et al. 2015</a>). [7 items, 33 KB]
- Directory ./Filters: Includes the transmission curves for 20 filters (FUV to W4). [24 items 138 KB]
- Directory ./Literature: Includes the individual tables with data from the literature. [33 items, 3.9 MB]
- Directory ./Spectra: Includes PHOENIX synthetic spectra. [1 item, 22 MB]
- Directory ./Stars: Includes data for individual stars to produce custom SEDs. [4 items, 20 KB]

The total size is 29 MB. If the `./Spectra` directory is not needed (only used in `cif20.plot_SED.py`, see below for more details), the total size shrinks to less than 7 MB.

### Files

All files are named `cif20.xxx_yyy_zzz.py`, where `xxx` defines the kind of output that it produces, `yyy` gives additional information about the output, and `zzz` enumerates the main variables involved. For example, the script `cif20.plot_literature_Mabs_SpT.py` produces an absolute magnitude vs. spectral type plot, and compares the values with those of the literature. The complete list of files and their description can be found below sorted by output, in alphabetical order.

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

**Notes:** 
- The output image format is `.eps` in all cases, which is the most suitable format for publication in journals. For a more manageable format, change `.eps` to `.png` in the code.

### Mother: the master table

`cif20.Mother.v01.csv` is the full version of Table A.3 (summary table) in <a href="https://arxiv.org/abs/2007.15077" target="_blank">Cifuentes et al. 2020</a>, available at the CDS via anonymous ftp to cdsarc.u-strasbg.fr(130.79.128.5) or via [*link to CDS in progress*]. The complete version available here contains 2483 rows and 175 columns. It is stored in the root directory and can be manipulated separately with tabular data management software such as <a href="http://www.star.bris.ac.uk/~mbt/topcat/" target="_blank">TOPCAT</a>.

**Version history**

|	Version	|	Date	| Comments |
|	---	|	---	| --- |
| 01 | July 2020 | Current version |

Row-by-row description of `cif20.Mother.v01.csv` available below.

|	ID	|	Name	|	Units	|	Description	|	Annotations	|
|	---	|	---	|	---	|	---	|	---	|
|	1	|	Karmn	|	-	|	Star identifier (JHHMMm+DDdAAA)	|	"For the M and L dwarfs the identifier "JHHMMm+DDdAAA" is used - plus "N" - “S” - “E" or "W" if two stars (in a close binary system) have the same  “HHMMm+DDd" string. For the K dwarfs the SUPERBLINK catalogue identifier (Lepine & Shara 2005; Lepine et al. 2013) is used. For the ultracool dwarfs we tabulate the *Gaia* UltraCool Dwarf Catalogue identifier (Smart et al. 2017; 2019).	|
|	2	|	Name	|	-	|	Discovery or most common name 	|		|
|	3	|	RA_J2000	|	hms	|	Right ascension (J2000.0 epoch) 	|	In equinox J2000	|
|	4	|	DE_J2000	|	dms	|	Declination (J2000.0 epoch) 	|	In equinox J2000	|
|	5	|	RA_J2015	|	hms	|	Right ascension (J2015.5 epoch)	|	In equinox J2000	|
|	6	|	DE_J2015	|	dms	|	Declination (J2015.5 epoch)	|	In equinox J2000	|
|	7	|	SpType	|	-	|	Spectral type 	|		|
|	8	|	SpTnum	|	-	|	Spectral type in numerical format 	|	SpTnum = -2.0 for K5V ; -1.0 for K7V ; 0.0 for M0.0V ; 0.5 for M0.5V ; ... ; 10.0 for L0.0 ; 10.5 for L0.5 ; etc.	|
|	9	|	Ref_SpT	|	-	|	Reference for the spectral type 	|		|
|	10	|	Plx	|	mas	|	Parallax 	|		|
|	11	|	ePlx	|	mas	|	Parallax error 	|		|
|	12	|	Ref_Plx	|	-	|	Reference for the parallax 	|		|
|	13	|	d_pc	|	pc	|	Distance 	|		|
|	14	|	ed_pc	|	pc	|	Distance error 	|		|
|	15	|	Ref_d	|	-	|	Reference for the distance 	|		|
|	16	|	Lbol	|	solLum	|	Bolometric luminosity from VOSA 	|		|
|	17	|	Lberr	|	solLum	|	Bolometric luminosity error from VOSA 	|		|
|	18	|	Teff	|	K	|	Effective temperature from VOSA	|	BT-Settl CIFIST grid of synthetic spectra uncertainties are 50 K for Teff (25 K for Teff <= 2400 K).	|
|	19	|	logg	|	-	|	Surface gravity from VOSA	|	BT-Settl CIFIST grid of synthetic spectra uncertainties are 0.5 dex for logg.	|
|	20	|	Radius	|	solRad	|	Stellar radius	|	Radii calculated from Lbol and Teff using the Stefan-Boltzmann law under the Black Body approximation.	|
|	21	|	eRadius	|	solRad	|	Stellar radius error 	|		|
|	22	|	Mass	|	solMass	|	Stellar mass	|	Masses are derived from the radii using the radius-mass relation for eclipsing binaries by Schweitzer et al. 2019.	|
|	23	|	eMass	|	solMass	|	Stellar mass error 	|		|
|	24	|	FUV_mag	|	mag	|	*GALEX* far-UV magnitude 	|		|
|	25	|	eFUV_mag	|	mag	|	*GALEX* far-UV magnitude 	|		|
|	26	|	FUV_mag	|	mag	|	*GALEX* near-UV magnitude 	|		|
|	27	|	eFUV_mag	|	mag	|	*GALEX* near-UV magnitude 	|		|
|	28	|	Qf_FUV	|	-	|	*GALEX* far-UV magnitude quality flag	|		|
|	29	|	Qf_NUV	|	-	|	*GALEX* near-UV magnitude quality flag	|		|
|	30	|	Ref_FNUV	|	-	|	*GALEX* far-UV and near-UV magnitudes reference 	|		|
|	31	|	u_mag	|	mag	|	*u* magnitude 	|		|
|	32	|	eu_mag	|	mag	|	*u* magnitude error 	|		|
|	33	|	Qf_u	|	-	|	*u* magnitude quality flag 	|		|
|	34	|	Ref_u	|	-	|	*u* magnitude reference 	|		|
|	35	|	BT_mag	|	mag	|	Tycho-2 *BT* magnitude	|		|
|	36	|	eBT_mag	|	mag	|	Tycho-2 *BT* magnitude error	|		|
|	37	|	Qf_BT	|	-	|	Tycho-2 *BT* magnitude quality flag	|		|
|	38	|	Ref_BT	|	-	|	Tycho-2 *BT* magnitude reference	|		|
|	39	|	B_mag	|	mag	|	Johnson *B* magnitude 	|		|
|	40	|	eB_mag	|	mag	|	Johnson *B* magnitude error 	|		|
|	41	|	Qf_B	|	-	|	Johnson *B* magnitude quality flag 	|		|
|	42	|	Ref_B	|	-	|	Johnson *B* magnitude reference 	|		|
|	43	|	g_mag	|	mag	|	*g* magnitude 	|		|
|	44	|	eg_mag	|	mag	|	*g* magnitude error 	|		|
|	45	|	Qf_g	|	-	|	*g* magnitude quality flag 	|		|
|	46	|	Ref_g	|	-	|	*g* magnitude reference 	|		|
|	47	|	BP_mag	|	mag	|	*Gaia* DR2 *BP* magnitude	|		|
|	48	|	eBP_mag	|	mag	|	*Gaia* DR2 *BP* magnitude error	|		|
|	49	|	Qf_BP	|	-	|	*Gaia* DR2 *BP* magnitude quality flag	|		|
|	50	|	Ref_BP	|	-	|	*Gaia* DR2 *BP* magnitude reference	|		|
|	51	|	VT_mag	|	mag	|	Tycho-2 *VT* magnitude	|		|
|	52	|	eVT_mag	|	mag	|	Tycho-2 *VT* magnitude error	|		|
|	53	|	Qf_VT	|	-	|	Tycho-2 *VT* magnitude quality flag	|		|
|	54	|	Ref_VT	|	-	|	Tycho-2 *VT* magnitude reference	|		|
|	55	|	V_mag	|	mag	|	Johnson *V* magnitude 	|		|
|	56	|	eV_mag	|	mag	|	Johnson *V* magnitude error 	|		|
|	57	|	Qf_V	|	-	|	Johnson *V* magnitude quality flag 	|		|
|	58	|	Ref_V	|	-	|	Johnson *V* magnitude reference 	|		|
|	59	|	r_mag	|	mag	|	*r* magnitude 	|		|
|	60	|	er_mag	|	mag	|	*r* magnitude error 	|		|
|	61	|	Qf_r	|	-	|	*r* magnitude quality flag 	|		|
|	62	|	Ref_r	|	-	|	*r* magnitude reference 	|		|
|	63	|	GG_mag	|	mag	|	*Gaia* DR2 *G* magnitude	|		|
|	64	|	eGG_mag	|	mag	|	*Gaia* DR2 *G* magnitude error	|		|
|	65	|	Qf_GG	|	-	|	*Gaia* DR2 *G* magnitude quality flag	|		|
|	66	|	Ref_GG	|	-	|	*Gaia* DR2 *G* magnitude reference	|		|
|	67	|	i_mag	|	mag	|	*i* magnitude 	|		|
|	68	|	ei_mag	|	mag	|	*i* magnitude error 	|		|
|	69	|	Qf_i	|	-	|	*i* magnitude quality flag 	|		|
|	70	|	Ref_i	|	-	|	*i* magnitude reference 	|		|
|	71	|	RP_mag	|	mag	|	*Gaia* DR2 *RP* magnitude	|		|
|	72	|	eRP_mag	|	mag	|	*Gaia* DR2 *RP* magnitude error	|		|
|	73	|	Qf_RP	|	-	|	*Gaia* DR2 *RP* magnitude quality flag	|		|
|	74	|	Ref_RP	|	-	|	*Gaia* DR2 *RP* magnitude reference	|		|
|	75	|	J_mag	|	mag	|	2MASS *J* magnitude	|		|
|	76	|	eJ_mag	|	mag	|	2MASS *J* magnitude error	|		|
|	77	|	H_mag	|	mag	|	2MASS *H* magnitude	|		|
|	78	|	eH_mag	|	mag	|	2MASS *H* magnitude error	|		|
|	79	|	J_mag	|	mag	|	2MASS *Ks* magnitude	|		|
|	80	|	eJ_mag	|	mag	|	2MASS *Ks* magnitude error	|		|
|	81	|	Qf_2M	|	-	|	2MASS *JHKs* three-character quality flag	|		|
|	82	|	Qf_J	|	-	|	2MASS *J* magnitude quality flag	|		|
|	83	|	Qf_H	|	-	|	2MASS *H* magnitude quality flag	|		|
|	84	|	Qf_Ks	|	-	|	2MASS *Ks* magnitude quality flag	|		|
|	85	|	Ref_JHK	|	-	|	2MASS *JHKs* magnitudes reference	|		|
|	86	|	W1_mag	|	mag	|	WISE *W1* magnitude	|		|
|	87	|	eW1_mag	|	mag	|	WISE *W1* magnitude error	|		|
|	88	|	W2_mag	|	mag	|	WISE *W2* magnitude	|		|
|	89	|	eW2_mag	|	mag	|	WISE *W2* magnitude error	|		|
|	90	|	W3_mag	|	mag	|	WISE *W3* magnitude	|		|
|	91	|	eW3_mag	|	mag	|	WISE *W3* magnitude error	|		|
|	92	|	W4_mag	|	mag	|	WISE *W4* magnitude	|		|
|	93	|	eW4_mag	|	mag	|	WISE *W4* magnitude error	|		|
|	94	|	Qf_Ws	|	-	|	WISE *W1W2W3W4* four-character quality flag	|		|
|	95	|	Qf_W1	|	-	|	WISE *W1* magnitude quality flag	|		|
|	96	|	Qf_W2	|	-	|	WISE *W2* magnitude quality flag	|		|
|	97	|	Qf_W3	|	-	|	WISE *W3* magnitude quality flag	|		|
|	98	|	Qf_W4	|	-	|	WISE *W4* magnitude quality flag	|		|
|	99	|	Ref_Ws	|	-	|	WISE *W1W2W3W4* magnitude  reference	|	'WISE': AllWISE (Cutri et al. 2014); 'WISE*': WISE (Cutri et al. 2012).	|
|	100	|	Teff_meta	|	K	|	Effective temperature from VOSA	|		|
|	101	|	logg_meta	|	dex	|	Surface gravity from VOSA (BT Settl)	|		|
|	102	|	Meta_meta	|	dex	|	Metallicity from VOSA (BT Settl)	|		|
|	103	|	Lbol_meta	|	solLum	|	Bolometric luminosity from VOSA (BT Settl)	|		|
|	104	|	Lberr_meta	|	solLum	|	Bolometric luminosity error from VOSA (BT Settl)	|		|
|	105	|	Lbol_lit	|	solLum	|	Bolometric luminosity from the literature	|		|
|	106	|	Lberr_lit	|	solLum	|	Bolometric luminosity error from the literature	|		|
|	107	|	M_lit	|	solMass	|	Mass from the literature	|		|
|	108	|	eM_lit	|	solMass	|	Mass error from the literature	|		|
|	109	|	R_lit	|	solRad	|	Radius from the literature	|		|
|	110	|	eR_lit	|	solRad	|	Radius error from the literature	|		|
|	111	|	Teff_lit	|	K	|	Effective temperature from the literature	|		|
|	112	|	eTeff_lit	|	K	|	Effective temperature error from the literature	|		|
|	113	|	Ref_LMRT	|	-	|	Literature values reference	|		|
|	114	|	Teff_Pas19	|	K	|	Effective temperature from Passegger+19	|		|
|	115	|	eTeff_Pas19	|	K	|	Effective temperature error from Passegger+19	|		|
|	116	|	FeH_lit	|	dex	|	Metallicity from the literature	|		|
|	117	|	eFeH_lit	|	dex	|	Metallicity error from the literature	|		|
|	118	|	FeH_lit_ref 	|	-	|	Metallicity reference	|		|
|	119	|	Binarity	|	-	|	Kind of multiplicity for this object	|	'New’ are candidates not tabulated by the WDS; ‘Known’ are binaries tabulated by the WDS; ‘Background’ are candidates for visual binaries without physical binding	|
|	120	|	rho	|	arcsec	|	Angular separation with the closest component in multiple systems	|		|
|	121	|	theta	|	deg	|	Position angle	|		|
|	122	|	RUWE	|	-	|	Re-normalised Unit Weight Error	|	See https://www.cosmos.esa.int/web/gaia/dr2-known-issues	|
|	123	|	DeltaG	|	mag	|	G magnitude difference with the closest component in multiple systems	|		|
|	124	|	gaia_id_1	|	-	|	*Gaia* DR2 identifier of single or primary star	|		|
|	125	|	ra_1	|	hms	|	Right ascension (J2015.5 epoch)	|		|
|	126	|	ra_error_1	|	hms	|	Right ascension error (J2015.5 epoch) 	|		|
|	127	|	dec_1	|	dms	|	Declination (J2015.5 epoch)	|		|
|	128	|	dec_error_1	|	dms	|	Declination error (J2015.5 epoch) 	|		|
|	129	|	parallax_1	|	mas	|	Parallax	|		|
|	130	|	parallax_error_1	|	mas	|	Parallax error 	|		|
|	131	|	pmra_1	|	mas yr-1	|	Proper motion in Right Ascension	|		|
|	132	|	pmra_error_1	|	mas yr-1	|	Proper motion error in Right Ascension	|		|
|	133	|	pmdec_1	|	mas yr-1	|	Proper motion in Declination	|		|
|	134	|	pmdec_error_1	|	mas yr-1	|	Proper motion error in Declination	|		|
|	135	|	pmtotal_1	|	mas yr-1	|	Total proper motion	|		|
|	136	|	pmtotal_error_1	|	mas yr-1	|	Total proper motion error	|		|
|	137	|	phot_bp_rp_excess_factor_1	|	dex	|	*BP*/*RP* excess factor	|		|
|	138	|	radial_velocity_1	|	km s-1	|	Radial velocity	|		|
|	139	|	radial_velocity_error_1	|	km s-1	|	Radial velocity error	|		|
|	140	|	l_1	|	deg	|	Galactic longitude for B	|		|
|	141	|	b_1	|	deg	|	Galactic latitude for B	|		|
|	142	|	teff_val_1	|	K	|	Stellar effective temperature for B	|		|
|	143	|	lum_val_1	|	solLum	|	Stellar luminosity for B	|		|
|	144	|	gaia_id_2	|	-	|	*Gaia* DR2 identifier for B	|	The second Gaia DR2 identifier is provided for physically bound binary systems visual binaries separated less than 5 arcsec	|
|	145	|	ra_2	|	hms	|	Right ascension (J2015.5 epoch) for B	|	In equinox J2000	|
|	146	|	ra_error_2	|	hms	|	Right ascension error (J2015.5 epoch) for B	|	In equinox J2000	|
|	147	|	dec_2	|	dms	|	Declination (J2015.5 epoch) for B	|	In equinox J2000	|
|	148	|	dec_error_2	|	dms	|	Declination error (J2015.5 epoch) for B	|	In equinox J2000	|
|	149	|	parallax_2	|	mas	|	Parallax for B	|		|
|	150	|	parallax_error_2	|	mas	|	Parallax error for B	|		|
|	151	|	pmra_2	|	mas yr-1	|	Proper motion in Right Ascension for B	|		|
|	152	|	pmra_error_2	|	mas yr-1	|	Proper motion error in Right Ascension for B	|		|
|	153	|	pmdec_2	|	mas yr-1	|	Proper motion in Declination for B	|		|
|	154	|	pmdec_error_2	|	mas yr-1	|	Proper motion error in Declination for B	|		|
|	155	|	pmtotal_2	|	mas yr-1	|	Total proper motion for B	|		|
|	156	|	pmtotal_error_2	|	mas yr-1	|	Total proper motion error for B	|		|
|	157	|	phot_g_mean_mag_2	|	mag	|	G-band mean magnitude for B	|		|
|	158	|	phot_g_mean_mag_error_2	|	mag	|	G-band mean magnitude error for B	|		|
|	159	|	phot_bp_mean_mag_2	|	mag	|	BP-band mean magnitude for B	|		|
|	160	|	phot_bp_mean_mag_error_2	|	mag	|	BP-band mean magnitude error for B	|		|
|	161	|	phot_rp_mean_mag_2	|	mag	|	RP-band mean magnitude for B	|		|
|	162	|	phot_rp_mean_mag_error_2	|	mag	|	RP-band mean magnitude error for B	|		|
|	163	|	phot_bp_rp_excess_factor_2	|	mag	|	*BP*/*RP* excess factor for B	|		|
|	164	|	radial_velocity_2	|	km s-1	|	Radial velocity for B	|		|
|	165	|	radial_velocity_error_2	|	km s-1	|	Radial velocity error for B	|		|
|	166	|	Young_YMG	|	-	|	Young Moving Group identificator	|	‘Active’ are stars with literature references for activity; ‘Checked’ are overluminous stars in HR diagrams in which a conclusive reference for youth or activity has not been found	|
|	167	|	Young_ref	|	-	|	Young Moving Group reference	|		|
|	168	|	Class	|	-	|	Carmencita class	|	See Alonso Floriano et al. 2015 for further information	|
|	169	|	Sample_source	|	-	|	Source for the star	|	AF15: Alonso-Floriano et al. 2015b; Karmn: Carmencita (Alonso-Floriano et al. 2015b; Caballero et al. 2016); Kir91: Kirkpatrick et al. 1991; Lep13: Lépine et al. 2013; REC100: RECONS100 (Henry et al. 2006) ; Sma17: Smart et al. 2017	|
|	170	|	Bool_contam	|	Bool	|	Boolean index for contaminated stars	|	‘true’ for contaminated (DeltaG > 5 mag)	|
|	171	|	Bool_delta	|	Bool	|	Boolean index for close multiple stars	|	‘true’ for close multiplicity (rho <= 5 arcsec)	|
|	172	|	Bool_young	|	Bool	|	Boolean index for overluminous young stars	|	 ‘true’ for stars with youth features identified in the literature	|
|	173	|	Bool_dphot	|	Bool	|	Boolean index for stars with non-parallactic distances	|	 ‘true’ for photometric distances	|
|	174	|	Bool_RUWE	|	Bool	|	Boolean index for stars with high *Gaia* DR2 RUWE	|	‘true’ for RUWE > 1.41 (see Arenou et al. 2018; Lindegren et al. 2018) 	|
|	175	|	Bool_excess	|	Bool	|	Boolean index for stars with photometric flux excess in *Gaia* DR2	|	'true’ for stars with excess (see Evans et al. 2018)	|

**Notes:** 
- In all cases, quality flags are 'true' for poor photometric quality, and 'false' for good or acceptable photometric quality. For *2MASS* and *WISE*/*AllWISE* photometry, the alphabetic quality flags (of the form NNN and NNNN, respectively) are provided. In the *Gaia* DR2 data, the indicators `RUWE` and `phot_bp_rp_excess_factor determine the quality of their astrometric and photometric solution.
- Uncertainties in *Gaia* DR2 photometry are computed from the associated bolometric flux and its error.
- In the description, 'B' (not *B*) should be read 'the closest component in multiple systems'.

---

## Support

Reach out to me at <a href="mailto:ccifuentes@cab.inta-csic.es">`ccifuentes@cab.inta-csic.es`</a>.

---

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**

---

## Suggested Resources

- <a href="https://www.python.org/dev/peps/pep-0008/" target="_blank">Style Guide for Python Code (PEP 8)</a>
- <a href="https://carmenes.caha.es" target="_blank">CARMENES Website</a>
