# CARMENES input catalogue of M dwarfs. 
## V. Luminosities, colours, and spectral energy distributions.

<arXiv link available soon>
  
> This repository contains all pieces of code necessary to produce all figures and tables in Cifuentes et al. 2020 (<link>). 

[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![Publication](https://img.shields.io/badge/Published%3F-waiting-orange.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Ask Me Anything !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://GitHub.com/ccifuentesr)

---

## Table of Contents

- [Installation](#installation)
- [Structure](#structure)
- [Features](#features)
- [Contributing](#contributing)
- [Team](#team)
- [FAQ](#faq)
- [Support](#support)
- [License](#license)

---

## Installation

> The files are self-contained, self-consistent, homogenoeusly formatted, fairly self-explanatory.

- The code is provided as `*.py` files meant to be run individually.
- They can be run as Python Notebooks. The symbol `# %%` starts a cell that can be run separately.
- Some code requires of the data stored in folders stored in the repository.

### Clone

- Clone this repo to your local machine using `https://github.com/ccifuentesr/CARMENES-V`


## Structure

### Directories

- Directory ./<filename>.py: Includes all the code files classified as detailed below.
- Directory ./Baraffe: (see https://ui.adsabs.harvard.edu/abs/2015A%26A...577A..42B/)
- Directory ./Filters: Includes the data to reproduce the transmission curves for 20 filters (FUV to W4). 
- Directory ./Literature: Includes the individual tables with data from the literature.
- Directory ./Spectra: Includes synthetic spectra...
- Directory ./Stars: Includes data for individual stars to produce SED in some of the diagrams.

### Files

All files are named `cif20.xxx_yyy_zzz.py`, where `xxx` define the kind of output that it produces, `yyy` gives additional information about the output, and `zzz` enumerates the main variables involved. For example, the script `cif20.plot_literature_Mabs_SpT.py` produces an absolute magnitude vs. spectral type plot, and compares the values with those of the literature. The complete list of files and their description can be found below.

| File | Description | Input<sup id="a1">[1](#f1)</sup>| Output | 
| --- | --- | --- | 
| cif20.calculator_averagecolors.py | Average colours and standard deviations for each spectral type from K5V to L8. | Magnitudes from FUV to W4. | Table A.2 |
| cif20.calculator_Karmn.py | A generator of Carmencita identification name (Karmn JHHMMm+DDdAAA). | Simbad valid name. | Karmn ID |
|	cif20.calculator_MR.py	|	Calculator of mass and radius in solar units.	|	Lbol, Teff.	|	... |	
|	cif20.calculator_completeness.py	|	Estimator of the completeness in volume of the sample.	|	d	|	-	|	Completeness of the sample (ratio).
| ... | ... | ... | ... |
| ... | ... | ... | ... |
| ... | ... | ... | ... |
| ... | ... | ... | ... |
| ... | ... | ... | ... |
| ... | ... | ... | ... |

1. <small id="f1"> Lbol = bolometric luminosity (Lsol), Teff = effective temperature (K), M = mass (Msol), R = radius (Rsol), d = distance (pc). Uncertainties are almost always used as an input, but they are omitted here for simplicity. </small> [↩](#a1) 

Note: 
---

## FAQ

- **FAQ questions?**
    - Under construction.

---

## Support

Reach out to me:

- E-mail at <a href="mailto:ccifuentes@cab.inta-csic.es">`ccifuentes@cab.inta-csic.es`</a>
- Website at <a href="http://aplaceformyhead.es" target="_blank">`aplaceformyhead.es`</a>

---

## License

[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

- **[MIT license](http://opensource.org/licenses/mit-license.php)**

## Suggested resources

- <a href="https://www.python.org/dev/peps/pep-0008/" target="_blank">Style Guide for Python Code (PEP 8)</a>

---

[![INSERT YOUR GRAPHIC HERE](http://i.imgur.com/dt8AUb6.png)]()
> GIF Tools

- Use <a href="http://recordit.co/" target="_blank">**Recordit**</a> to create quicks screencasts of your desktop and export them as `GIF`s.
- For terminal sessions, there's <a href="https://github.com/chjj/ttystudio" target="_blank">**ttystudio**</a> which also supports exporting `GIF`s.

**Recordit**

![Recordit GIF](http://g.recordit.co/iLN6A0vSD8.gif)

**ttystudio**

![ttystudio GIF](https://raw.githubusercontent.com/chjj/ttystudio/master/img/example.gif)


- For all the possible languages that support syntax highlithing on GitHub (which is basically all of them), refer <a href="https://github.com/github/linguist/blob/master/lib/linguist/languages.yml" target="_blank">here</a>.

