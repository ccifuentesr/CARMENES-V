# %%

# CALCULATOR: KARMN ID
# Cifuentes et al. 2020

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u


def karmnID(name):
    """Produces ID name for Carmencita stars.

    Args:
        name (string): Any name that can be resolved by Simbad 
        (e.g., 'GJ 51' or 'Proxima').

    Returns:
        string:  Identifier "JHHMMm+DDdAAA". 

    "N", "S", "E" or "W" is to be added manually if two stars 
    (in a close binary system) have the same "HHMMm+DDd" string. 
    The user is encouraged to visually inspect the star and its
    surroundings using Aladin. Gaia DR2 catalogue is capable of
    resolving close binaries, although spectroscopic binaries
    require of a more precise look.

    """
    target = SkyCoord.from_name(name)
    ra_hh = str(int(target.ra.hms[0])).rjust(2, '0')
    ra_mm = str(int(target.ra.hms[1])).rjust(2, '0')
    ra_ss = str(int(float(target.ra.hms[2])/6))
    de_mm = str(abs(int(float(target.dec.dms[1])/6)))
    if target.dec.dms[1] < 0:
        de_dd = f'{int(target.dec.dms[0]):03}'
        print('J%s%s%s%s%s' % (ra_hh, ra_mm, ra_ss, de_dd, de_mm))
    else:
        de_dd = int(target.dec.dms[0])
        print('J%s%s%s%+d%s' % (ra_hh, ra_mm, ra_ss, de_dd, de_mm))
