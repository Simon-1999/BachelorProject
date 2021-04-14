# main.py
#
# Simon van Eeden
#
# Running this code will excecute lightcurve analysis and fitting

# Packages
from astropy.time import Time
import json

# Internal
import processing
import vis
from classes import Lightcurve

# Loading light curve
lc_name = 'XTEJ1817-155'
lc_telescope = 'RXTE'
lc1 = Lightcurve(lc_name, lc_telescope)

# Get outburst region
ob1 = lc1.get_fraction('54230', '54540')

# Fit outburst
ob1.fit_gaussian(amplitude=70, amplitude_fixed=False, mean=54350, save=False)

# Bin lightcurve
# lc2 = lc1.binning(bin_days=10)

# Plot
# vis.lc_modes(lc1)
# vis.lcs_multi([ob1, ob2])
# vis.lcs_fill([lc1, lc2])
