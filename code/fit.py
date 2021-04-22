# fit.py
#
# Simon van Eeden
#
# Running this code will excecute lightcurve fitting

# Packages
from astropy.time import Time
import json
import numpy as np

# Internal imports
import plots
from classes import Lightcurve

# Loading light curve
lc_name = 'XTEJ1728-295'
lc_telescope = 'RXTE'
lc1 = Lightcurve(lc_name, lc_telescope)
# plots.lc_modes(lc1)

# Get outburst region
ob1 = lc1.get_fraction('55470', '55605')

# Fit outburst
# ob1.fit_gaussian(amplitude=.1, amplitude_fixed=True, mean=57450, save=True)
ob1.fit_exp(amplitude=1., tau=-1, save=False)

# Bin lightcurve
# lc2 = lc1.binning(bin_days=10)

# Plot multiple lightcurves
# plots.lcs_multi([lc1, lc1])
# plots.lcs_fill([lc1, lc1])
