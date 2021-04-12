# main.py
#
# Simon van Eeden
#
# Running this code

# Packages
from astropy.time import Time
import json

# Internal
import processing
import vis
from classes import Lightcurve

# Loading light curve
lc_name = 'XMMJ174457-2850.3'
lc_telescope = 'SwiftGC'
lc1 = Lightcurve(lc_name, lc_telescope)

# Get outburst region
# ob1 = lc1.get_fraction(start_time=Time('51250', format='mjd'), end_time=Time('51600', format='mjd'))

# Fit outburst
# ob1.fit_gaussian(amplitude=30, mean=51425)

# Plot lc
# vis.single_lc(lc1)

# Bin lightcurve
lc2 = lc1.binning(bin_days=10)

vis.multiple_lcs([lc1, lc2], save=True)
