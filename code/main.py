# main.py
#
# Simon van Eeden
#
# Running this code

# Packages
from astropy.time import Time

# Internal
import processing
import vis
from classes import Lightcurve

lc1 = Lightcurve('XTEJ1734-234', 'RXTE', 'data/raw/RXTE/XTEJ1734-234/XTEJ1734_lc.txt')
ob1 = lc1.fraction(start_time=Time('51250', format='mjd'), end_time=Time('51600', format='mjd'))
ob1.gaussian_fit(amplitude=30, mean=51400)

# Plot lc
# vis.single_lc(ob1)

# Bin lightcurve
# lc2 = copy.deepcopy(lc1)
# lc2 = lc1.binning(bin_days=10)

# vis.multiple_lcs([lc1, lc2])
