# main.py
#
# Simon van Eeden
#
# Running this code will excecute lightcurve analysis and fitting

# Packages
from astropy.time import Time
import json

# Internal imports
import vis
from classes import Lightcurve

# Loading light curve
lc_name = 'GRS1741-2853'
lc_telescope = 'SwiftGC'
lc1 = Lightcurve(lc_name, lc_telescope)

# Get outburst region
ob1 = lc1.get_fraction('57503', '57520')

# Fit outburst
# ob1.fit_gaussian(amplitude=.2, amplitude_fixed=False, mean=58050, save=True)
# ob1.fit_exp(amplitude=1., tau=-1, save=True)

# Duration distrubution
# durations = [10, 1.21, 1.35, 1.56, 2.93, 3.36, 8.94, 21.33, 6.57, 5.86, 6.25, 5.87, 9, 4, 3.41, 25.47, 2.27, 3.36, 12.30, 4.13, 5.83, 3.31, 1.77, 2.80, 2.40, 2.70, 46.08, 8.05, 2.62, 3.15, 8.94, 6.89, 13.26, 5.55, 3.06, 5.77, 4.10]
# for i, duration in enumerate(durations):
#     durations[i] = durations[i] * 6
# vis.hist_duration(durations)

# Decay distribution
decays = [7.79, 1.20, 3.26, 1.47, 2.28, 5.21, 3.38, 8.94, 3.60, 4.19, 2.62, 8.46, 42.39, 4.51, 1.50, 2.57, 2.56, 1.68, 2.72, 36.93, 31.49, 3.90, 4.31, 11.26, 11.06, 6.48, 30.58, 4.66, 5.67, 7.83, 5.15, 6.71, 83.91, 13.08, 31.36, 12.12, 5.96, 4.11, 2.13, 2.76, 19.57]
vis.hist_decay(decays)

# Bin lightcurve
# lc2 = lc1.binning(bin_days=10)

# Plot
# vis.lc_modes(lc1)

# vis.lcs_multi([lc1, lc1])
# vis.lcs_fill([lc1, lc1])
