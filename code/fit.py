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
lc_name = 'SAXJ1753.5-2349'
lc_telescope = 'RXTE'
lc_std = 4.70637
lc_avg = -0.25799
lc1 = Lightcurve(lc_name, lc_telescope)
# plots.lc_modes(lc1)

# Get background
# lc1.calc_background('55555', '56440', save=False)

# Get outburst region
ob1_peaktime = '51392'
ob1 = lc1.get_outburst(ob1_peaktime, std=lc_std)

# Remove background
ob1.correct_background(lc_avg)
lc1.correct_background(lc_avg)

# Fits
gaussian_fit, fit_stddev, fit_mean, fit_dur = ob1.fit_gaussian(ob1_peaktime)
plots.lc_gaussian(lc1, ob1, fit_dur, gaussian_fit, fit_stddev)

exponential_fit, fit_tdecay, ob1_decay = ob1.fit_decay()
plots.lc_exponential(lc1, ob1, ob1_decay, exponential_fit, fit_tdecay)

# ob1 = lc1.get_fraction('53500', '55000')
# print(f"Stddev = {np.std(ob1.ts['rate'])}")
# print(f"peak rate = {max(ob1.ts['rate'])}")
# print(f"upperlimit F10 = {0.2*max(ob1.ts['rate']):2}")
# print(f"F10 = {0.1*max(ob1.ts['rate']):2}")
# print(f"lowerlimit F10 = {0.01*max(ob1.ts['rate']):2}")

# Fit outburst
# ob1.fit_gaussian(amplitude=40, amplitude_fixed=False, mean=52510)

# ob1.fit_exp(amplitude=1., tau=-1, save=False)

# Bin lightcurve
# lc2 = lc1.binning(bin_days=60)
# lc3 = lc1.get_outbursts(std=std)
# print(f"Average = {np.average(lc1.ts['rate'])}")
# print(f"Median = {np.median(lc1.ts['rate'])}")
# print(f"Stddev = {np.std(lc1.ts['rate'])}")

# Plot multiple lightcurves
# plots.lcs_multi([lc1, lc3])
# plots.lcs_fill([lc1, lc1])
