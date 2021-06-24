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
lc_std = 5.10
lc_avg = 4.78
lc1 = Lightcurve(lc_name, lc_telescope)
lc_telescope = 'Swift'
lc_std = 0.00
lc_avg = 0.00
lc2 = Lightcurve(lc_name, lc_telescope)
# lc1 = lc1.get_fraction(54644, 54657)
plots.lc_modes(lc1, styling=True)
# plots.lcs_modes([lc1, lc2], styling=True)

# Get background
# lc1.calc_background('55555', '56440', save=False)

# Get outburst region
# ob1_peaktime = '53897'
# ob1 = lc1.get_outburst(ob1_peaktime, std=lc_std)

# # Remove background
# ob1.correct_background(lc_avg)
# lc1.correct_background(lc_avg)

# # Gaussian fit
# save = False
# gaussian_fit, fit_stddev, fit_mean, fit_dur, fit_y = ob1.fit_gaussian()
# plots.lc_gaussian(lc1, ob1, fit_dur, gaussian_fit, fit_stddev, fit_y, save=save, styling=True)

# Exponential fit
# exponential_fit, fit_tdecay, ob1_decay, fit_y = ob1.fit_decay()
# plots.lc_exponential(lc1, ob1, ob1_decay, exponential_fit, fit_tdecay, fit_y, save=save, ylog=False, styling=True)

# Linear fit
# linear_fit, fit_tdecay, ob1_decay, fit_y = ob1.fit_ld()
# plots.lc_linear(lc1, ob1, ob1_decay, linear_fit, fit_tdecay, fit_y, save=save)

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
