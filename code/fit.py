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
lc_std = 0.00
lc_avg = 0.00
lc1 = Lightcurve(lc_name, lc_telescope)

# Loading another light curve
lc_telescope = 'Swift'
lc_std = 0.00
lc_avg = 0.00
lc2 = Lightcurve(lc_name, lc_telescope)

# Plot lightcurves from same source with different telescopes
plots.lcs_modes([lc1, lc2], styling=True, save=False)

# Plot lightcurves from any source
plots.lcs_multi([lc1, lc2], styling=True, save=False)

# Plot one lightcurve
plots.lc_modes(lc1, styling=True, save=False)

# Calculate background
lc1.calc_background('54000', '55000', save=False)

# Get fraction of lightcurve
lc1_frac = lc1.get_fraction(55200, 55500)

# Bin lightcurve
lc1_binned = lc1.binning(bin_days=20)

# Get outburst region
ob1_peaktime = '55462'
ob1 = lc1.get_outburst(ob1_peaktime, std=lc_std)

# Remove background
ob1.correct_background(lc_avg)
lc1.correct_background(lc_avg)

# Gaussian fit
gaussian_fit, fit_stddev, fit_mean, fit_dur, fit_y = ob1.fit_gaussian()
plots.lc_gaussian(lc1, ob1, fit_dur, gaussian_fit, fit_stddev, fit_y, save=False, styling=True)
print(f'ob_dur = {fit_stddev*6}')

# Exponential fit
exponential_fit, fit_tdecay, ob1_decay, fit_y = ob1.fit_decay()
plots.lc_exponential(lc1, ob1, ob1_decay, exponential_fit, fit_tdecay, fit_y, save=False, ylog=False, styling=True)
print(f'ob_tdecay = {fit_tdecay}')

# Linear fit
linear_fit, fit_tdecay, ob1_decay, fit_y = ob1.fit_ld()
plots.lc_linear(lc1, ob1, ob1_decay, linear_fit, fit_tdecay, fit_y, save=False, styling=True)
