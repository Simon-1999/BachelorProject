# lc_fits.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurve with multiple fits

import matplotlib.pyplot as plt
from .plotstyles import *

def lc_fits(lc, ob, fits, show=True, save=False, styling=False):
    print('lcs_multi: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title('Lightcurves')

    # Axis
    plt.xlabel('Time (MJD)')
    plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')

    # Light curve
    xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']
    plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, color='grey', elinewidth=.5, label=f'Background')  

    # Outburst
    xerr = ob.ts['time_err_pos'], ob.ts['time_err_pos']
    yerr = ob.ts['rate_err_neg'], ob.ts['rate_err_neg']
    plt.errorbar(ob.ts.time.mjd, ob.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, color='k', elinewidth=.5, label=f'Outburst')  

    plt.xlim(ob.ts.time.mjd[0] - 25, ob.ts.time.mjd[-1] + 25)

    # Fits
    for fit in fits:
        plt.plot(fit.ts['time'].mjd, fit.ts['rate'], '--', color='gray', label=fit.name)

    # Legend
    plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)

    if save:
        plt.savefig(f'output/analysis/lcs_multi.png', dpi=200)

    if show:
        plt.show()

   