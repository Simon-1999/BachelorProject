# lcs_multi.py
# 
# Simon van Eeden
#
# This file contains functions for plotting multiple lightcurves including binned curves

import matplotlib.pyplot as plt
from .plotstyles import *

def lcs_multi(lcs, show=True, save=False, styling=False):
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

    # Light curves
    for lc in lcs:
        if lc.ts.__class__.__name__ == 'BinnedTimeSeries':
            plt.plot(lc.ts['time_bin_start'].mjd, lc.ts['rate'], '-', drawstyle='steps-post', label=f'{lc.name} - {lc.telescope}')

        if lc.ts.__class__.__name__ == 'TimeSeries':
            # Query errors
            xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
            yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']
            plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, elinewidth=.5, label=f'{lc.name} - {lc.telescope}')  

    # Legend
    plt.legend(shadow=False, edgecolor='k')

    if save:
        plt.savefig(f'output/analysis/lcs_multi.png', dpi=200)

    if show:
        plt.show()

   