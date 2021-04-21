# lc_fill.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurves filling errorbars

import matplotlib.pyplot as plt
from .plotstyles import *

def lcs_fill(lcs, show=True, save=False, styling=False):
    print('lcs_fill: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title('Lightcurves')

    # Labels
    plt.xlabel('Time (MJD)')
    plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')

    # Plot lcs
    ax = plt.gca()
    for lc in lcs: 
        color = next(ax._get_lines.prop_cycler)['color']
        if lc.ts.__class__.__name__ == 'BinnedTimeSeries':
            plt.plot(lc.ts['time_bin_start'].mjd, lc.ts['rate'], '-',drawstyle='steps-post', label=lc.name)

        if lc.ts.__class__.__name__ == 'TimeSeries':
            plt.plot(lc.ts.time.mjd, lc.ts['rate'], '-', color=color, lw=0.5, label=lc.name)
            plt.fill_between(lc.ts.time.mjd, lc.get_rate_upper(), lc.get_rate_lower(), color=color, alpha=.3)  

    # Legend
    plt.legend(shadow=False, edgecolor='k')

    if save:
        plt.savefig(f'output/analysis/lcs_fill.png', dpi=150)

    if show:
        plt.show()