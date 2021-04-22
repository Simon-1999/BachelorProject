# lc_modes.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurve modes

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from .plotstyles import *


def lc_modes(lc, show=True, save=False, styling=False):
    print('lc_modes: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title(f'Lightcurve {lc.name} - {lc.telescope}')

    # Labels
    plt.xlabel('Time (MJD)')
    plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')

    # Plot lc
    CM = {
        'PC': 'r',
        'PCA': 'b',
        'ASM': 'r',
        'WT': 'b',
        'SwiftGC': 'g',
    }

    # Plot time vs rate
    for row in lc.ts:  
        plt.plot(row['time'].mjd, row['rate'], 's', markersize=3, color=CM[row['mode']])

    # Add labels
    patches = []
    for mode in CM.keys():
        if mode in lc.ts['mode']:
            patches.append(mpatches.Patch(color=CM[mode], label=f'{mode}'))
    plt.legend(handles=patches, shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)

    # Query errors
    xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']

    # Plot errorbars
    plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt=' ', elinewidth=.5)   

    if save:
        plt.savefig(f'output/lightcurves/single_lc_{lc.name}_{lc.telescope}_1.png', dpi=250)

    if show:
        plt.show()