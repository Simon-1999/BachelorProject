# lc_modes.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurve modes

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
from .plotstyles import *


def lc_modes(lc, show=True, save=False, styling=False, color=False):
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

    MM = {
        'PC': 's',
        'PCA': 'D',
        'ASM': 'D',
        'WT': 'v',
        'SwiftGC': 'o',
    }

    # Plot time vs rate
    for row in lc.ts:  
        if color:
            plt.plot(row['time'].mjd, row['rate'], 's', markersize=3, color=CM[row['mode']])
        else:
            plt.plot(row['time'].mjd, row['rate'], MM[row['mode']], markersize=4, color='k', fillstyle='none', alpha=0.5)

    # Add labels
    patches = []
    for mode in CM.keys():
        if mode in lc.ts['mode']:
            if color:
                plt.plot([],[], 's', ms=3, color=CM[mode], label=mode)
            else:
                plt.plot([],[], MM[mode], ms=4, color='k', fillstyle='none', label=mode)
            # patches.append(mpatches.Patch(color=CM[mode], marker=MM[mode], label=f'{mode}'))
    plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)

    # Query errors
    xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']

    # Plot errorbars
    plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt=' ', elinewidth=.5)   

    if save:
        plt.savefig(f'output/lightcurves/single_lc_{lc.name}_{lc.telescope}_1.png', dpi=250)

    if show:
        plt.show()