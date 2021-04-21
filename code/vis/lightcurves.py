# lightcurves.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurves

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
            patches.append(mpatches.Patch(color=CM[mode], label=f'{mode} mode'))
    plt.legend(handles=patches, shadow=False, edgecolor='k')

    # Query errors
    xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']

    # Plot errorbars
    plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt=' ', elinewidth=.5)   

    if save:
        plt.savefig(f'output/lightcurves/single_lc_{lc.name}_{lc.telescope}_1.png', dpi=250)

    if show:
        plt.show()


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

   
