# lc_modes.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurve modes

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import os
from .plotstyles import *


def lc_modes(lc, show=True, save=False, styling=False, color=False):
    print('lc_modes: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    fig = plt.figure(figsize=(5, 3.75))

    # Title
    plt.title(f'{lc.name}')

    # Labels
    plt.xlabel('Time (MJD)')
    plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')

    MODES = {
        'PC': {
            'color': 'k',
            'marker': 'D',
            'rows': [],
            'label': 'Swift PC'
        },
        'PCA': {
            'color': 'k',
            'marker': 'D',
            'rows': [],
            'label': 'RXTE PCA'
        },
        'ASM': {
            'color': 'k',
            'marker': 'D',
            'rows': [],
            'label': 'RXTE ASM'
        },
        'WT': {
            'color': 'b',
            'marker': 'v',
            'rows': [],
            'label': 'Swift WT'
        },
        'SwiftGC': {
            'color': 'k',
            'marker': 'D',
            'rows': [],
            'label': 'SwiftGC'
        },
    }

    # Classify datapoints to mode
    for row in lc.ts:
        mode = row['mode']
        MODES[mode]['rows'].append(row)

    # Plot used modes
    for mode_name, mode in MODES.items():
        # Search for used modes
        if len(mode['rows']) > 0:
            plt.errorbar([], [], 
                xerr=1, 
                yerr=1, 
                color=mode['color'], 
                fmt=mode['marker'], 
                fillstyle='none',
                markeredgewidth=.5,
                markersize=4,
                elinewidth=.5, 
                label=mode['label']) 

            # Plot lightcurve of mode
            for row in mode['rows']:
                xerr = [row['time_err_neg']], [row['time_err_pos']]
                yerr = [row['rate_err_neg']], [row['rate_err_pos']]
                plt.errorbar([row['time'].mjd], [row['rate']], 
                    xerr=xerr, 
                    yerr=yerr, 
                    color=mode['color'], 
                    fmt=mode['marker'], 
                    fillstyle='none',
                    markeredgewidth=.5,
                    markersize=4,
                    elinewidth=.5)

    # # Plot time vs rate
    # for row in lc.ts:  
    #     if color:
    #         plt.plot(row['time'].mjd, row['rate'], 's', markersize=3, color=CM[row['mode']])
    #     else:
    #         plt.plot(row['time'].mjd, row['rate'], MM[row['mode']], markersize=4, color='k', fillstyle='none', alpha=0.5)

    # # Add labels
    # patches = []
    # for mode in CM.keys():
    #     if mode in lc.ts['mode']:
    #         if color:
    #             plt.plot([],[], 's', ms=3, color=CM[mode], label=mode)
    #         else:
    #             plt.plot([],[], MM[mode], ms=4, color='k', fillstyle='none', label=mode)
    #         # patches.append(mpatches.Patch(color=CM[mode], marker=MM[mode], label=f'{mode}'))
    plt.legend(shadow=False, edgecolor='white', fancybox=False)

    # # Query errors
    # xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    # yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']

    # # Plot errorbars
    # plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt=' ', elinewidth=.5)   
    plt.minorticks_on()
    if save:
        filename = f'lc_{lc.name}_{lc.telescope}.png'
        SAVE_FOLDER = f'output/report/{lc.name}'

        # If saving directory does not exists make one
        if not os.path.isdir(SAVE_FOLDER):
            os.makedirs(SAVE_FOLDER)
            # print(f"NOTE: created this folder {SAVE_FOLDER}")

        plt.savefig(f'{SAVE_FOLDER}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    if show:
        plt.tight_layout()
        plt.show()