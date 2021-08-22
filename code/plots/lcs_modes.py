# lc_modes.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurves modes

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import os
from .plotstyles import *
from mpl_axes_aligner import align


def lcs_modes(lcs, show=True, save=False, styling=False, color=False):
    print('lcs_modes: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    fig, ax1 = plt.subplots(figsize=(10, 2.5))

    # Title
    ax1.set_title(f'{lcs[0].name}')

    # Labels
    ax1.set_xlabel('Time (MJD)')
    ax1.set_ylabel('RXTE Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')
    
    MODES = {
        'PC': {
            'color': 'b',
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

    for row in lcs[0].ts:
        mode = row['mode']
        MODES[mode]['rows'].append(row)

    # SCALING_FACTOR = lcs[0].get_peak_row()['rate']/lcs[1].get_peak_row()['rate']
    # for i, row in enumerate(MODES['PC']['rows']):

    #     MODES['PC']['rows'][i]['rate'] = MODES['PC']['rows'][i]['rate'] * SCALING_FACTOR
    #     MODES['PC']['rows'][i]['rate_err_pos'] = MODES['PC']['rows'][i]['rate_err_pos'] * SCALING_FACTOR
    #     MODES['PC']['rows'][i]['rate_err_neg'] = MODES['PC']['rows'][i]['rate_err_neg'] * SCALING_FACTOR    

    # Plot used modes
    for mode_name, mode in MODES.items():
        # Search for used modes
        if len(mode['rows']) > 0:
            ax1.errorbar([], [], 
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
                ax1.errorbar([row['time'].mjd], [row['rate']], 
                    xerr=xerr, 
                    yerr=yerr, 
                    color=mode['color'], 
                    fmt=mode['marker'], 
                    fillstyle='none',
                    markeredgewidth=.5,
                    markersize=4,
                    elinewidth=.5)

    MODES = {
        'PC': {
            'color': 'b',
            'marker': 's',
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

    for row in lcs[1].ts:
        mode = row['mode']
        MODES[mode]['rows'].append(row)

    # SCALING_FACTOR = lcs[0].get_peak_row()['rate']/lcs[1].get_peak_row()['rate']
    # for i, row in enumerate(MODES['PC']['rows']):

    #     MODES['PC']['rows'][i]['rate'] = MODES['PC']['rows'][i]['rate'] * SCALING_FACTOR
    #     MODES['PC']['rows'][i]['rate_err_pos'] = MODES['PC']['rows'][i]['rate_err_pos'] * SCALING_FACTOR
    #     MODES['PC']['rows'][i]['rate_err_neg'] = MODES['PC']['rows'][i]['rate_err_neg'] * SCALING_FACTOR    

    # Plot used modes
    ax2 = ax1.twinx()
    for mode_name, mode in MODES.items():
        # Search for used modes
        if len(mode['rows']) > 0:
            ax1.errorbar([], [], 
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
                ax2.errorbar([row['time'].mjd], [row['rate']], 
                    xerr=xerr, 
                    yerr=yerr, 
                    color=mode['color'], 
                    fmt=mode['marker'], 
                    fillstyle='none',
                    markeredgewidth=.5,
                    markersize=4,
                    elinewidth=.5)

    ax2.set_ylabel('Swift Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')
  
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
    ax1.legend(shadow=False, edgecolor='white', fancybox=False, loc='upper left')
    align.yaxes(ax1, 0, ax2, 0, 0.22)
    # # Query errors
    # xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    # yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']

    # # Plot errorbars
    # plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt=' ', elinewidth=.5)   
    ax1.minorticks_on()
    ax2.minorticks_on()
    if save:
        filename = f'lc_joined_{lcs[0].name}.png'
        SAVE_FOLDER = f'output/report/lightcurves'

        # If saving directory does not exists make one
        if not os.path.isdir(SAVE_FOLDER):
            os.makedirs(SAVE_FOLDER)
            # print(f"NOTE: created this folder {SAVE_FOLDER}")

        plt.savefig(f'{SAVE_FOLDER}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    if show:
        plt.tight_layout()
        plt.show()