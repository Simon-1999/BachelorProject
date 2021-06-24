# lc_gaussian.py
# 
# Simon van Eeden

import matplotlib.pyplot as plt
import numpy as np
import os
from .plotstyles import *

def lc_gaussian(lc, ob, fit_dur, fit, std, fit_y, show=True, save=False, styling=False):
    # print('lc_gaussian: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Setup figure
    fig = plt.figure(figsize=(6, 4.5))

    # Avoid minus sign from font
    plt.rc('axes', unicode_minus=False)

    # ------- Upper figure -------
    ax1 = plt.subplot2grid((4,1), (0,0), rowspan = 3, fig=fig)
    ax1.set_xticklabels([])
    ax1.minorticks_on()

    # Title
    ax1.set_title(f'{lc.name}')
    # plt.text(0.05, .9, f'{lc.name}', transform=plt.gca().transAxes) 

    # Peak of outburst
    ob_peak = ob.get_peak_row()['time'].mjd

    # Axis
    # ax1.xlabel('Time (MJD)')
    ax1.set_ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')

    # Fit of length
    # ax1.plot([fit_dur.ts['time'].mjd[0] - ob_peak], [fit_dur.ts['rate'][0]], '>', color='b', fillstyle='none', lw=1, label='Start', zorder=100, ms=7)
    # ax1.plot([fit_dur.ts['time'].mjd[-1] - ob_peak], [fit_dur.ts['rate'][-1]], '<', color='b', fillstyle='none', lw=1, label='End', zorder=100, ms=7)

    MODES = {
        'PC': {
            'color': 'k',
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
            'marker': 'o',
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
                ax1.errorbar([row['time'].mjd-ob_peak], [row['rate']], 
                    xerr=xerr, 
                    yerr=yerr, 
                    color=mode['color'], 
                    fmt=mode['marker'], 
                    fillstyle='none',
                    markeredgewidth=.5,
                    markersize=4,
                    elinewidth=.5,
                    alpha=.3)

    MODES = {
        'PC': {
            'color': 'k',
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
            'marker': 'o',
            'rows': [],
            'label': 'SwiftGC'
        },
    }

    # Classify datapoints to mode
    for row in ob.ts:
        mode = row['mode']
        MODES[mode]['rows'].append(row)

    # Plot used modes
    for mode_name, mode in MODES.items():
        # Search for used modes
        if len(mode['rows']) > 0:

            # Plot lightcurve of mode
            for row in mode['rows']:
                xerr = [row['time_err_neg']], [row['time_err_pos']]
                yerr = [row['rate_err_neg']], [row['rate_err_pos']]
                ax1.errorbar([row['time'].mjd-ob_peak], [row['rate']], 
                    xerr=xerr, 
                    yerr=yerr, 
                    color=mode['color'], 
                    fmt=mode['marker'], 
                    fillstyle='none',
                    markeredgewidth=.5,
                    markersize=4,
                    elinewidth=.5)

    # # Light curve
    # xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    # yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']
    # ax1.errorbar(lc.ts.time.mjd - ob_peak, lc.ts['rate'], 
    #     xerr=xerr, 
    #     yerr=yerr, 
    #     fmt='o', ms=5, 
    #     color='grey', 
    #     elinewidth=.5,
    #     fillstyle='none')  

    # Outburst
    # xerr = ob.ts['time_err_pos'], ob.ts['time_err_pos']
    # yerr = ob.ts['rate_err_neg'], ob.ts['rate_err_neg']
    # ax1.errorbar(ob.ts.time.mjd - ob_peak, ob.ts['rate'],
    #     xerr=xerr, 
    #     yerr=yerr, 
    #     fmt='o', 
    #     ms=5, 
    #     color='k', 
    #     elinewidth=.5, 
    #     markeredgewidth=.5,
    #     label=f'Fit-data',
    #     fillstyle='none')  

    # Fit full
    ax1.plot(fit.ts['time'].mjd - ob_peak, fit.ts['rate'], '-', color='k', lw=1)

    # Add durations to legend
    dur_gauss = std*6

    
    # dur_region = ob.ts.time.mjd[-1] - ob.ts.time.mjd[0]
    # ax1.plot([], [], ' ',label='$t_{dur}$ = ' + f'{dur_gauss:.2f}' + ' days')
    # plt.plot([], [], ' ',label='$t_{ob}$ outburst region = ' + f'{dur_region:.2f}')

    # Set limit
    length_ob_region = ob.ts.time.mjd[-1] - ob.ts.time.mjd[0]
    ax1.set_xlim(min([ob.ts.time.mjd[0], fit.ts['time'].mjd[0]]) - length_ob_region*.15 - ob_peak, max([ob.ts.time.mjd[-1], fit.ts['time'].mjd[-1]]) + length_ob_region*.15 - ob_peak)   
    ax1.set_ylim(min(ob.ts['rate']) - 1.1*max(ob.ts['rate_err_neg']), 1.1*max(ob.ts['rate']) + max(ob.ts['rate_err_pos'])) 

    # Legend
    ax1.legend(shadow=False, 
        edgecolor='white', 
        fancybox=False)

    # ------- Lower figure -------
    ax2 = plt.subplot2grid((4, 1), (3, 0), fig=fig)
    sig = np.ones(len(ob.ts.time.mjd))
    resid = ob.ts['rate'] - fit_y

    # Error on fit
    ax2.errorbar(ob.ts.time.mjd - ob_peak, resid/(ob.ts['rate_err_pos']), sig, 
        fmt='ok', 
        elinewidth=0.5, 
        ms=3)
    ax2.hlines(0, fit.ts['time'].mjd[0] - ob_peak, fit.ts['time'].mjd[-1] - ob_peak, lw=1, linestyle='-' ,color='k')

    # Axis labels
    ax2.set_xlabel(f'Time (days) - MJD = {ob_peak:.2f}')
    ax2.set_ylabel(r'$\frac{\mathregular{F}_{\mathregular{Obs}}-\mathregular{F}_{\mathregular{Model}}}{\sigma_\mathregular{F}}$', fontsize=16)


    # Axis limits
    # ax2.set_xlim(ob.ts.time.mjd[0] - 25 - ob_peak, ob.ts.time.mjd[-1] + 25 - ob_peak)
    ax2.set_xlim(min([ob.ts.time.mjd[0], fit.ts['time'].mjd[0]]) - length_ob_region*.15 - ob_peak, max([ob.ts.time.mjd[-1], fit.ts['time'].mjd[-1]]) + length_ob_region*.15 - ob_peak)   
    ax2.set_ylim(-5.6, 5.6)

    # Style
    ax2.minorticks_on()
    ax2.ticklabel_format(useOffset=False)

    # ax2.tick_params(axis='both', which='major', labelsize=12)
    # ax2.tick_params(axis='both', which='major', length=5)
    # ax2.tick_params(axis='both', which='minor', length=2.5)
    # ax2.tick_params(axis='both', which='both', direction='in', right=True, top=True)

    # Adjust whitespace between subplots
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.0)

    if save:
        filename = f'gauss_{ob.ts.time.mjd[0]:.0f}_{lc.name}_{lc.telescope}.png'
        SAVE_FOLDER = f'output/report/{lc.name}'

        # If saving directory does not exists make one
        if not os.path.isdir(SAVE_FOLDER):
            os.makedirs(SAVE_FOLDER)
            # print(f"NOTE: created this folder {SAVE_FOLDER}")
        
        # Check if filename is already used
        i = 1
        while os.path.exists(SAVE_FOLDER + filename + ".png"):
            # print(f"WARNING: The file {filename} already exists, auto filename has been used")
            filename = f'gauss_{ob_peaktime:.0f}_{lc.name}_{lc.telescope}({i}).png'
            i += 1

        plt.savefig(f'{SAVE_FOLDER}/{filename}', dpi=300, bbox_inches='tight')
        plt.close()

    if show:
        plt.show()
   