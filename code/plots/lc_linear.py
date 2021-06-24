# lc_linear.py
# 
# Simon van Eeden

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from .plotstyles import *

def lc_linear(lc, ob, ob_decay, fit, fit_tdecay, fit_y, show=True, save=False, styling=False, ylog=False):
    print('lc_exponential: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Setup figure
    fig = plt.figure(figsize=(8, 6))

    # Avoid minus sign from font
    plt.rc('axes', unicode_minus=False)

    # ------- Upper figure -------
    ax1 = plt.subplot2grid((4,1), (0,0), rowspan = 3, fig=fig)
    # ax1.set_ylabel(r'Count rate (c s$^{-1}$)')
    ax1.set_xticklabels([])
    ax1.minorticks_on()
    # ax1.tick_params(axis='both', which='major', labelsize=12)
    # ax1.tick_params(axis='both', which='major', length=5)
    # ax1.tick_params(axis='both', which='minor', length=2.5)
    # ax1.tick_params(axis='both', which='both', direction='in', right=True, top=True)

    # Peak of outburst
    ob_peak = ob_decay.ts.time.mjd[0]

    # Title
    ax1.set_title(f'Outburst of {lc.name} - {lc.telescope}')

    # Axis label
    # ax1.set_xlabel('Time (MJD)')
    ax1.set_ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')

    # Light curve
    xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']
    ax1.errorbar(lc.ts.time.mjd - ob_peak, lc.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, color='grey', elinewidth=.5, capsize=1.5, label=f'Quiescent state', alpha=.5)  

    # Outburst
    xerr = ob.ts['time_err_pos'], ob.ts['time_err_pos']
    yerr = ob.ts['rate_err_neg'], ob.ts['rate_err_neg']
    ax1.errorbar(ob.ts.time.mjd - ob_peak, ob.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, color='k', elinewidth=.5, capsize=1.5, label=f'Outburst region')  

    # Decay region
    # ax1.plot(ob_decay.ts.time.mjd - ob_peak, ob_decay.ts['rate'], 's', ms=3, color='b', label=f'Decay region')  

    # Linear fit and decay time
    ax1.plot(fit.ts['time'].mjd - ob_peak, fit.ts['rate'], '--', color='grey', lw=1, label=fit.name)
    ax1.plot([], [], ' ', label='$t_{decay}$ = ' + f'{abs(fit_tdecay):.2f}')

    # Axis limit
    ax1.set_xlim(ob.ts.time.mjd[0] - 2 - ob_peak, ob.ts.time.mjd[-1] + 5 - ob_peak)

    # Log scale
    if ylog:
        ax1.set_yscale('log')
    # plt.ylim(min(ob.ts['rate']) - 2, max(ob.ts['rate']) + 2)  

    # ------- Lower figure -------
    ax2 = plt.subplot2grid((4, 1), (3, 0), fig=fig)
    sig = np.ones(len(ob_decay.ts.time.mjd[10:-1]))
    resid = ob_decay.ts['rate'][10:-1] - fit_y

    # Error on fit
    ax2.errorbar(ob_decay.ts.time.mjd[10:-1] - ob_peak, resid/(ob_decay.ts['rate_err_pos'][10:-1]), sig, fmt='.k', elinewidth=0.5, capsize=1.5, ms=3)
    ax2.hlines(0, 0, fit.ts['time'].mjd[-1] - ob_peak, linewidth=1, linestyle='--', color='grey')

    # Axis labels
    ax2.set_xlabel(f'Time after outburst peak (days) - MJD = {ob_peak:.2f}')
    ax2.set_ylabel(r'$\frac{F_{Obs}-F_{Model}}{\sigma_F}$', fontsize=14)

    # Axis limits
    ax2.set_xlim(ob.ts.time.mjd[0] - 2 - ob_peak, ob.ts.time.mjd[-1] + 5 - ob_peak)
    ax2.set_ylim(-5.6, 5.6)

    # Style
    ax2.minorticks_on()
    ax2.ticklabel_format(useOffset=False)
    # ax2.tick_params(axis='both', which='major', labelsize=12)
    # ax2.tick_params(axis='both', which='major', length=5)
    # ax2.tick_params(axis='both', which='minor', length=2.5)
    # ax2.tick_params(axis='both', which='both', direction='in', right=True, top=True)

    # Legend
    ax1.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)

    # Adjust whitespace between subplots
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.0)

    if save:
        plt.savefig(f'output/analysis/exponential fits obr/{lc.name}_{lc.telescope}_{ob.ts.time.mjd[0]:.0f}.png', dpi=150)

    if show:
        plt.show()

   