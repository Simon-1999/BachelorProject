# lc_exponential.py
# 
# Simon van Eeden

import matplotlib.pyplot as plt
from .plotstyles import *

def lc_exponential(lc, ob, ob_decay, fit, fit_tdecay, show=True, save=False, styling=False):
    print('lc_gaussian: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title(f'Outburst of {lc.name} - {lc.telescope}')

    # Axis
    plt.xlabel('Time (MJD)')
    plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')

    # Light curve
    xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
    yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']
    plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, color='grey', elinewidth=.5, label=f'Quiescent state')  

    # Outburst
    xerr = ob.ts['time_err_pos'], ob.ts['time_err_pos']
    yerr = ob.ts['rate_err_neg'], ob.ts['rate_err_neg']
    plt.errorbar(ob.ts.time.mjd, ob.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, color='k', elinewidth=.5, label=f'Outburst region')  

    # Decay region
    plt.plot(ob_decay.ts.time.mjd, ob_decay.ts['rate'], 's', color='b', label=f'Decay region')  

    # Fit of length
    plt.plot(fit.ts['time'].mjd, fit.ts['rate'], '-', color='red', lw=1, label=fit.name)

    # Add decay to legend
    plt.plot([], [], ' ', label='$t_{decay}$ = ' + f'{abs(fit_tdecay):.2f}')

    # Set x limit
    plt.xlim(ob.ts.time.mjd[0] - 25, ob.ts.time.mjd[-1] + 25)   

    # Legend
    plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)

    if save:
        plt.savefig(f'output/analysis/exponential fits obr/{self.name}_{self.telescope}_{ob.ts.time.mjd[0]:.0f}.png', dpi=150)

    if show:
        plt.show()

   