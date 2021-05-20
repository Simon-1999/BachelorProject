# lc_gaussian.py
# 
# Simon van Eeden

import matplotlib.pyplot as plt
from .plotstyles import *

def lc_gaussian(lc, ob, fit_dur, fit, std, show=True, save=False, styling=False):
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

    # Fit of length
    plt.plot(fit_dur.ts['time'].mjd, fit_dur.ts['rate'], '-', color='red', lw=1, label=fit_dur.name)

    # Fit full
    plt.plot(fit.ts['time'].mjd, fit.ts['rate'], '--', color='red', lw=.5, label=fit.name)

    # Add durations to legend
    dur_gauss = std*6
    dur_region = ob.ts.time.mjd[-1] - ob.ts.time.mjd[0]
    plt.plot([], [], ' ',label='$t_{ob}$ gauss = ' + f'{dur_gauss:.2f}')
    plt.plot([], [], ' ',label='$t_{ob}$ outburst region = ' + f'{dur_region:.2f}')

    # Set x limit
    plt.xlim(ob.ts.time.mjd[0] - 25, ob.ts.time.mjd[-1] + 25)   

    # Legend
    plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)

    if save:
        plt.savefig(f'output/analysis/gaussian fits obr/{self.name}_{self.telescope}_{ob.ts.time.mjd[0]:.0f}.png', dpi=150)

    if show:
        plt.show()
   