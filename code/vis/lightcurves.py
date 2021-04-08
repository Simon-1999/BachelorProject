# lightcurves.py
# 
# Simon van Eeden
#
# This file contains functions for plotting lightcurves

import matplotlib.pyplot as plt

def single_lc(lc, show=True, save=False, styling=False):
    print('single_lc: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title(f'Lightcurve {lc.name} - {lc.telescope}')

    # Labels
    plt.xlabel('$Time$ $(MJD)$')
    plt.ylabel('$Rate$')

    # Plot errors
    lc.plot_lc()
     
    if save:
        plt.savefig(f'output/lightcurves/single_lc_{lc.name}_{lc.telescope}.png', dpi=150)

    if show:
        plt.show()


def multiple_lcs(lcs, show=True, save=False, styling=False):
    print('multiple_lcs: Plotting...')

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title(f'Lightcurves')

    # Axis
    plt.xlabel('$Time$ $(MJD)$')
    plt.ylabel('$Rate$')

    # Light curves
    for lc in lcs:
        lc.plot_lc_clean()

    # Legend
    plt.legend()

    if save:
        plt.savefig(f'output/lightcurves/multiple_lcs.png', dpi=200)

    if show:
        plt.show()


def default_style():

    plt.style.use('default')
    plt.rcParams['xtick.major.size'] = 5.0
    plt.rcParams['xtick.minor.size'] = 3.0
    plt.rcParams['ytick.major.size'] = 5.0
    plt.rcParams['ytick.minor.size'] = 3.0
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'


def science_style():

    plt.style.use('default')
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.size'] = 15
    plt.rcParams['legend.fontsize'] = 18
    plt.rcParams['xtick.major.size'] = 5.0
    plt.rcParams['xtick.minor.size'] = 3.0
    plt.rcParams['ytick.major.size'] = 5.0
    plt.rcParams['ytick.minor.size'] = 3.0
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
