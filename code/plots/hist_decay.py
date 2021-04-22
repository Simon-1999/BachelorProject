# hist_decay.py
# 
# Simon van Eeden
#
# This file contains functions for plotting distribution of decay times

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from .plotstyles import *
import matplotlib.ticker as ticker

def hist_decay(decays, bins=20, styling=False, save=False, show=True, log=True):

    if styling:
        science_style()
    else:
        default_style()

    # One ticker per count
    axes = plt.gca()
    axes.yaxis.set_minor_locator(ticker.MultipleLocator(1))

    # Title
    plt.title('Outburst decay time distribution')

    # Histogram
    if log:
        bin_list = 10**np.linspace(0, len(str(int(max(decays)))), bins)
    else:
        bin_list = bins
        
    plt.hist(decays, 
        bins=bin_list, 
        fill=False, 
        color='k')

    # Axis
    if log:
        plt.xscale('log')

    # Labels
    plt.ylabel('Number of outbursts')
    plt.xlabel('Decay timescale (days)')

    if save:
        plt.savefig(f'output/analysis/distribution/Histogram_decaytime_bins{bins}.png', dpi=250)
    
    if show:
        plt.show()
   
