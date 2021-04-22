# hist_duration.py
# 
# Simon van Eeden
#
# This file contains functions for plotting distrubtion of duration

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from .plotstyles import *
import matplotlib.ticker as ticker

def hist_duration(durations, bins=20, styling=False, save=False, show=True, log=True):

    if styling:
        science_style()
    else:
        default_style()
        
    # One ticker per count
    axes = plt.gca()
    axes.yaxis.set_minor_locator(ticker.MultipleLocator(1))

    # Title
    plt.title('Outburst duration distribution')

    # Histogram
    if log:
        bin_list = 10**np.linspace(0, len(str(int(max(durations)))), bins)
    else:
        bin_list = bins
        
    plt.hist(durations, 
        bins=bin_list, 
        fill=False, 
        color='k')

    # Axis
    if log:
        plt.xscale('log')

    # Labels
    plt.ylabel('Number of outbursts')
    plt.xlabel('Duration (days)')

    if save:
        plt.savefig(f'output/analysis/distribution/Histogram_duration_bins{bins}.png', dpi=250)

    if show:
        plt.show()
   
