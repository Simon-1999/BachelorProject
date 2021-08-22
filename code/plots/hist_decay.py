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

def hist_decay(decays, labels, bins=20, styling=False, save=False, show=True, log=True):

    if styling:
        science_style()
    else:
        default_style()

    fig = plt.figure(figsize=(6, 4.5))

    # One ticker per count
    axes = plt.gca()
    axes.yaxis.set_minor_locator(ticker.MultipleLocator(1))

    # Title
    plt.title('Outburst decay time distribution')

    colors = []
    color_idx = np.linspace(0, 1, len(decays))
    for i in color_idx:
        colors.append(plt.cm.jet(i))

    # Histogram
    if log:
        bin_list = 10**np.linspace(0, 2.1, bins)
    else:
        bin_list = bins
        
    plt.hist(decays, 
        bins=bin_list, 
        stacked=True,
        color=colors,
        label=labels)
    plt.legend(shadow=False, edgecolor='white', fancybox=False)

    # Axis
    if log:
        plt.xscale('log')

    # Labels
    plt.ylabel('Number of outbursts')
    plt.xlabel('Decay timescale (days)')

    if save:
        plt.savefig(f'output/analysis/distribution/Histogram_decaytime_sources_bins{bins}.png', dpi=250)
    
    if show:
        plt.show()
   
