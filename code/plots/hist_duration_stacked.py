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

def hist_duration_stacked(durations, labels, bins=20, styling=False, save=False, show=True, log=True):

    if styling:
        science_style()
    else:
        default_style()
        
    fig = plt.figure(figsize=(6, 4.5))

    # One ticker per count
    axes = plt.gca()
    axes.yaxis.set_minor_locator(ticker.MultipleLocator(1))

    # Title
    plt.title('Outburst duration distribution')

    colors = []
    color_idx = np.linspace(0, 1, len(durations))
    for i in color_idx:
        colors.append(plt.cm.jet(i))
    colors = ['b', 'r', 'lightgrey']

    # Histogram
    if log:
        bin_list = 10**np.linspace(0.7, 2.8, bins)
    else:
        bin_list = bins
        
    # Similar to decays
    decays = durations

    # print("BH")
    # print(f"Average = {np.average(decays[1])}")
    # print(f"Median = {np.median(decays[1])}")
    # print("NS")
    # print(f"Average = {np.average(decays[0])}")
    # print(f"Median = {np.median(decays[0])}")

    plt.hist(decays[0],
        histtype='step', 
        bins=bin_list, 
        edgecolor='b',
        linestyle='--',
        lw=2,
        facecolor='none',
        label='NS',
        alpha=.9)

    plt.hist(decays[1], 
        histtype='step', 
        bins=bin_list, 
        edgecolor='r',
        linestyle=':',
        lw=2,
        facecolor='none',
        label='BH',
        alpha=.9)

    plt.hist(decays[2], 
        histtype='step', 
        bins=bin_list, 
        edgecolor='k',
        linestyle='-',
        lw=2,
        facecolor='none',
        label='?',
        alpha=.3)

    plt.hist(decays[0]+decays[1]+decays[2], 
        histtype='step', 
        bins=bin_list, 
        edgecolor='k',
        linestyle='-',
        lw=1,
        facecolor='none',
        label='All',
        alpha=1)

    plt.legend(shadow=False, edgecolor='white', fancybox=False)

    # Axis
    if log:
        plt.xscale('log')

    # Labels
    plt.ylabel('Number of outbursts')
    plt.xlabel('Duration (days)')

    if save:
        plt.savefig(f'output/analysis/distribution/Histogram_duration_BH-NS_bins{bins}.png', dpi=250)

    if show:
        plt.show()   
