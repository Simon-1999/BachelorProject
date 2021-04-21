# hist_duration.py
# 
# Simon van Eeden
#
# This file contains functions for plotting distrubtion of duration

import matplotlib.pyplot as plt
import matplotlib as mpl
from .plotstyles import *

def hist_duration(durations, styling=False, save=False, show=True):

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title('Outburst duration distribution')

    # Histogram
    plt.hist(durations, 
        bins=20, 
        fill=False, 
        color='k',
        histtype = 'step')

    # Labels
    plt.ylabel('Number of outbursts')
    plt.xlabel('Duration (days)')

    if save:
        plt.savefig('test.png', dpi=250)

    if show:
        plt.show()
   
