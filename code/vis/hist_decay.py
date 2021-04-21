# hist_decay.py
# 
# Simon van Eeden
#
# This file contains functions for plotting distribution of decay times

import matplotlib.pyplot as plt
import matplotlib as mpl
from .plotstyles import *

def hist_decay(decays, styling=False, save=False, show=True):

    if styling:
        science_style()
    else:
        default_style()

    # Title
    plt.title('Outburst decay time distribution')

    # Histogram
    plt.hist(decays, 
        bins=50, 
        fill=False, 
        color='k',
        histtype ='step')

    # Labels
    plt.ylabel('Number of outbursts')
    plt.xlabel('Decay time (days)')

    if save:
        plt.savefig('test.png', dpi=250)
    
    if show:
        plt.show()
   
