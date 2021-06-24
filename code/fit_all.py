# fit_all.py
#
# Simon van Eeden
#
# Running this code will excecute gaussian fit and exponential fits

# Packages
from astropy.time import Time
import json
import numpy as np

# Internal imports
import plots
from classes import Lightcurve

# Open outburst file
file = 'data/outbursts/Outbursts_fitparam.tsv'
read_file = open(file,'r')
plotted_telescopes = []

# Loop over file
for l in read_file:

    # Remove spaces
    l = l.strip('\n')

    # Unread header        
    if not l[0] == '#':

        # Get columns
        row = l.split('\t')
        read_source = row[0] != ''
        read_obs = row[6] != ''

        if read_source:  
            # Read source name
            plotted_telescopes = []
            lc_name = row[0]  
            print(f'### {lc_name} ###')
            print("Telescope \tOutburst time \tstd \t exp") 
        else:
            lc_telescope = row[3].strip(' ')
            lc_avg = float(row[4].replace(',', '.'))
            lc_std = float(row[5].replace(',', '.'))         
            lc = Lightcurve(lc_name, lc_telescope)  

            if lc_telescope not in plotted_telescopes:
                plots.lc_modes(lc, styling=True, show=False, save=True)
                plotted_telescopes.append(lc_telescope)

        # if read_obs:
            # Load lightcurve
            # lc_telescope = row[3]
            # lc_avg = float(row[4].replace(',', '.'))
            # lc_std = float(row[5].replace(',', '.'))         
            # lc = Lightcurve(lc_name, lc_telescope)  

            # if lc_telescope not in plotted_telescopes:
            #     plots.lc_modes(lc, styling=True, show=False, save=True)
            #     plotted_telescopes.append(lc_telescope)

            # # Read mean of outburst
            # ob_peaktime = int(row[6])
            # ob = lc.get_outburst(ob_peaktime, std=lc_std)

            # # Remove background
            # ob.correct_background(lc_avg)
            # lc.correct_background(lc_avg)

            # save = True
            # show = False

            # # Gaussian fit
            # gaussian_fit, fit_stddev, fit_mean, fit_dur, fit_y = ob.fit_gaussian()
            # plots.lc_gaussian(lc, ob, fit_dur, gaussian_fit, fit_stddev, fit_y, save=save, show=show, styling=True)

            # # Exponential fit
            # exponential_fit, fit_tdecay, ob_decay, fit_y = ob.fit_decay()
            # plots.lc_exponential(lc, ob, ob_decay, exponential_fit, fit_tdecay, fit_y, save=save, ylog=False, show=show, styling=True)

            # print(f'{ob.telescope}\t\t{ob_decay.ts.time.mjd[0]:.0f}\t{fit_stddev*6:.2f}\t{fit_tdecay:.2f}')

            