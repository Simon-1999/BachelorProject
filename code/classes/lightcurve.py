# lightcurve.py
# 
# Simon van Eeden
#
# This file contains the lightcurve class

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy.timeseries import aggregate_downsample
from astropy.modeling import models, fitting
from astropy import units as u
import copy

class Lightcurve():

    def __init__(self, name, telescope, file):
 
        self.name = name
        self.telescope = telescope

        if telescope == 'RXTE':
            self.ts = self.load_rxte_ts(file)
        elif telescope == 'Swift':
            self.ts = self.load_swift_ts(file)
        elif telescope == 'Swift GC':
            self.ts = self.load_swiftgc_ts(file)
        else:
            print('[Error] Telescope name unkown, could not load...')


    def load_rxte_ts(self, file):
        print(f'Lightcurve: loading rxte file {file} ...')

        # Initialize temporary storage frame
        times = []
        data = {        
            'time_err_pos': [],
            'time_err_neg': [],
            'rate': [],
            'rate_err_pos': [],
            'rate_err_neg': [],
            'mode': [],
        }

        # Open read file
        read_file= open(file,'r')

        # Loop over file
        for l in read_file:

            # Remove leading or trailing spaces
            l_strip = l.strip()
        
            # Process data
            if l_strip[0] != ";":

                # Make list of columns
                l_split = l_strip.split()

                # Temporary save
                times.append(l_split[0])
                data['time_err_pos'].append(0)
                data['time_err_neg'].append(0)
                data['rate'].append(float(l_split[1]))
                data['rate_err_pos'].append(0)
                data['rate_err_neg'].append(0)
                data['mode'].append('std')    

        # Close file
        read_file.close()

        # Make timeseries object
        ts = TimeSeries(time=Time(times, format='mjd'),
                        data=data)
        return ts


    def load_swift_ts(self, file):
        print(f'Lightcurve: loading swift file {file} ...')

        # Initialize temporary storage frame
        times = []
        data = {        
            'time_err_pos': [],
            'time_err_neg': [],
            'rate': [],
            'rate_err_pos': [],
            'rate_err_neg': [],
            'mode': [],
        }

        # Loop variables
        read = False
        mode_type = 'Unkown'

        # Open read and write file
        read_file= open(file,'r')    

        for l in read_file:
            l_strip = l.strip()
        
            # Process data
            if read and not l_strip.startswith('!'):

                l_split = l_strip.split()

                # Query data
                time = l_split[0]
                time_error_up = l_split[1]
                time_error_down = l_split[2]
                rate = l_split[3]
                rate_error_up = l_split[4]
                rate_error_down = l_split[5]

                # Temporary save
                times.append(l_split[0])
                data['time_err_pos'].append(float(l_split[1]))
                data['time_err_neg'].append(abs(float(l_split[2])))
                data['rate'].append(float(l_split[3]))
                data['rate_err_pos'].append(float(l_split[4]))
                data['rate_err_neg'].append(abs(float(l_split[5])))
                data['mode'].append(mode_type)  

            # Mode type
            if l_strip.startswith('! WT data'):
                read = True   
                mode_type = 'WT'

            if l_strip.startswith('! PC data'):
                read = True
                mode_type = 'PC'

            if l_strip.startswith('! PC Upper limit'):
                read = False

        read_file.close()

        # Make timeseries object
        ts = TimeSeries(time=Time(times, format='mjd'),
                        data=data)
        return ts


    def load_swiftgc_ts(self, file):
        print(f'Lightcurve: loading Swift GC file {file} ...')

        # Initialize temporary storage frame
        times = []
        data = {        
            'time_err_pos': [],
            'time_err_neg': [],
            'rate': [],
            'rate_err_pos': [],
            'rate_err_neg': [],
            'mode': [],
        }

        # Open read file
        read_file= open(file,'r')

        # Loop over file
        for l in read_file:

            # Remove leading or trailing spaces
            l_strip = l.strip()
        
            # Process data
            if l_strip[0] != "#":

                # Make list of columns
                l_split = l_strip.split()

                # Temporary save
                times.append(l_split[2])
                data['time_err_pos'].append(0)
                data['time_err_neg'].append(0)
                data['rate'].append(float(l_split[13]))
                data['rate_err_pos'].append(abs(float(l_split[14])))
                data['rate_err_neg'].append(abs(float(l_split[14])))
                data['mode'].append('std')    

        # Close file
        read_file.close()

        # Make timeseries object
        ts = TimeSeries(time=Time(times, format='mjd'),
                        data=data)
        return ts


    def plot_lc(self):

        CM = {
            'PC': 'r',
            'WT': 'b',
            'std': 'k',
        }

        # Plot time vs rate
        for row in self.ts:  
            plt.plot(row['time'].mjd, row['rate'], 's', markersize=3, color=CM[row['mode']])

        # Add labels
        patches = []
        for mode in CM.keys():
            if mode in self.ts['mode']:
                patches.append(mpatches.Patch(color=CM[mode], label=f'{mode} mode'))
        plt.legend(handles=patches)

        # Query errors
        xerr = self.ts['time_err_pos'], self.ts['time_err_pos']
        yerr = self.ts['rate_err_neg'], self.ts['rate_err_neg']

        # Plot errorbars
        plt.errorbar(self.ts.time.mjd, self.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt=' ')   


    def plot_lc_clean(self):

        if self.ts.__class__.__name__ == 'BinnedTimeSeries':
            plt.plot(self.ts['time_bin_start'].mjd, self.ts['rate'], '-', drawstyle='steps-post', label=self.name)

        if self.ts.__class__.__name__ == 'TimeSeries':
            # Query errors
            xerr = self.ts['time_err_pos'], self.ts['time_err_pos']
            yerr = self.ts['rate_err_neg'], self.ts['rate_err_neg']
            plt.errorbar(self.ts.time.mjd, self.ts['rate'], xerr=xerr, yerr=yerr, fmt='s', ms=3, label=self.name)  


    def binning(self, bin_days=10):

        lc = copy.deepcopy(self)
        
        # Delete mode column
        del lc.ts['mode']

        # Bin data
        ts_binned = aggregate_downsample(lc.ts, time_bin_size=bin_days * u.day) 
        lc.ts = ts_binned

        # Change name
        lc.name = lc.name + f"_binned{bin_days}"

        return lc


    def fraction(self, start_time, end_time):

        lc = copy.deepcopy(self)

        # Select fraction
        lc.ts = lc.ts.loc[start_time:end_time]

        # Change name
        lc.name = lc.name + "_frac"

        return lc


    def gaussian_fit(self, amplitude, mean):
        print(f'gaussian_fit: fitting {self.name}...')

        # Query x and y
        x = self.ts.time.mjd
        y = self.ts['rate']

        # Fit gaussian
        g_init = models.Gaussian1D(amplitude=amplitude, mean=mean, stddev=1.)
        fit = fitting.LevMarLSQFitter()
        g = fit(g_init, x, y)

        # Print parameters
        print("--- Fit parameters ---")
        print(f'Amplitude = {g.amplitude.value}')
        print(f'Mean = {g.mean.value}')
        print(f'Stddev = {g.stddev.value}')
        print("----------------------")

        # Plot the data with the best-fit model
        print('gaussian_fit: plotting...')
        plt.plot(x, y, 'ks')
        plt.plot(x, g(x), label='Gaussian')
        plt.title(f'{self.name} - {self.telescope}')
        plt.xlabel('$Time$ $(MJD)$')
        plt.ylabel('$Rate$')
        plt.legend()
        plt.show()

        return 
                