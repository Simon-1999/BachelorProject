# lightcurve.py
# 
# Simon van Eeden
#
# This file contains the lightcurve class

from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy.timeseries import aggregate_downsample
from astropy.modeling import models, fitting
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib as mpl
import copy
import json
import numpy as np

class Lightcurve():

    def __init__(self, name, telescope):
 
        self.name = name
        self.telescope = telescope
        self.ts = self.load_ts()
        self.ts.sort('time')

    def load_ts(self):

        # Search for filepath
        with open(f'{self.telescope}_paths.json') as json_file:
            paths = json.load(json_file)
        file = paths[self.name]['path']
        read_type = paths[self.name]['read_type']

        # Select reading method
        if read_type == 'RXTE':
            return self.load_rxte_ts(file)
        
        elif read_type == "Swift":
            return self.load_swift_ts(file)

        elif read_type == "Swiftsec":
            return self.load_swiftsec_ts(file)

        elif read_type == 'SwiftGC':
            return self.load_swiftgc_ts(file)

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
            if l_strip[0] != ";" and l_strip[0] != "%":

                # Make list of columns
                l_split = l_strip.split()

                # Temporary save
                times.append(l_split[0])
                data['time_err_pos'].append(0)
                data['time_err_neg'].append(0)
                data['rate'].append(float(l_split[1]))
                data['rate_err_pos'].append(abs(float(l_split[2])))
                data['rate_err_neg'].append(abs(float(l_split[2])))
                data['mode'].append('RXTE')    

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


    def load_swiftsec_ts(self, file):
        print(f'Loading file {file} ...')

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
                times.append(str(54368.68333046 + float(l_split[0])/60/60/24))
                data['time_err_pos'].append(float(l_split[1])/60/60/24)
                data['time_err_neg'].append(abs(float(l_split[2])/60/60/24))
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
                data['mode'].append('SwiftGC')    

        # Close file
        read_file.close()

        # Make timeseries object
        ts = TimeSeries(time=Time(times, format='mjd'),
                        data=data)
        return ts


    def get_rate_lower(self):

        rate_lower = []

        for i in range(len(self.ts)):

            rate_lower.append(self.ts['rate'][i] - self.ts['rate_err_neg'][i])

        return rate_lower


    def get_rate_upper(self):

        rate_upper = []

        for i in range(len(self.ts)):

            rate_upper.append(self.ts['rate'][i] + self.ts['rate_err_pos'][i])

        return rate_upper


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


    def get_fraction(self, start_time, end_time):

        lc = copy.deepcopy(self)

        # Select fraction
        lc.ts = lc.ts.loc[Time(start_time, format='mjd'):Time(end_time, format='mjd')]

        return lc


    def plot_lc(self):

        # Query errors
        xerr = self.ts['time_err_pos'], self.ts['time_err_pos']
        yerr = self.ts['rate_err_neg'], self.ts['rate_err_neg']

        # Plot errorbars
        plt.errorbar(self.ts.time.mjd, self.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt='s', ms=3, elinewidth=.5)  


    def fit_gaussian(self, amplitude, mean, amplitude_fixed=False, save=True):
        print(f'gaussian_fit: fitting {self.name}...')

        # Query x and y
        x = self.ts.time.mjd
        y = self.ts['rate']

        # Fit gaussian
        g_init = models.Gaussian1D(
            amplitude=amplitude, 
            mean=mean, 
            stddev=3.)
        g_init.amplitude.fixed = amplitude_fixed
        # g_init.stddev.fixed=True
        fit = fitting.LevMarLSQFitter()
        g = fit(g_init, x, y, weights=1/self.ts['rate_err_pos'])

        # Print parameters
        print("--- Fit parameters ---")
        print(f'Amplitude = {g.amplitude.value}')
        print(f'Mean = {g.mean.value}')
        print(f'Stddev = {g.stddev.value}')
        print("----------------------")

        # Plot the data with the best-fit model
        print('gaussian_fit: plotting...')
        self.plot_lc()
        x = np.arange(x[0]-10, x[-1]+10, 0.1)
        plt.plot(x, g(x), label='Gaussian fit')
        props = dict(boxstyle='square', facecolor='white', alpha=0)
        textstr = f'Amplitude = {g.amplitude.value:.2f}\nMean = {g.mean.value:.2f}\nStddev = {g.stddev.value:.2f}'
        plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
            verticalalignment='top', bbox=props)
        plt.title(f'{self.name} - {self.telescope}')
        plt.xlabel('Time (MJD)')
        plt.ylabel('Rate ($c$ $s^{-1}$)')
        plt.legend(shadow=False, edgecolor='k')
        if save:
            plt.savefig(f'output/analysis/gaussian fits/{self.name}_{self.telescope}_{g.mean.value:.0f}.png', dpi=150)
        plt.show()

        return 

