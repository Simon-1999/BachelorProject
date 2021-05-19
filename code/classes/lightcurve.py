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
            return self.load_rxte_pca(file)

        elif read_type == 'RXTE_asm':
            return self.load_rxte_asm(file)
        
        elif read_type == "Swift":
            return self.load_swift(file)

        elif read_type == "Swiftsec":
            return self.load_swift_sec(file, 54368.68333046)

        elif read_type == "Swiftsec_2":
            return self.load_swift_sec(file, 55588.34722681)

        elif read_type == 'SwiftGC':
            return self.load_swiftgc(file)

        else:
            print('[Error] Telescope name unkown, could not load...')


    def load_rxte_asm(self, file):
        print(f'load_rxte_asm: loading rxte file...')

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
            if l_strip[0] != "%":

                # Make list of columns
                l_split = l_strip.split()

                # Temporary save
                times.append(l_split[0])
                data['time_err_pos'].append(0)
                data['time_err_neg'].append(0)
                data['rate'].append(float(l_split[1]))
                data['rate_err_pos'].append(abs(float(l_split[2])))
                data['rate_err_neg'].append(abs(float(l_split[2])))
                data['mode'].append('ASM')    

        # Close file
        read_file.close()

        # Make timeseries object
        ts = TimeSeries(time=Time(times, format='mjd'),
                        data=data)
        return ts


    def load_rxte_pca(self, file):
        print(f'load_rxte_pca: loading rxte file...')

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
                data['rate_err_pos'].append(abs(float(l_split[2])))
                data['rate_err_neg'].append(abs(float(l_split[2])))
                data['mode'].append('PCA')    

        # Close file
        read_file.close()

        # Make timeseries object
        ts = TimeSeries(time=Time(times, format='mjd'),
                        data=data)
        return ts


    def load_swift(self, file):
        print(f'load_swift: loading Swift file...')

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


    def load_swift_sec(self, file, mjd_offset):
        print(f'load_swift_sec: loading Swift file...')

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
                times.append(str(mjd_offset + float(l_split[0])/60/60/24))
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


    def load_swiftgc(self, file):
        print(f'load_swiftgc: loading Swift GC file...')

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

                if float(l_split[13]) > 0:
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


    def get_outbursts(self, std):

        lc = copy.deepcopy(self)
        lc.name = lc.name + '_outbursts'

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

        # Loop over file
        for i in range(1, len(lc.ts) - 2):

            if float(lc.ts[i]['rate']) > std and float(lc.ts[i + 1]['rate']) > std:

                # Temporary save
                times.append(lc.ts[i]['time'].mjd)
                data['time_err_pos'].append(float(lc.ts[i]['time_err_pos']))
                data['time_err_neg'].append(float(lc.ts[i]['time_err_neg']))
                data['rate'].append(float(lc.ts[i]['rate']))
                data['rate_err_pos'].append(float(lc.ts[i]['rate_err_pos']))
                data['rate_err_neg'].append(float(lc.ts[i]['rate_err_neg']))
                data['mode'].append(lc.ts[i]['mode'])  

                if float(lc.ts[i + 1]['rate']) > std and lc.ts[i + 1]['time'].mjd not in times:

                    # Temporary save
                    times.append(lc.ts[i + 1]['time'].mjd)
                    data['time_err_pos'].append(float(lc.ts[i + 1]['time_err_pos']))
                    data['time_err_neg'].append(float(lc.ts[i + 1]['time_err_neg']))
                    data['rate'].append(float(lc.ts[i + 1]['rate']))
                    data['rate_err_pos'].append(float(lc.ts[i + 1]['rate_err_pos']))
                    data['rate_err_neg'].append(float(lc.ts[i + 1]['rate_err_neg']))
                    data['mode'].append(lc.ts[i + 1]['mode'])  

                # One datapoint after
                if lc.ts[i + 2]['time'].mjd not in times:

                    # Temporary save
                    times.append(lc.ts[i + 2]['time'].mjd)
                    data['time_err_pos'].append(float(lc.ts[i + 2]['time_err_pos']))
                    data['time_err_neg'].append(float(lc.ts[i + 2]['time_err_neg']))
                    data['rate'].append(float(lc.ts[i + 2]['rate']))
                    data['rate_err_pos'].append(float(lc.ts[i + 2]['rate_err_pos']))
                    data['rate_err_neg'].append(float(lc.ts[i + 2]['rate_err_neg']))
                    data['mode'].append(lc.ts[i + 2]['mode']) 

                # One datapoint before 
                if lc.ts[i - 1]['time'].mjd not in times:

                    # Temporary save
                    times.append(lc.ts[i - 1]['time'].mjd)
                    data['time_err_pos'].append(float(lc.ts[i - 1]['time_err_pos']))
                    data['time_err_neg'].append(float(lc.ts[i - 1]['time_err_neg']))
                    data['rate'].append(float(lc.ts[i - 1]['rate']))
                    data['rate_err_pos'].append(float(lc.ts[i - 1]['rate_err_pos']))
                    data['rate_err_neg'].append(float(lc.ts[i - 1]['rate_err_neg']))
                    data['mode'].append(lc.ts[i - 1]['mode']) 

        # Make timeseries object
        ts = TimeSeries(time=Time(times, format='mjd'),
                        data=data)

        lc.ts = ts

        return lc


    def get_outburst(self, peak_time, std):

        lc = copy.deepcopy(self)

        # Split into left and right from peak time
        ob_l = lc.ts.loc[Time('50000', format='mjd'):Time(peak_time, format='mjd')]
        ob_r = lc.ts.loc[Time(peak_time, format='mjd'):Time('59000', format='mjd')]

        # Get outburst borders
        ob_l.sort('time', reverse=True)
        ob_start_time = self.get_border_time(ob_l, std) - 0.01
        ob_end_time = self.get_border_time(ob_r, std) + 0.01

        # Select ob region
        lc.ts = self.ts.loc[Time(ob_start_time, format='mjd'):Time(ob_end_time, format='mjd')]
        
        ##### Plot outburst #####

        # Style
        plt.style.use('default')
        plt.rc('axes', unicode_minus=False)
        mpl.rcParams['font.family'] = "Tw Cen MT"
        mpl.rcParams['patch.linewidth'] = .8
        plt.rcParams['xtick.major.size'] = 5.0
        plt.rcParams['xtick.minor.size'] = 3.0
        plt.rcParams['xtick.top'] = True
        plt.rcParams['ytick.major.size'] = 5.0
        plt.rcParams['ytick.minor.size'] = 3.0
        plt.rcParams['ytick.right'] = True
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.minorticks_on()

        # Title
        plt.title(f'{self.name} - {self.telescope}')

        # Full lichtcurve
        xerr = self.ts['time_err_pos'], self.ts['time_err_pos']
        yerr = self.ts['rate_err_neg'], self.ts['rate_err_neg']
        plt.errorbar(self.ts.time.mjd, self.ts['rate'], xerr=xerr, yerr=yerr, color='gray', fmt='s', ms=3, elinewidth=.5, label='Background')  

        # Outburst region
        xerr = lc.ts['time_err_pos'], lc.ts['time_err_pos']
        yerr = lc.ts['rate_err_neg'], lc.ts['rate_err_neg']
        plt.errorbar(lc.ts.time.mjd, lc.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt='s', ms=3, elinewidth=.5, label='Outburst region')  

        # Std level
        plt.plot([self.ts.time.mjd[0], self.ts.time.mjd[-1]], [std, std], 'k--', label='1$\sigma$')

        # Set xlimit to outburst region
        plt.xlim(ob_start_time - 25, ob_end_time + 25)

        # Legend
        plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)

        plt.show()

        return lc


    def get_border_time(self, ts, std):

        SIGNAL_TRESSHOLD = 2

        i = 0
        for row in ts:

            below_std = row['rate'] < std

            # Track adjecent datapoints below std
            if below_std:
                i += 1
            else:
                i = 0

            # End of outburst 
            if i == SIGNAL_TRESSHOLD:
                return row['time'].mjd        


    def calc_background(self, start_time, end_time, save=False):

        ts = self.ts.loc[Time(start_time, format='mjd'):Time(end_time, format='mjd')]

        # Get std
        std = np.std(ts['rate'])
        ave = np.average(ts['rate'])
        med = np.median(ts['rate'])

        # Style
        plt.style.use('default')
        plt.rc('axes', unicode_minus=False)
        mpl.rcParams['font.family'] = "Tw Cen MT"
        mpl.rcParams['patch.linewidth'] = .8
        plt.rcParams['xtick.major.size'] = 5.0
        plt.rcParams['xtick.minor.size'] = 3.0
        plt.rcParams['xtick.top'] = True
        plt.rcParams['ytick.major.size'] = 5.0
        plt.rcParams['ytick.minor.size'] = 3.0
        plt.rcParams['ytick.right'] = True
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        # plt.minorticks_on()

        # Plot lightcurve
        self.plot_lc()

        # Borders of extraction region
        bge_x_lower = ts['time'].mjd[0]
        bge_x_upper = ts['time'].mjd[-1]
        bge_y_lower = min(ts['rate']) - 1.1 * max(ts['rate_err_neg'])
        bge_y_upper = max(ts['rate']) + 1.1 * max(ts['rate_err_pos'])

        plt.plot([bge_x_lower, bge_x_lower, bge_x_upper, bge_x_upper, bge_x_lower],
            [bge_y_lower, bge_y_upper, bge_y_upper, bge_y_lower, bge_y_lower], 
            'k--', 
            lw=1, 
            label="Extraction region")

        # Title
        plt.title('Non outburst level extraction')

        # Properties
        plt.plot([], [], ' ', label=f"Std = {std:.5f}")
        plt.plot([], [], ' ', label=f"Average = {ave:.5f}")
        plt.plot([], [], ' ', label=f"Median = {med:.5f}")

        # Labels
        plt.xlabel('Time (MJD)')
        plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')
        
        # Legend
        plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)
        if save:
            plt.savefig(f'output/analysis/non outburst level/{self.name}_{self.telescope}.png', dpi=150)
        
        plt.show()


    def correct_background(self, background_level):

        self.ts['rate'] = self.ts['rate'] - background_level

        print(f'correct_background: corrected for {background_level}')


    def plot_lc(self):

        # Query errors
        xerr = self.ts['time_err_pos'], self.ts['time_err_pos']
        yerr = self.ts['rate_err_neg'], self.ts['rate_err_neg']

        # Plot errorbars
        plt.errorbar(self.ts.time.mjd, self.ts['rate'], xerr=xerr, yerr=yerr, color='k', fmt='s', ms=3, elinewidth=.5)  


    def fit_gaussian(self, mean, amplitude_fixed=False, save=False):
        print(f'gaussian_fit: fitting {self.name}...')

        lc = copy.deepcopy(self)

        # Query x and y
        x = self.ts.time.mjd
        y = self.ts['rate']

        amplitude = max(self.ts['rate'])

        # Fit gaussian
        g_init = models.Gaussian1D(
            amplitude=amplitude, 
            mean=mean, 
            stddev=6.)

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

        # Style
        plt.style.use('default')
        plt.rc('axes', unicode_minus=False)
        mpl.rcParams['font.family'] = "Tw Cen MT"
        mpl.rcParams['patch.linewidth'] = .8
        plt.rcParams['xtick.major.size'] = 5.0
        plt.rcParams['xtick.minor.size'] = 3.0
        plt.rcParams['xtick.top'] = True
        plt.rcParams['ytick.major.size'] = 5.0
        plt.rcParams['ytick.minor.size'] = 3.0
        plt.rcParams['ytick.right'] = True
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.minorticks_on()

        # Plot lightcurve
        self.plot_lc()

        # Plot gauss
        x_gauss = np.arange(x[0]-10, x[-1]+10, 0.1)
        y_gauss = g(x_gauss)
        # plt.plot(x_gauss, y_gauss, 'r-', lw=1, label='Gaussian fit')

        data = {        
            'rate': y_gauss,
        }

        # Make timeseries object
        ts = TimeSeries(time=Time(x_gauss, format='mjd'),
                        data=data)

        lc.ts = ts
        lc.name = 'Gaussian fit'
        return lc

        # # Plot 1 std deviation

        # # Plot 10 percent line
        # # plt.plot(x, [max(y)*0.1]*len(x), 'k-', label='10%')
        # # plt.plot(x, [max(y)*0.2]*len(x), 'k-', label='20%')
        # # plt.plot(x, [max(y)*0.01]*len(x), 'k-', label='1%')

        # # Plot 6*sigma area below guass
        # lower_x = g.mean.value - 3*g.stddev.value
        # upper_x = g.mean.value + 3*g.stddev.value
        # # x = np.arange(lower_x, upper_x, 0.1)
        # plt.plot([lower_x]*2, [-10, max(y)+10], '--', lw=1, color='r', label="6$\sigma$")
        # plt.plot([upper_x]*2, [-10, max(y)+10], '--', lw=1, color='r')
        # # plt.plot(x, g(x), '--', lw=1, color = 'k', label="6$\sigma$")

        # # Find left datapoints within range
        # ts_range_left_rate = []
        # ts_range_left_time = []

        # self.ts.sort('time')
        # ts_left = self.ts.loc[Time(f'{x[0]}', format='mjd'):Time(f'{g.mean.value}', format='mjd')]
        # ts_left.sort('time', reverse=True)
        # F_10_rate_left = 0
        # F_10_time_left = 0
        # F_10_dif_left = max(y)
        # eor = False
        # first_low = False
        # for row in ts_left:
        #     if first_low and row['rate'] < 0.05 * max(y):
        #         eor = True
        #         print('done')
        #         break

        #     if row['rate'] < 0.05 * max(y):
        #         first_low = True
        #     else:
        #         first_low = False

        #     if row['rate'] > 0.05 * max(y) and row['rate'] < 0.2 * max(y) and not eor:
        #         dif_10 = abs(0.10 * max(y) - row['rate'])

        #         # Save closest to 10% value
        #         if dif_10 < F_10_dif_left:
        #             F_10_rate_left = row['rate']
        #             F_10_time_left = row['time'].mjd
        #             F_10_dif_left = dif_10

        #         ts_range_left_rate.append(row['rate'])
        #         ts_range_left_time.append(row['time'].mjd)
                
        # plt.plot(ts_range_left_time, ts_range_left_rate, 's', color='b', fillstyle='none')
        # # plt.plot(F_10_time_left, F_10_rate_left, 'ko', ms=10)
        # plt.plot([F_10_time_left]*2, [-10, max(y)+10], '--', lw=1, color='b', label='F_10%')

        # # Find right datapoints within range
        # ts_range_right_rate = []
        # ts_range_right_time = []

        # # self.ts.sort('time', reverse=True)
        # ts_right = self.ts.loc[Time(f'{g.mean.value}', format='mjd'):Time(f'{x[-1]}', format='mjd')]
        # # ts_right.sort('time', reverse=True)
        # F_10_rate_right = 0
        # F_10_time_right = 0
        # F_10_dif_right = max(y)
        # eor = False
        # first_low = False
        # for row in ts_right:
        #     if row['rate'] < 0.2 * max(y):
        #         below_upper = True
        #     else:
        #         below_upper = False

        #     if row['rate'] > 0.05 * max(y):
        #         above_lower = True
        #     else:
        #         above_lower = False

        #     if first_low and not above_lower:
        #         eor = True
        #         print('done')
        #         break

        #     if not above_lower:
        #         first_low = True
        #     else:
        #         first_low = False

        #     if above_lower and below_upper and not eor:
        #         dif_10 = abs(0.10 * max(y) - row['rate'])

        #         # Save closest to 10% value
        #         if dif_10 < F_10_dif_right:
        #             F_10_rate_right = row['rate']
        #             F_10_time_right = row['time'].mjd
        #             F_10_dif_right = dif_10
        #         ts_range_right_rate.append(row['rate'])
        #         ts_range_right_time.append(row['time'].mjd)
                
        # plt.plot(ts_range_right_time, ts_range_right_rate, 's', color='b', fillstyle='none', label="Within 20-1%")
        # # plt.plot(F_10_time_right, F_10_rate_right, 'ko', ms=10)
        # plt.plot([F_10_time_right]*2, [-10, max(y)+10],'--', lw=1, color='b')

        # # props = dict(boxstyle='square', facecolor='white', alpha=0)
        # # textstr = f'Amplitude = {g.amplitude.value:.2f}\nMean = {g.mean.value:.2f}\nStddev = {g.stddev.value:.2f}'
        # # plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
        # #     verticalalignment='top', bbox=props)

        # # Plot fit parameters
        # # plt.plot([], [], ' ', label=f"Stddev = {g.stddev.value:.2f}")
        # # plt.plot([], [], ' ', label=f"Mean = {g.mean.value:.2f}")
        # # plt.plot([], [], ' ', label=f"Amplitude = {g.amplitude.value:.2f}")
        # plt.plot([], [], ' ', label=f"Duration Gauss = {(g.stddev.value*6):.2f}")
        # plt.plot([], [], ' ', label=f"Duration F_10% = {(F_10_time_right - F_10_time_left):.2f}")

        # # Title
        # plt.title(f'{self.name} - {self.telescope}')

        # # Labels
        # plt.xlabel('Time (MJD)')
        # plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')
        
        # # Legend
        # plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)
        
        # if save:
        #     plt.savefig(f'output/analysis/gaussian fits/{self.name}_{self.telescope}_{g.mean.value:.0f}.png', dpi=150)
        # plt.show()

        # return 


    def fit_exp(self, amplitude=5, tau=-1, save=False):
        print(f'fit_exp: fitting {self.name}...')

        # Query x and y
        x = self.ts.time.mjd
        t_start = float(x[0])
        for i, t in enumerate(x):
            x[i] = str(float(x[i]) - t_start)
        y = self.ts['rate']

        # Fit exponent
        g_init = models.Exponential1D(amplitude=amplitude, tau=tau)
        fit = fitting.LevMarLSQFitter()
        g = fit(g_init, x, y, weights=1/self.ts['rate_err_pos'])

        # Print parameters
        print("--- Fit parameters ---")
        print(f'Amplitude = {g.amplitude.value}')
        print(f'Tau = {g.tau.value}')
        print("----------------------")

        # Plot the data with the best-fit model
        print('fit_exp: plotting...')
        plt.style.use('default')
        plt.rc('axes', unicode_minus=False)
        mpl.rcParams['font.family'] = "Tw Cen MT"
        mpl.rcParams['patch.linewidth'] = .8
        plt.rcParams['xtick.major.size'] = 5.0
        plt.rcParams['xtick.minor.size'] = 3.0
        plt.rcParams['xtick.top'] = True
        plt.rcParams['ytick.major.size'] = 5.0
        plt.rcParams['ytick.minor.size'] = 3.0
        plt.rcParams['ytick.right'] = True
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        plt.minorticks_on()
        self.plot_lc()
        x = np.arange(x[0], x[-1]+1, 0.1)
        plt.plot(x, g(x), 'k--', label='Exponential fit')
        props = dict(boxstyle='square', facecolor='white', alpha=0)
        plt.plot([], [], ' ', label=f"Tau = {g.tau.value:.2f}")
        plt.plot([], [], ' ', label=f"Amplitude = {g.amplitude.value:.2f}")
        # textstr = f'Amplitude = {g.amplitude.value:.2f}\nTau = {g.tau.value:.2f}'
        # plt.text(0.05, 0.05, textstr, transform=plt.gca().transAxes,
        #     verticalalignment='bottom', bbox=props)
        plt.title(f'{self.name} - {self.telescope}')
        plt.yscale('log')
        plt.xlabel('Time (days)')
        plt.ylabel('Rate ($\mathregular{c}$ $\mathregular{s}^{\mathregular{-1}}$)')
        plt.legend(shadow=False, edgecolor='k', fancybox=False, borderaxespad=1)
        # lg.get_frame().set_boxstyle('square', pad=0.0)
        if save:
            plt.savefig(f'output/analysis/exponential fits/{self.name}_{self.telescope}_{t_start:.0f}.png', dpi=150)
        plt.show()

        return 

    

