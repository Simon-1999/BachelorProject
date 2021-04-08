# RXTE.py
# 
# Simon van Eeden
#
# This file contains functions for processing RXTE light curve data

from astropy.time import Time
from astropy.timeseries import TimeSeries

def read_RXTE(file, name):
    """Opens RXTE file and returns dictionary with data"""

    data = {
        'telescope': 'rxte',
        'timestamps': [],
        'name': name,
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

            # Query column values
            time = l_split[0]
            rate = l_split[1]

            # Convert into dictionary
            timestamp = {
                'time': float(time),
                'time_error_up': 0,
                'time_error_down': 0,
                'rate': float(rate),
                'rate_error_up': 0,
                'rate_error_down': 0,
                'mode': 'std',
            }

            # Save timestamp
            data['timestamps'].append(timestamp)           

    # Close file
    read_file.close()

    return data

def load_RXTE(file, name):

    times = []
    data = {        
        'time_err_pos': [],
        'time_err_neg': [],
        'rate': [],
        'rate_err_pos': [],
        'rate_err_neg': [],
        'mode': [],
    }
    time_list, time_err_pos, time_err_neg = [], [], []
    rate, rate_err_pos, rate_err_neg = [], [], []
    mode = []

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

            # Convert into dictionary
            times.append(l_split[0])
            data['time_err_pos'].append(0)
            data['time_err_neg'].append(0)
            data['rate'].append(float(l_split[1]))
            data['rate_err_pos'].append(0)
            data['rate_err_neg'].append(0)
            data['mode'].append('std')    

    # Close file
    read_file.close()

    ts = TimeSeries(time=Time(times, format='mjd'),
                    data=data)

    print(ts.info)
    return ts

def process_RXTE(file, object_name):
    """Opens RXTE file, processes data and saves it in text file"""

    # Open read and write file
    original_file= open(file,'r')
    processed_file = open(f'data/processed/RXTE/{object_name}.txt','w')

    first_row_passed = False

    for l in original_file:

        l_strip = l.strip()
        print(l_strip)
       
        # Process data
        if l_strip[0] != ";":

            if not first_row_passed:
                processed_file.write(f'#Time (MJD)\tRate\n')              
                first_row_passed = True

            l_split = l_strip.split()

            time = l_split[0]
            rate = l_split[1]

            processed_file.write(f'{time}\t{rate}\n')

    original_file.close()
    processed_file.close()
           