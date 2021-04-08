# swift.py
# 
# Simon van Eeden
#
# This file contains functions for processing Swift light curve data

def read_swift(file, name):
    """Opens swift file and returns dictionary with data"""

    data = {
        'telescope': 'swift',
        'timestamps': [],
        'name': name,
    }

    read = False
    mode_type = 'Unkown'

    # Open read and write file
    original_file= open(file,'r')    

    for l in original_file:
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

            # Convert into dictionary
            timestamp = {
                'time': float(time),
                'time_error_up': float(time_error_up),
                'time_error_down': abs(float(time_error_down)),
                'rate': float(rate),
                'rate_error_up': float(rate_error_up),
                'rate_error_down': abs(float(rate_error_down)),
                'mode': mode_type,
            }

            # Save timestamp
            data['timestamps'].append(timestamp)

        # Query 
        if l_strip.startswith('! WT data'):
            read = True   
            mode_type = 'WT'

        if l_strip.startswith('! PC data'):
            read = True
            mode_type = 'PC'

        if l_strip.startswith('! PC Upper limit'):
            read = False

    original_file.close()

    return data
           