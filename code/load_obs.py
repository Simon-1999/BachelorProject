import matplotlib.pyplot as plt

def load_obs(file):

    GAUSS_DURATION_RATIO = 6

    obs = {}
    obs['gauss_std'] = [[]]
    obs['duration'] = [[]]
    obs['exponent_tau'] = [[]]   
    obs['labels'] = ['1 or 2 outbursts']

    types = {
        'NS': 0,
        'BH': 1,
        'LMXB': 2,
        'SFXT': 3,
        '?': 4,
        '': 4,
    }

    # Open outburst file
    read_file = open(file,'r')

    # Loop over file
    for l in read_file:

        # Remove spaces
        l = l.strip()

        # Unread header        
        if not l[0] == '#':

            # Get columns
            l_split = l.split(',')
            read_obs = l_split[0] == ''

            if not read_obs:  
                lc_name = l_split[0] 
                lc_type = l_split[2]
                lc_num_outbursts = int(l_split[3])                

                if lc_num_outbursts < 3:
                    i = 0
                else:
                    obs['gauss_std'].append([])
                    obs['duration'].append([])
                    obs['exponent_tau'].append([])
                    obs['labels'].append(lc_name)
                    i = -1

            if read_obs:
                # Save data
                if l_split[5] != '':
                    obs['gauss_std'][i].append(float(l_split[5]))
                    obs['duration'][i].append(float(l_split[5])*GAUSS_DURATION_RATIO)
                if l_split[6] != '':
                    obs['exponent_tau'][i].append(float(l_split[6]))

    # Close file
    read_file.close()

    return obs