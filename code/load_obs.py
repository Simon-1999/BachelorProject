def load_obs(file):

    GAUSS_DURATION_RATIO = 6

    obs = {}
    obs['gauss_std'] = []
    obs['duration'] = []
    obs['exponent_tau'] = []    

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

            # Save data
            if l_split[4] != '':
                # print(l_split[4])
                obs['gauss_std'].append(float(l_split[4]))
                obs['duration'].append(float(l_split[4])*GAUSS_DURATION_RATIO)
            if l_split[5] != '':
                obs['exponent_tau'].append(float(l_split[5]))

    # Close file
    read_file.close()

    return obs