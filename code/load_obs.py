import matplotlib.pyplot as plt

def load_obs(file):

    obs = {}
    obs['t_dur'] = [[]]
    obs['t_dec'] = [[]]   
    obs['labels'] = ['1 or 2 outbursts']

    # Open outburst file
    read_file = open(file,'r')

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

            # Initiate space for outburst properties per source
            if read_source:  
                lc_name = row[0] 
                lc_type = row[1]
                lc_num_outbursts = int(row[2])         

                if lc_num_outbursts < 2:
                    i = 0
                else:
                    obs['t_dur'].append([])
                    obs['t_dec'].append([])
                    obs['labels'].append(lc_name)
                    i = -1

            # Save outburst properties
            if read_obs:
                # Save data
                if row[7] != '':
                    # print(f"=dur: {row[7].strip(' ').replace(',', '.')}=")
                    obs['t_dur'][i].append(float(row[7].replace(',', '.')))
                if row[8] != '':
                    obs['t_dec'][i].append(float(row[8].replace(',', '.')))
                    # print(f"=dec: {row[8].strip(' ').replace(',', '.')}=")

    # Close file
    read_file.close()

    return obs

def load_obs_type(file):

    obs = {}
    obs['t_dur'] = [[], [], []]
    obs['t_dec'] = [[], [], []]   
    obs['labels'] = ['NS', 'BH', '?']

    types = {
        'NS': 0,
        'BH': 1,
        '?': 2,
    }

    # Open outburst file
    read_file = open(file,'r')

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
                lc_name = row[0] 
                lc_type = row[1]
                i = types[lc_type]

            if read_obs:
                # Save data
                if row[7] != '':
                    obs['t_dur'][i].append(float(row[7].replace(',', '.')))
                if row[8] != '':
                    obs['t_dec'][i].append(float(row[8].replace(',', '.')))
                    
    # Close file
    read_file.close()

    return obs
    