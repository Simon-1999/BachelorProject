# analysis.py
#
# Simon van Eeden
#
# Running this code will analyse outburst properties

import plots
from load_obs import *
import numpy as np

# Load outbursts file
obs_file = 'data/outbursts/Outbursts_fitparam.tsv'
obs = load_obs_type(obs_file)

# Duration distrubution
plots.hist_duration_stacked(obs['t_dur'], obs['labels'], bins=10, log=True, styling=True, save=False)

# Calculate average and median of duration
durations = []
for lis in obs['t_dur']:
    for dur in lis:
        durations.append(dur)

# print(f"Average = {np.average(durations)}")
# print(f"Median = {np.median(durations)}")

# Decay distribution
plots.hist_decay_stacked(obs['t_dec'], obs['labels'], bins=8, log=True, styling=True, save=False)

# Calculate average and median of decay time
decays = []
for lis in obs['t_dec']:
    for dec in lis:
        decays.append(dec)

# print(f"Average = {np.average(decays)}")
# print(f"Median = {np.median(decays)}")

# Load obs per source
obs = load_obs(obs_file)

# Distribution of duration and decay time per source
plots.hist_decay(obs['t_dec'], obs['labels'], bins=8, log=True, styling=True, save=False)
plots.hist_duration(obs['t_dur'], obs['labels'], bins=10, log=True, styling=True, save=False)


