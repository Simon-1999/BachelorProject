# analysis.py
#
# Simon van Eeden
#
# Running this code will analyse outburst properties

import plots
from load_obs import *
import numpy as np

# Load outbursts file
obs_file = 'data/outbursts/Outbursts.csv'
obs = load_obs(obs_file)

# Duration distrubution
print('----- Duration -----')
print(f"Average = {np.average(obs['duration'])}")
print(f"Median = {np.median(obs['duration'])}")
plots.hist_duration(obs['duration'], bins=20)

# Decay distribution
print('----- Decay timescale -----')
print(f"Average = {np.average(obs['exponent_tau'])}")
print(f"Median = {np.median(obs['exponent_tau'])}")
plots.hist_decay(obs['exponent_tau'], bins=11)


