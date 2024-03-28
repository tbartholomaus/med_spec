#!/opt/anaconda/envs/seisenv/bin python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Initial version created on Fri Nov 10 19:35:50 2017 by @author: Tim Bartholomaus. 

Most recent update by Yoram Terleth, February 2024. 
General clean-up, and improved user-friendliness. Added saving of the produced GHT timeseries to pickle and mat files.
"""

# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
from scipy.ndimage import median_filter
import statsmodels.api as sm
import pickle
import os, sys, time
from scipy.io import savemat
import pandas as pd 

# Set environment time to UTC
os.environ['TZ'] = 'UTC'
time.tzset()

# USER INPUT ############################################################################
# Path to pickle files
path_to_pickle_files = '/data/stor/proj/Turner/med_spec_output/complete_spec_pickles/'

# Stations to analyze
stations = ['YG_SG14']

# dict variable
BB = dict() # BB is a dict to contain multiple versions of the class BBseis

# Frequency range for integration
fGHT = [3,10]
#######################################################################################

# Class for storing seismic data
class BBseis:
    def __init__(self, name, t, GHT_powdB):
        self.name = name
        self.t = t
        self.powdB = GHT_powdB
        self.pdB_LF = np.array([0])
        self.T_amp = np.array([0])
        self.Ta_LF = np.array([0])

# Function for smoothing tremor data
def smooth_tremor(filt_a, filt_b, station, t, GHT_powdB):
    GHT_pow = 10**(GHT_powdB/10)
    T_amp = np.sqrt(GHT_pow)

    # Pre-processing
    first_ind = np.where(GHT_powdB > -300)[0][0]
    last_ind = np.where(GHT_powdB > -300)[0][-1]
    GHT_powdB[:first_ind] = np.nan
    GHT_powdB[last_ind+1:] = np.nan
    T_amp[:first_ind] = np.nan
    T_amp[last_ind+1:] = np.nan

    first = GHT_powdB[first_ind+5]
    last  = GHT_powdB[last_ind-5 ]
    first_ind = first_ind +1
    last_ind  = last_ind  -1
    padded = GHT_powdB.copy()
    padded[:first_ind] = first
    padded[last_ind+1:] = last

    # Low-pass filtering
    pdB_LF = signal.filtfilt(filt_b, filt_a, padded)
    pdB_LF[:first_ind] = np.nan
    pdB_LF[last_ind+1:-1] = np.nan

    padded = T_amp.copy()
    first = padded[first_ind+5]
    last  = padded[last_ind-5 ]
    padded[:first_ind] = first
    padded[last_ind+1:] = last
    
    IQR = np.nanpercentile(T_amp, 75) - np.nanpercentile(T_amp, 25)
    pctile_thresh = np.nanmedian(T_amp) + 10*IQR
    try:
        bad_ind = np.where(T_amp > pctile_thresh)
    except RuntimeWarning:
        pass
    T_amp[bad_ind] = np.nan
    GHT_powdB[bad_ind] = np.nan
    
    lowess = sm.nonparametric.lowess
    frac = .5 * 86400/( (t[-1] - t[0]).astype('int64') )
    Ta_LF = lowess(T_amp, t, frac=frac, it=5, return_sorted=False)

    Ta_LF[:first_ind] = np.nan
    Ta_LF[last_ind+1:-1] = np.nan
    Ta_LF[np.where(GHT_powdB < -300)[0]] = np.nan
    GHT_powdB[np.where(GHT_powdB < -300)[0]] = np.nan
    
    return GHT_powdB, pdB_LF, T_amp, Ta_LF

# Iterate over each station
for station in stations:
    # Load seismic data from pickle file
    with open(path_to_pickle_files + 'mp' + station + '.pickle', 'rb') as f:
        t, t_dt64, freqs, Pdb_array, pp, data_dir, station, date = pickle.load(f, encoding='latin1')

    # Calculate start and end times
    t_start = np.min(t_dt64)
    t_end = np.max(t_dt64)
    
    # Calculate Nyquist frequency and low-pass filter parameters
    sample_period = pp['coarse_duration'] * pp['coarse_overlap'] / 86400
    Nyq = 1/sample_period/2
    lp_cutoff_period = 6
    lp_cutoff_freq = 1/(lp_cutoff_period/24)
    filt_b, filt_a = signal.butter(6, lp_cutoff_freq/Nyq, 'low')

    # Convert power to linear scale
    Pow = 10**(Pdb_array/10)
    ind = np.where(np.all([freqs>fGHT[0], freqs<fGHT[1]], axis=0))[0] 
    GHT_freqs = freqs[ind]

    # Integrate power over frequency range
    GHT_pow = np.sum(Pow[ind,:], axis=0) * np.unique(np.diff(freqs))
    
    # Store seismic data in BBseis object
    BB[ station ] = BBseis(station, t_dt64, 10*np.log10(GHT_pow))
    
    # Smooth tremor data
    BB[ station ].GHT_powdB, BB[ station ].pdB_LF, BB[ station ].T_amp, BB[ station ].Ta_LF = smooth_tremor(
            filt_a, filt_b, BB[station].name, BB[station].t, BB[station].powdB)

# Interpolate data to regular time intervals
t_start = np.datetime64(t_start, 'D')
t_end = np.datetime64(t_end, 'D') + 1
t_interp = np.arange(t_start, t_end, 
                     np.timedelta64(15, 'm'), dtype='datetime64')

# Plotting
fig, ax = plt.subplots(2,1, sharex=True)
for station in BB:
    interp_func = interp1d(BB[station].t.astype('datetime64[m]').astype('int64'), 
                           BB[station].Ta_LF, fill_value='extrapolate', 
                           kind='linear')
    BB[ station ].Ta_int = interp_func(t_interp.astype('int64'))

    ax[1].plot(BB[station].t, pd.DataFrame(BB[station].powdB).rolling(1).mean(), label=station, linewidth=0.3)
        
    ax[0].plot(BB[station].t, BB[station].T_amp, linewidth=0.3)

ax[1].set_ylabel( 'Power (dB rel. 1 m^2/s^2)' )
ax[0].set_ylabel('Tremor amp (m/s)')




# SAVING: 

# save the figure: 
fig_path = '/data/stor/proj/Turner/med_spec_2024/output_figs/'

# Check if the directory exists, if not, create it
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

plt.savefig(fig_path+'GHT_'+ station + '', dpi=600)


# save the GHT timereries as a pickle file 
GHT_OUT_path = '/data/stor/proj/Turner/med_spec_2024/output_GHT/'

# Check if the directory exists, if not, create it
if not os.path.exists(GHT_OUT_path):
    os.makedirs(GHT_OUT_path)

with open(GHT_OUT_path + str(fGHT) + 'all.pickle', 'wb') as f:
    pickle.dump([BB, fGHT, t_interp], f)


# save the timeseries to a MAT file: same path as pickle GHT file. 
savemat(GHT_OUT_path+'mat_file'+ station +str(fGHT)+'.mat',{'a':[BB, fGHT, t_interp]},appendmat=True)
