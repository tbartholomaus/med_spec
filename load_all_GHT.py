#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 19:35:50 2017

@author: timb
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
from scipy.ndimage.filters import median_filter

#import matplotlib.mlab as mlab
from matplotlib import dates as mdates

#from obspy.core import read
#from obspy import UTCDateTime
#from obspy import signal

#import datetime as dt
import pickle

#import glob
import os, sys, time
#import fnmatch
#import imp
#imp.reload(get_med_spectra)

#import get_med_spectra
#from med_spectra.UTCDateTime_funcs import UTCfloor, UTCceil, UTC2dn # custom functions for working with obspy UTCDateTimes

os.environ['TZ'] = 'UTC' # Set the time to be UTC
time.tzset()
#%%
class BBseis:
    def __init__(self, name, t, GHT_powdB):
        self.name = name
        self.t = t
        self.powdB = GHT_powdB    # dB rel. 1 (m/s)^2
        self.pdB_LF = np.array([0]) # low frequency filtered version of dB
        self.T_amp = np.array([0]) # m/s : Tremor amplitude
        self.Ta_LF = np.array([0]) # m/s : low frequency filtered version of Tremor amplitude



#%%
stations = ['BBWU', 'BBEU', 'BBGU', 'BBWL', 'BBEL', 'BBGL']
#stations = ['BBWU', 'BBEU', 'BBWL', 'BBEL', 'BBGL']
BB = dict()

fGHT = [1.5, 10]
fGHT = [15, 35]

lp_cutoff_period = 6 # hrs
t_interp = np.arange('2017-06-29T00:00Z', '2017-09-26T00:00Z', 
                     np.timedelta64(15, 'm'), dtype='Datetime64')

for station in stations:
    print(station)
#    station = 'BBWU'#TWLV'


# %% Getting back the objects:
    with open('output_results/mp' + station + '.pickle', 'rb') as f:  # Python 3: open(..., 'rb')
        t, t_dt64, freqs, Pdb_array, pp, data_dir, station = pickle.load(f, encoding='latin1')

    sample_period = pp['coarse_duration'] * pp['coarse_overlap'] / 86400 # days  Rate at which the GHT had been sampled
    Nyq = 1/sample_period/2 # 1/days
    lp_cutoff_freq = 1/(lp_cutoff_period/24) # 1/days
    filt_b, filt_a = signal.butter(6, lp_cutoff_freq/Nyq, 'low') # Construct the filter terms used for lowpass filtering the GHT signal

#sys.exit()

##%%
#    mask_val1 = Pdb_array<=-300
#    mask_val2 = np.isinf(Pdb_array)
#    Pdb_array_mask = np.ma.masked_where(np.any([mask_val1, mask_val2], axis=0), Pdb_array)
#
##xx = np.amin(Pdb_array_mask, axis=1)
###xx = np.median(Pdb_array_mask, axis=1)
###Pdb_array_mask - xx[:,np.newaxis]
##Pdb_array_mask = Pdb_array_mask - xx[:,np.newaxis]
    if station == 'BBGU':
        t         = t[:2244]
        t_dt64    = t_dt64[:2244]
        Pdb_array = Pdb_array[:,:2244]

##%%
#plt.clf()
#station = 'BBWL'
#bad_ind = np.where(BB[station].T_amp>3e-8)
#BB[station].T_amp[bad_ind] = np.nan
#plt.plot(BB[station].t, BB[station].T_amp)
#
##%%
#import statsmodels.api as sm
#lowess = sm.nonparametric.lowess
#
##f = 1/(doy[-1]-doy[0]) # span equals 1 day; note that 'doy' is time in day of year format
#frac = .5 * 86400/( (BB[station].t[-1] - BB[station].t[0]).astype('int64') )
#n_est = lowess(BB[station].T_amp, BB[station].t, frac=frac, it=5, return_sorted=False)# delta=0.01)
#
#plt.plot(BB[station].t, n_est)
#%%

#    t = np.arange('2017-06-29T00:00Z', '2017-09-25T23:00Z', np.timedelta64(15, 'm'), dtype='Datetime64')
#    t.shape
#    t = t[:Pdb_array.shape[1]]
    
    
    Pow = 10**(Pdb_array/10) # Power not in dB, but vel squared/Hz
    ind = np.where(np.all([freqs>fGHT[0], freqs<fGHT[1]], axis=0))[0]
    GHT_freqs = freqs[ind]

    # Integrate the power per frequency over the range of frequencies identified by "ind"
    #   Integration is the sum of power (not in dB) over some frequency range (fGHT) times
    #   delta F, the spacing between frequency bins.
    GHT_pow = np.sum(Pow[ind,:], axis=0) * np.unique(np.diff(freqs)) # (m/s)^2  Vel squared
    
    BB[ station ] = BBseis(station, t_dt64, 10*np.log10(GHT_pow) ) # Add to a dictionairy item
    
    # Create padded timeseries with first and last values so as not to throw
    #    off the filtered time series too much with end effects.  Then discard padded data.
    pad_amount = 20
    first_ind = np.where(BB[station].powdB > -300)[0][0]
    first = BB[station].powdB[first_ind+5]
    last_ind = np.where(BB[station].powdB > -300)[0][-1]
    last  = BB[station].powdB[last_ind-5 ]
    first_ind = first_ind +1 # Be a little extra conservative- cut out the first and last actual, good points.
    last_ind  = last_ind  -1
    padded = BB[station].powdB.copy()
    padded[:first_ind] = first
    padded[last_ind+1:] = last
#    padded = np.append(np.repeat(first, pad_amount) , BB[station].GHT_powdB, 
#                       np.repeat(last, pad_amount))

    # Define and nan out the tremor amplitude:
    BB[ station ].T_amp = np.sqrt(GHT_pow)
    BB[ station ].T_amp[:first_ind] = np.nan
    BB[ station ].T_amp[last_ind+1:-1] = np.nan

    if station == 'BBGU':
        BB['BBGU'].powdB[1833] = -153
        BB['BBGU'].T_amp[1831:1836] = 8.7e-9

    # Create _power dB Low-pass Filter_ version of the GHT power
    BB[ station ].pdB_LF = signal.filtfilt(filt_b, filt_a, padded)
    BB[ station ].pdB_LF[:first_ind] = np.nan
    BB[ station ].pdB_LF[last_ind+1:-1] = np.nan

    # Nan out the unfiltered power
    BB[ station ].powdB[:first_ind] = np.nan
    BB[ station ].powdB[last_ind+1:-1] = np.nan

    # Filter the tremor amplitude:
    padded = BB[ station ].T_amp.copy()
    first = padded[first_ind+5]
    last  = padded[last_ind-5 ]
    padded[:first_ind] = first
    padded[last_ind+1:] = last
#    BB[ station ].Ta_LF = signal.filtfilt(filt_b, filt_a, padded)
    
    # Implement the lowess smoothing
    bad_ind = np.where(BB[station].T_amp>3e-8)
    BB[station].T_amp[bad_ind] = np.nan
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess
    frac = .5 * 86400/( (BB[station].t[-1] - BB[station].t[0]).astype('int64') )
    BB[ station ].Ta_LF = lowess(BB[station].T_amp, BB[station].t, frac=frac, it=5, return_sorted=False)# delta=0.01)

#    BB[ station ].Ta_LF = median_filter(BB[station].T_amp, size=19)
    
    
    BB[ station ].Ta_LF[:first_ind] = np.nan
    BB[ station ].Ta_LF[last_ind+1:-1] = np.nan

    
    interp_func = interp1d(BB[station].t.astype('datetime64[m]').astype('int64'), 
                           BB[station].Ta_LF, fill_value='extrapolate', 
                           kind='linear')#, bounds_error=False,)
    BB[ station ].Ta_int = interp_func(t_interp.astype('int64'))
#    # PLOTTING
#    fig, ax = plt.subplots(2, sharex=True)
#    ax[0].plot(t_dt64, BB[station].powdB)
#    ax[0].plot(t_dt64, BB[ station ].pdB_LF)
##    ax[0].plot(t, padded)
#    ax[0].set_ylim([-180, -135])
#    ax[0].set_ylabel( 'Power (dB rel. 1 m^2/s^2)' )
#    ax[0].set_title(station)
#    
#    ax[1].plot(t_dt64, BB[station].T_amp)
#    ax[1].plot(t_dt64, BB[station].Ta_LF)
#    ax[1].set_ylim([0, 10e-9])
#    ax[1].set_ylabel('Tremor amp (m/s)')
#    fig.autofmt_xdate()


#sys.exit()

# Saving the objects:
with open('BB_GHT' + str(fGHT) + '.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([BB, fGHT, t_interp], f)

## %% Getting back the objects:
#with open('BB_GHT.pickle', 'rb') as f:  # Python 3: open(..., 'rb')
#    BB, fGHT = pickle.load(f, encoding='latin1')



# %%
#fig, ax = plt.subplots(6, sharex=True, sharey=True)
#ax[0].plot(BB['BBEU'].t, BB['BBEU'].powdB)
#ax[0].plot(BB['BBEU'].t, BB['BBEU'].pdB_LF)
#ax[1].plot(BB['BBWU'].t, BB['BBWU'].powdB)
#ax[1].plot(BB['BBWU'].t, BB['BBWU'].pdB_LF)
##ax[2].plot(BB['BBGU'].t, BB['BBGU'].powdB)
##ax[2].plot(BB['BBGU'].t, BB['BBGU'].pdB_LF)
#ax[3].plot(BB['BBEL'].t, BB['BBEL'].powdB)
#ax[3].plot(BB['BBEL'].t, BB['BBEL'].pdB_LF)
#ax[4].plot(BB['BBWL'].t, BB['BBWL'].powdB)
#ax[4].plot(BB['BBWL'].t, BB['BBWL'].pdB_LF)
#ax[5].plot(BB['BBGL'].t, BB['BBGL'].powdB)
#ax[5].plot(BB['BBGL'].t, BB['BBGL'].pdB_LF)
#ax[0].set_ylim([-170, -150])
#fig.autofmt_xdate()

# %%
# %%


fig = plt.figure()
plt.plot(t_interp, BB['BBEU'].Ta_int - BB['BBWU'].Ta_int, label='BBEU')
plt.plot(t_interp, BB['BBEL'].Ta_int - BB['BBWU'].Ta_int, label='BBEL')
plt.plot(t_interp, BB['BBWL'].Ta_int - BB['BBWU'].Ta_int, label='BBWL')
plt.plot(t_interp, BB['BBGL'].Ta_int - BB['BBWU'].Ta_int, label='BBGL')
plt.legend()
fig.autofmt_xdate()

# %%Look at the relative amplitudes
fig, ax = plt.subplots(ncols=3, sharey=True)
for i in range(3):
    ax[i].plot(BB['BBEU'].t, BB['BBEU'].Ta_LF, label='BBEU', color='darkred')
    ax[i].plot(BB['BBWU'].t, BB['BBWU'].Ta_LF, label='BBWU', color='darkblue')
    ax[i].plot(BB['BBGU'].t, BB['BBGU'].Ta_LF, label='BBGU', color='darkgreen')
    ax[i].plot(BB['BBEL'].t, BB['BBEL'].Ta_LF, label='BBEL', color='lightcoral')
    ax[i].plot(BB['BBWL'].t, BB['BBWL'].Ta_LF, label='BBWL', color='lightblue')
    ax[i].plot(BB['BBGL'].t, BB['BBGL'].Ta_LF, label='BBGL', color='lightgreen')
ax[2].legend()
fig.autofmt_xdate()
ax[0].set_ylim([0e-9, 16e-9])
ax[0].set_xlim(np.array(['2017-07-02', '2017-07-08'],
                        dtype='datetime64' ))
ax[1].set_xlim(np.array(['2017-08-22', '2017-08-28'],
                        dtype='datetime64' ))
ax[2].set_xlim(np.array(['2017-09-10', '2017-09-16'],
                        dtype='datetime64' ))
ax[0].set_ylabel('LP Filtered Tremor Amplitude (m/s)')

#%%
plt.figure(11)
plt.clf()

plt.plot(BB['BBEU'].t, BB['BBEU'].Ta_LF, label='BBEU', color='darkred')
plt.plot(BB['BBWU'].t, BB['BBWU'].Ta_LF, label='BBWU', color='darkblue')
plt.plot(BB['BBGU'].t, BB['BBGU'].Ta_LF, label='BBGU', color='darkgreen')
plt.plot(BB['BBEL'].t, BB['BBEL'].Ta_LF, label='BBEL', color='lightcoral')
plt.plot(BB['BBWL'].t, BB['BBWL'].Ta_LF, label='BBWL', color='lightblue')
plt.plot(BB['BBGL'].t, BB['BBGL'].Ta_LF, label='BBGL', color='lightgreen')
fig.autofmt_xdate()

#%%
fig.clf(10)
fig, ax = plt.subplots(6, sharex=True, sharey=True, num=10)
ax[0].plot(BB['BBEU'].t, BB['BBEU'].T_amp)
ax[0].plot(BB['BBEU'].t, BB['BBEU'].Ta_LF)
#ax[0].plot(t_interp, BB['BBEU'].Ta_int)
ax[1].plot(BB['BBWU'].t, BB['BBWU'].T_amp)
ax[1].plot(BB['BBWU'].t, BB['BBWU'].Ta_LF)
#ax[1].plot(t_interp, BB['BBWU'].Ta_int)
ax[2].plot(BB['BBGU'].t, BB['BBGU'].T_amp)
ax[2].plot(BB['BBGU'].t, BB['BBGU'].Ta_LF)
ax[3].plot(BB['BBEL'].t, BB['BBEL'].T_amp)
ax[3].plot(BB['BBEL'].t, BB['BBEL'].Ta_LF)
ax[4].plot(BB['BBWL'].t, BB['BBWL'].T_amp)
ax[4].plot(BB['BBWL'].t, BB['BBWL'].Ta_LF)
ax[5].plot(BB['BBGL'].t, BB['BBGL'].T_amp)
ax[5].plot(BB['BBGL'].t, BB['BBGL'].Ta_LF)
ax[0].set_ylim([0, 15e-9])
ax[0].set_ylabel('BBEU')
ax[1].set_ylabel('BBWU')
ax[2].set_ylabel('BBGU')
ax[3].set_ylabel('BBEL')
ax[4].set_ylabel('BBWL')
ax[5].set_ylabel('BBGL')

fig.subplots_adjust(hspace=.0001)
fig.autofmt_xdate()
plt.tight_layout()
