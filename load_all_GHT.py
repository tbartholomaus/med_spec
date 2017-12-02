#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 19:35:50 2017

@author: timb
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
#import matplotlib.mlab as mlab
from matplotlib import dates as mdates

#from obspy.core import read
#from obspy import UTCDateTime
#from obspy import signal

import datetime as dt
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
        self.GHT_powdB = GHT_powdB
        self.GHT_pdLF = np.array([0])



#%%
stations = ['BBWU', 'BBEU', 'BBGU', 'BBWL', 'BBEL', 'BBGL']
BB = dict()

lp_cutoff_period = 6 # hrs


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

#%%

    t = np.arange('2017-06-29T00:00Z', '2017-09-25T23:00Z', np.timedelta64(30, 'm'), dtype='Datetime64')
    t.shape
    t = t[:Pdb_array.shape[1]]
    
    
    Pow = 10**(Pdb_array/10) # Power not in dB, but vel squared/Hz
    fGHT = [1.5, 10]
    ind = np.where(np.all([freqs>fGHT[0], freqs<fGHT[1]], axis=0))[0]
    GHT_freqs = freqs[ind]

    # Integrate the power per frequency over the range of frequencies identified by "ind"
    #   Integration is the sum of power (not in dB) over some frequency range (fGHT) times
    #   delta F, the spacing between frequency bins.
    GHT_pow = np.sum(Pow[ind,:], axis=0) * np.unique(np.diff(freqs))
    
    BB[ station ] = BBseis(station, t_dt64, 10*np.log10(GHT_pow) ) # Add to a dictionairy item
    
    if station == 'BBGU':
        BB['BBGU'].GHT_powdB[1833] = -153

    # Create padded timeseries with first and last values so as not to throw
    #    off the filtered time series too much with end effects.  Then discard padded data.
    pad_amount = 20
    first_ind = np.where(BB[station].GHT_powdB > -300)[0][0]
    first = BB[station].GHT_powdB[first_ind+5]
    last_ind = np.where(BB[station].GHT_powdB > -300)[0][-1]
    last  = BB[station].GHT_powdB[last_ind-5 ]
    padded = BB[station].GHT_powdB.copy()
    padded[:first_ind] = first
    padded[last_ind+1:] = last
#    padded = np.append(np.repeat(first, pad_amount) , BB[station].GHT_powdB, 
#                       np.repeat(last, pad_amount))
    # Create _power dB Low-pass Filter_ version of the GHT power
    BB[ station ].GHT_pdLF = signal.filtfilt(filt_b, filt_a, padded)
    BB[ station ].GHT_pdLF[:first_ind] = np.nan
    BB[ station ].GHT_pdLF[last_ind+1:-1] = np.nan
    
    # PLOTTING
    plt.figure()
    plt.plot(t_dt64, BB[station].GHT_powdB)
    plt.plot(t_dt64, BB[ station ].GHT_pdLF)
#    plt.plot(t, padded)
    plt.ylim([-180, -135])
    plt.ylabel( 'Power (dB rel. 1 m^2/s^2)' )
    plt.title(station)
    



# Saving the objects:
with open('BB_GHT.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([BB], f)

## %% Getting back the objects:
#with open('BB_GHT.pickle', 'rb') as f:  # Python 3: open(..., 'rb')
#    BB = pickle.load(f, encoding='latin1')[0]



# %%
fig, ax = plt.subplots(6, sharex=True, sharey=True)
ax[0].plot(BB['BBEU'].t, BB['BBEU'].GHT_powdB)
ax[0].plot(BB['BBEU'].t, BB['BBEU'].GHT_pdLF)
ax[1].plot(BB['BBWU'].t, BB['BBWU'].GHT_powdB)
ax[1].plot(BB['BBWU'].t, BB['BBWU'].GHT_pdLF)
#ax[2].plot(BB['BBGU'].t, BB['BBGU'].GHT_powdB-10)
#ax[2].plot(BB['BBGU'].t, BB['BBGU'].GHT_pdLF-10)
ax[3].plot(BB['BBEL'].t, BB['BBEL'].GHT_powdB)
ax[3].plot(BB['BBEL'].t, BB['BBEL'].GHT_pdLF)
ax[4].plot(BB['BBWL'].t, BB['BBWL'].GHT_powdB)
ax[4].plot(BB['BBWL'].t, BB['BBWL'].GHT_pdLF)
ax[5].plot(BB['BBGL'].t, BB['BBGL'].GHT_powdB-10)
ax[5].plot(BB['BBGL'].t, BB['BBGL'].GHT_pdLF-10)
ax[0].set_ylim([-170, -150])

