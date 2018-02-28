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
#from matplotlib import dates as mdates

import statsmodels.api as sm

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

stations = ['BBWU', 'BBEU', 'BBGU', 'BBWL', 'BBEL', 'BBGL']
stations = ['XF_BOOM', 'XF_DOST', 'XF_GRAP']
stations = ['XH_FX01', 'XH_FX03', 'XH_FX06', 'XH_FX10', 'XH_FX11', 'XH_FX12']
stations = ['XH_FX01', 'XH_FX03', 'XH_FX06', 'XH_FX10', 'XH_FX11', 'XH_FX12']
#stations = ['7E_DL1']
BB = dict() # BB is a dict to contain multiple versions of the class BBseis

fGHT = [1.5, 10]
fGHT = [0.8, 1.5]
#fGHT = [1, 2.0]

#%%
class BBseis:
    def __init__(self, name, t, GHT_powdB):
        self.name = name
        self.t = t
        self.powdB = GHT_powdB    # dB rel. 1 (m/s)^2
        self.pdB_LF = np.array([0]) # low frequency filtered version of dB
        self.T_amp = np.array([0]) # m/s : Tremor amplitude
        self.Ta_LF = np.array([0]) # m/s : low frequency filtered version of Tremor amplitude


def smooth_tremor(filt_a, filt_b, station, t, GHT_powdB):
    GHT_pow = 10**(GHT_powdB/10)
    T_amp = np.sqrt(GHT_pow) # T_amp is (m/s) because GHT_pow is (m/s)^2

# CLEAN UP DATA AND CREATE T_amp
    # Create padded timeseries with first and last values so as not to throw
    #    off the filtered time series too much with end effects.  Then discard padded data.
            #    pad_amount = 20
    first_ind = np.where(GHT_powdB > -300)[0][0]
    last_ind = np.where(GHT_powdB > -300)[0][-1]

    # Nan out the unfiltered power
    GHT_powdB[:first_ind] = np.nan
    GHT_powdB[last_ind+1:] = np.nan
    # Define and nan out the tremor amplitude:
    T_amp[:first_ind] = np.nan
    T_amp[last_ind+1:] = np.nan

    if station == 'BBGU':
        GHT_powdB[1833] = -153
        T_amp[1831:1836] = 8.7e-9

# LOW PASS FILTERING
    # Create _power dB Low-pass Filter_ version of the GHT power
    # First, extend the data, so that nan's don't play a role.
    first = GHT_powdB[first_ind+5]
    last  = GHT_powdB[last_ind-5 ]
    first_ind = first_ind +1 # Be a little extra conservative- cut out the first and last actual, good points.
    last_ind  = last_ind  -1
    padded = GHT_powdB.copy()
    padded[:first_ind] = first
    padded[last_ind+1:] = last
#    padded = np.append(np.repeat(first, pad_amount) , BB[station].GHT_powdB, 
#                       np.repeat(last, pad_amount))
    pdB_LF = signal.filtfilt(filt_b, filt_a, padded)
    pdB_LF[:first_ind] = np.nan
    pdB_LF[last_ind+1:-1] = np.nan

    # Filter the tremor amplitude:
    padded = T_amp.copy()
    first = padded[first_ind+5]
    last  = padded[last_ind-5 ]
    padded[:first_ind] = first
    padded[last_ind+1:] = last
#    BB[ station ].Ta_LF = signal.filtfilt(filt_b, filt_a, padded)
    
    # Implement the lowess smoothing
    IQR = np.nanpercentile(T_amp, 75) - np.nanpercentile(T_amp, 25)
    pctile_thresh = np.nanmedian(T_amp)+ 5*IQR#np.nanpercentile(T_amp, 98)
    try:
        bad_ind = np.where(T_amp>pctile_thresh)
    except RuntimeWarning:
        pass
    T_amp[bad_ind] = np.nan
    GHT_powdB[bad_ind] = np.nan
    lowess = sm.nonparametric.lowess
    frac = .5 * 86400/( (t[-1] - t[0]).astype('int64') )
    Ta_LF = lowess(T_amp, t, frac=frac, it=5, return_sorted=False)# delta=0.01)

#    BB[ station ].Ta_LF = median_filter(BB[station].T_amp, size=19)
    
    
    Ta_LF[:first_ind] = np.nan
    Ta_LF[last_ind+1:-1] = np.nan
    Ta_LF[np.where(GHT_powdB < -300)[0]] = np.nan
    GHT_powdB[np.where(GHT_powdB < -300)[0]] = np.nan
    
    return GHT_powdB, pdB_LF, T_amp, Ta_LF

#%%

#fGHT = [15, 35]

lp_cutoff_period = 6 # hrs

t_start = np.datetime64('2018-02-01T00:00')
t_end = np.datetime64('2018-02-01T00:00')


for station in stations:
    print(station)
#    station = 'BBWU'#TWLV'


# %% Getting back the objects:
    with open('output_results/mp' + station + '.pickle', 'rb') as f:  # Python 3: open(..., 'rb')
        t, t_dt64, freqs, Pdb_array, pp, data_dir, station = pickle.load(f, encoding='latin1')

    t_start = np.min(np.append(t_dt64, t_start))
    t_end = np.max(np.append(t_dt64, t_end))
    
    sample_period = pp['coarse_duration'] * pp['coarse_overlap'] / 86400 # days  Rate at which the GHT had been sampled
    Nyq = 1/sample_period/2 # 1/days
    lp_cutoff_freq = 1/(lp_cutoff_period/24) # 1/days
    filt_b, filt_a = signal.butter(6, lp_cutoff_freq/Nyq, 'low') # Construct the filter terms used for lowpass filtering the GHT signal


    if station == 'BBGU':
        t         = t[:2244]
        t_dt64    = t_dt64[:2244]
        Pdb_array = Pdb_array[:,:2244]

    
    Pow = 10**(Pdb_array/10) # Power not in dB, but vel squared/Hz
    # Find those indices greater than fGHT[0] and less than fGHT[1]:
    ind = np.where(np.all([freqs>fGHT[0], freqs<fGHT[1]], axis=0))[0] 
    GHT_freqs = freqs[ind] # The list of frequencies over which power will be integrated.

    # Integrate the power per frequency over the range of frequencies identified by "ind"
    #   Integration is the sum of power (not in dB) over some frequency range (fGHT) times
    #   delta F, the spacing between frequency bins.
    GHT_pow = np.sum(Pow[ind,:], axis=0) * np.unique(np.diff(freqs)) # (m/s)^2  Vel squared
    
    # Create the dictionary item BB with class BBseis
    BB[ station ] = BBseis(station, t_dt64, 10*np.log10(GHT_pow) ) # Add to a dictionairy item
        # Dictionary contains station name, time, and the integrated tremor power in dB rel. 1 (m/s)^2
    
    BB[ station ].GHT_powdB, BB[ station ].pdB_LF, BB[ station ].T_amp, BB[ station ].Ta_LF = smooth_tremor(
            filt_a, filt_b, BB[station].name, BB[station].t, BB[station].powdB)
#            station, t, GHT_powdB, filt_a, filt_b)

    

#    # PLOTTING
##    fig, ax = plt.subplots(2, sharex=True)
#    ax[0].plot(t_dt64, BB[station].powdB)
##    ax[0].plot(t_dt64, BB[ station ].pdB_LF)
##    ax[0].plot(t, padded)
##    ax[0].set_ylim([-180, -135])
#    ax[0].set_ylabel( 'Power (dB rel. 1 m^2/s^2)' )
#    ax[0].set_title(station)
#    
#    ax[1].plot(t_dt64, BB[station].T_amp)
##    ax[1].plot(t_dt64, BB[station].Ta_LF)
##    ax[1].set_ylim([0, 10e-9])
#    ax[1].set_ylabel('Tremor amp (m/s)')
#    fig.autofmt_xdate()

#%%
t_start = np.datetime64(t_start, 'D')
t_end = np.datetime64(t_end, 'D') + 1
t_interp = np.arange(t_start, t_end, 
                     np.timedelta64(15, 'm'), dtype='Datetime64')



    # PLOTTING
fig, ax = plt.subplots(2, sharex=True)
    
    
for station in BB:
    interp_func = interp1d(BB[station].t.astype('datetime64[m]').astype('int64'), 
                           BB[station].Ta_LF, fill_value='extrapolate', 
                           kind='linear')#, bounds_error=False,)
    BB[ station ].Ta_int = interp_func(t_interp.astype('int64'))

    ax[0].plot(BB[station].t, BB[station].powdB, label=station)
#    ax[0].plot(t_dt64, BB[ station ].pdB_LF)
#    ax[0].plot(t, padded)
#    ax[0].set_title(station)
    
    ax[1].plot(BB[station].t, BB[station].T_amp)
#    ax[1].plot(t_dt64, BB[station].Ta_LF)
#    ax[1].set_ylim([0, 10e-9])

ax[0].set_ylim([-150, -100])
ax[0].set_ylabel( 'Power (dB rel. 1 m^2/s^2)' )
ax[1].set_ylabel('Tremor amp (m/s)')
fig.autofmt_xdate()
fig.legend()
ax[0].set_title(str(fGHT))

#sys.exit()

# Saving the objects:
with open('Foxx_GHT' + str(fGHT) + '.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
#with open(station + '_GHT' + str(fGHT) + '.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([BB, fGHT, t_interp], f)

## %% Getting back the objects:
#with open('Foxx_GHT' + str(fGHT) + '.pickle', 'rb') as f:  # Python 3: open(..., 'rb')
#    BB, fGHT, t_interp = pickle.load(f)#, encoding='latin1')


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
sys.exit()

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
