#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 19:35:50 2017

@author: timb
"""

import numpy as np
import matplotlib.pyplot as plt



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
#stations = ['BBWU', 'BBEU', 'BBWL', 'BBEL', 'BBGL']
BB = dict() # BB is a dict to contain multiple versions of the class BBseis

fGHT = [1.5, 10]
#fGHT = [0.8, 2]

#%%
class BBseis:
    def __init__(self, name, t, GHT_powdB):
        self.name = name
        self.t = t
        self.powdB = GHT_powdB    # dB rel. 1 (m/s)^2
        self.pdB_LF = np.array([0]) # low frequency filtered version of dB
        self.T_amp = np.array([0]) # m/s : Tremor amplitude
        self.Ta_LF = np.array([0]) # m/s : low frequency filtered version of Tremor amplitude


## %% Getting back the objects:
#with open('Foxx_GHT' + str(fGHT) + '.pickle', 'rb') as f:  # Python 3: open(..., 'rb')
#    BB, fGHT, t_interp = pickle.load(f)#, encoding='latin1')


#%%





    # PLOTTING
fig, ax = plt.subplots()#6, sharex=True)
#fig = plt.figure()

with open('Foxx_GHT[1.5, 10].pickle', 'rb') as f:  # Python 3: open(..., 'rb')
#with open('DL1_GHT[1.5, 10].pickle', 'rb') as f:  # Python 3: open(..., 'rb')
    BB, fGHT, t_interp = pickle.load(f)#, encoding='latin1')

stations = list(BB.keys())

for i in np.arange(len(BB)):
    station = stations[i]

    ax.plot(BB[station].t, BB[station].powdB, label='moulin level?')#, label=station)
#    ax[0].plot(t_dt64, BB[ station ].pdB_LF)
#    ax[0].plot(t, padded)
#    ax[0].set_title(station)
    
#    ax[1].plot(BB[station].t, BB[station].T_amp)
#    ax[1].plot(t_dt64, BB[station].Ta_LF)
#    ax[1].set_ylim([0, 10e-9])

    ax.set_ylabel(station)

#with open('Foxx_GHT[0.8, 1.5].pickle', 'rb') as f:  # Python 3: open(..., 'rb')
with open('DL1_GHT[0.8, 1.5].pickle', 'rb') as f:  # Python 3: open(..., 'rb')
    BB, fGHT, t_interp = pickle.load(f)#, encoding='latin1')
for i in np.arange(len(BB)):
    station = stations[i]
    ax.plot(BB[station].t, BB[station].powdB, label='water flow?')#, label=station)
#    ax[i].set_ylim([-145, -125])
    ax.set_ylabel( 'Power (dB rel. 1 m^2/s^2)' )#station)

with open('DL1_GHT[1, 2.0].pickle', 'rb') as f:  # Python 3: open(..., 'rb')
    BB, fGHT, t_interp = pickle.load(f)#, encoding='latin1')
for i in np.arange(len(BB)):
    station = stations[i]
    ax.plot(BB[station].t, BB[station].powdB, label='intermediate?')#, label=station)
#    ax[i].set_ylim([-145, -125])
    ax.set_ylabel( 'Power (dB rel. 1 m^2/s^2)' )#
    
#ax[0].set_ylabel( 'Power (dB rel. 1 m^2/s^2)' )
#ax[1].set_ylabel('Tremor amp (m/s)')
fig.autofmt_xdate()
fig.legend(loc='lower right')
#ax[0].set_title(str(fGHT))
ax.set_title('Tremor at DL1 in two frequency bands,\nreflecting different source mechanisms')







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
