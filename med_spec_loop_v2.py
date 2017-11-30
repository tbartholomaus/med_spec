# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 12:06:45 2016
Completed _v1 on March 27, 2017

@author: tbartholomaus

Wrapper script that calculates the median spectra from a timeseries of seismic
data, following the methods described in Bartholomaus et al., 2015, in GRL.

This wrapper script handles the data loading and manipulation, and passes a 
coarse_duration length of seismic data to the script get_med_spectra, to
calculate the median spectra of many small samples of seismic data with length 
fine_duration.  get_med_spectra.py returns the median spectra, which is stored
in an array, and finally plotted.

_v2 Nov 29, 2017: Modified to deconvolve the instrument response from the
    waveforms.

"""

#%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
from matplotlib import dates as mdates

import obspy
from obspy.core import read
from obspy import UTCDateTime
#from obspy import signal

import datetime as dt
import pickle

import glob
import os
import fnmatch
import sys
#import imp
#imp.reload(get_med_spectra)

import get_med_spectra_v1
from clone_inv import clone_inv
from UTCDateTime_funcs import UTCfloor, UTCceil, UTC2dn # custom functions for working with obspy UTCDateTimes

#sys.exit()
#%%
station = 'BBGL'#TWLV'
#data_dir = '/mnt/lfs2/tbartholomaus/Seis_data/day_vols/LEMON/'
data_dir = '/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/LemonCrk_GHT/Seis_analysis/wf_data/Moscow_Mtn/GB/'

#station = 'UI05'#TWLV'
##data_dir = '/mnt/lfs2/tbartholomaus/Seis_data/day_vols/LEMON/'
#data_dir = '/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/LemonCrk_GHT/Seis_analysis/wf_data/Moscow_Mtn/NM_together/'


resp_dir = '../RESP/'
#data_dir = '/mnt/gfs/tbartholomaus/Seis_data/day_vols/data_temp/LEMON/on_ice'

# A set of parameters that define how the script will be run
pp = {'coarse_duration': 3600.0,  # s
      'coarse_overlap' : 0.5,   # Ratio of overlap
      'fine_duration'  : 20.0,   # s
      'fine_overlap'   : 0.5}   # Ratio of overlap

#data_dir = '/Volumes/disk_staff/timb/Seis_data/day_vols/TAKU/day_volumes/' + station + '/'
#data_dir = '../mseed_files/'

#file_names = []
#for root, dirnames, fnames in os.walk('/mnt/gfs/tbartholomaus/Seis_data/day_vols/data_temp/LEMON/on_ice/BBGU'):
#    for fname_cnt in fnmatch.filter(fnames, '*.HHZ.mseed'):
#        file_names.append(os.path.join(root, fname_cnt))

file_names = glob.glob(data_dir + station + '/*HZ*')

file_names.sort()
if station == 'BBGL' and len(file_names)>5:
    file_names = file_names[11:] # For BBGL
elif station == 'BBGU':
    file_names = file_names[1:] # For BBGU
elif station == 'BBWU':
    file_names = file_names[1:] # For BBGU
elif station == 'BBEL':
    file_names = file_names[1:] # For BBEL
elif station == 'BBWL':
    file_names = file_names[1:] # For BBWL
file_counter = 0

# %% LOAD IN THE INSTRUMENT CORRECTIONS TO CREATE INV OBJECTS

if station[:3] == 'BBG':
    pre_filt = (0.1, .4, 80, 85.)
    resp_file_gb = resp_dir + "RESP_C100_SRi32_200sps.rsp"

    inv = obspy.read_inventory(resp_file_gb)

    inst_sens = inv[0][0][0].response.instrument_sensitivity.value
    # factor change is sqrt( 10**(dB/10) )
    # factor_change is a number greater than 1 because the GB sensor has a spectral
    #   power greater than the NM sensor.
    factor_change = 3.0
    inv[0][0][0].response.instrument_sensitivity.value = inst_sens*factor_change
    inv[0][0][0].response.response_stages[5].stage_gain = factor_change
    
    #    sta_chan_id = gb[0].get_id()
    # add another inventory object with the same response, but with the correct station names
    #    inv_gb = clone_inv(inv_gb, sta_chan_id[:2], sta_chan_id[3:7])
    inv = clone_inv(inv, 'XX', station)
    
elif station[:3] == 'BBW' or station[:3] == 'BBE' or station == 'UI05':
    pre_filt = (0.01, .02, 210, 225.)
    resp_file_nm = resp_dir + "XX.UI02.resp_171112/XX.UI02.HHZ.resp"
    inv = obspy.read_inventory(resp_file_nm)

    inv = clone_inv(inv, 'XX', station)

inv[0][0][0].start_date = UTCDateTime("2017-6-29") # All stations were installed on 6/29, after 8 local


    
#%%

# Read in and find the time of the last data record of the last miniseed file
st = read(file_names[-1])
last_time = st[-1].stats.endtime

#day_num = 247
#
#day_vol = './' + station + '.**..*HZ.2015.' + str(day_num)
print('Loading file ' + file_names[0])
# Read in and remove instrument response from first file
st = read(file_names[0])
st.merge(fill_value='interpolate').remove_response(inventory=inv, output="VEL")

#st1 = read('./' + station + '.**..*HZ.2015.248')

# create an np.array of times at which to calculate the median power.
#    The actual, calculated median powers will be calculated for the coarse_duration
#    beginning at this time.
t = np.arange( UTCfloor(st[0].stats.starttime), UTCceil(last_time), 
              pp['coarse_duration'] * pp['coarse_overlap'])#, dtype='datetime64[D]')
t = t[:-2]

Fs_old = 0 # initialize the sampling rate with zero, to ensure proper running in the first iteration

#%% Find the number of frequencies that the FFT will produce, and initialize the output array
L = int(st[0].stats.sampling_rate) * pp['fine_duration']
freq_nums = int(2**np.ceil(np.log2( L )) / 2) + 1#1025 #2049 # Number of frequencies output.  2049 for 200Hz data, 1025 for 100Hz data.
Pdb_array = np.ones( (freq_nums, len(t)) ) * -500.0 # initialize the final array of median powers with -500 db/Hz

#%% Start the big for loop that will go through each time and each miniseed 
#       file and calculate the median powers.
run_start_time = dt.datetime.now()
print('\n\n' + '===========================================')
print(station + ' run started: ' + '{:%b %d, %Y, %H:%M}'.format(run_start_time))
print('===========================================' + '\n\n')

for i in range(len(t)): # Loop over all the t's, however, the for loop will never complete
#     the loop ends when file_counter == len(file_names).  Perhaps a while loop would be more elegant 
#%%    
#    tr = st[0]
    tr_trim = st[0].copy() # copy the trace in st

    # the minus small number and False nearest_sample make the tr include the first data point, but stop just shy of the last data point
    tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - 0.00001, nearest_sample=False)
    
#    print()
#    print(i)
#    print(t[i])

    # If you're at the end of the day volume, and the duration of tr_trim is not coarse_duration long:
    while t[i] + pp['coarse_duration'] > tr_trim.stats.endtime + 0.01:
        file_counter += 1
        
        if file_counter >= len(file_names):
            break # break out of the for loop when there are no more files to load.
        
        print("{:>4.0%}".format(float(i)/len(t)) + ' complete.  Loading file: ' + file_names[file_counter])
         # Read in another day volume as another trace, pre_filter it (which is
         #      critical, otherwise spectra are crap, remove its instrument
         #      response at the same time, and then add that trace after the
         #      existing trace, into the stream "st"
        st = read(file_names[file_counter]).remove_response(inventory=inv, 
            output="VEL", pre_filt=pre_filt)
        st.merge(fill_value='interpolate')#method=0) # Merge the new and old day volumes

        tr_trim = st[0].copy() # copy the trace in st
        tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - 0.00001, nearest_sample=False)

        st.trim(starttime=t[i], endtime=st[0].stats.endtime ) # trim off the part of the merged stream that's already been processed.
    
#    print(tr_trim)
#    print(st)
    
    # If tr_trim does not contain a full coarse_duration of samples, then skip
    #     the present coarse_duration and go onto the next iteration of the for
    #     loop.
    if tr_trim.stats.npts < pp['coarse_duration'] * int(tr_trim.stats.sampling_rate) * 1:
        print('Incomplete coarse_duration at ' + UTCDateTime.strftime(t[i], "%d %b %y %H:%M") + ': Skipping')
        continue
    
#    print('Calculating median spectra for ' + UTCDateTime.strftime(t[i], "%d/%m/%y %H:%M"))
    # pass the seismic data with coarse_duration to the helper function for
    #     calculation of the median spectra.
    freqs, Pdb, Fs_old = get_med_spectra_v1.med_spec(tr_trim, pp, Fs_old)
    Pdb_array[:,i] = Pdb[:freq_nums] # Save the median spectra into an array
    
#makeerror # junk command that doesn't exist, and will throw an error

# At the end of the big for loop:
print('\n\n' + '===========================================')
print(station + ' run complete: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()))
print('Elapsed time: ' + str(dt.datetime.now() - run_start_time))
print('===========================================' + '\n\n')

t_datenum = UTC2dn(t) # Convert a np.array of obspy UTCDateTimes into datenums for the purpose of plotting
# %% Pickle the output of the big runs

# Saving the objects:
with open('mp' + station + '.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([t, freqs, Pdb_array, pp, data_dir, station], f)

# %% Getting back the objects:
#with open('output/mpBBWU.pickle', 'rb') as f:  # Python 3: open(..., 'rb')
#    t, freqs, Pdb_array, pp, data_dir, station = pickle.load(f)

#%% Plot the output of the big runs as median spectrograms

mask_val1 = Pdb_array<=-300
mask_val2 = np.isinf(Pdb_array)
Pdb_array_mask = np.ma.masked_where(np.any([mask_val1, mask_val2], axis=0), Pdb_array)

t_datenum = UTC2dn(t) # Convert a np.array of obspy UTCDateTimes into datenums for the purpose of plotting

#plt.imshow(np.log10(Pxx_vals[0:-2,]), extent = [0, len(file_names), freqs[1], freqs[freq_nums-1]])
fig, ax = plt.subplots()#figsize=(8, 4))
qm = ax.pcolormesh(t_datenum, freqs, Pdb_array_mask, cmap='YlOrRd')#, extent = [0, len(file_names), freqs[1], freqs[freq_nums-1]])
ax.set_yscale('log')
ax.set_ylim([0.1, 250])
#ax.set_ylim([0.5, 20])
ax.set_ylabel('Frequency (Hz)')

# Set the date limits for the plot, to make each station's plots consistent
#ax.set_xlim(mdates.date2num([dt.date(2017, 6, 25), dt.date(2017, 9, 30)]))
#ax.set_xlim(mdates.date2num([dt.date(2017, 7, 3), dt.date(2017, 7, 8)]))
ax.set_xlim(mdates.date2num([dt.datetime(2017, 10, 4, 23,0,0), dt.datetime(2017, 10, 7, 1, 0, 0)]))

# Format the xaxis of the pcolormesh to be dates
ax.xaxis.set_major_locator(mdates.AutoDateLocator())
ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
fig.autofmt_xdate()

qm.set_clim(vmin=-210, vmax=-140)
cb = plt.colorbar(qm, ticks=np.arange(-210,-140, 20))
cb.set_label('Power (dB, rel. 1 (m/s)^2/Hz)')
plt.title(station)

# #%%
plt.savefig('Spec_' + station + '_fld', dpi=150) # _ght, _fld
#plt.show()

#%% Obsolete junk down here

#==============================================================================
# ind3Hz = np.where(freqs>3)[0][0]
# 
# Pxx_atFreq = Pxx_vals[ind3Hz,:]
# plt.plot(np.log10(Pxx_atFreq))
# plt.title(station)
# plt.savefig(station+'_at3Hz')
# plt.show()
#==============================================================================

# %%
#    
##Pxx, freqs = signal.psd(tr, NFFT=256, Fs=st[0].stats.sampling_rate, detrend='linear', 
##                        window=mlab.window_hanning)
#Pxx, freqs = mlab.psd(tr.data, NFFT=2048, Fs=st[0].stats.sampling_rate, detrend='linear', 
#                        scale_by_freq=True, window=mlab.window_hanning)
#plt.plot(freqs[0:50], Pxx[0:50])
##%%
## ticker+0
#                        
## tr = st[0]
## tr = tr.copy()
#
## tr.trim(ticker, ticker+process_window)
## st[0].stats.endtime
## plt.plot(tr.data)
#tr.plot()
## plt.plot(freqs, Pxx)
## plt.plot(med_val)
## mlab.window_hanning(tr)
## window = mlab.window_hanning(tr)
## tr.filter('bandpass', freqmin = 1.5, freqmax = 10)
#
#med_val
#
#plt.plot(med_val)
