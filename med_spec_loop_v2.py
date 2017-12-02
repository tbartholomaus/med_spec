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
    
Can read input station from commandline using:
    >> nohup python -u ./med_spec_loop_v2.py BBWL > BBWL.log &    

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
import os, time
import fnmatch
import sys
#import imp
#imp.reload(get_med_spectra)

import get_med_spectra_v1
from clone_inv import clone_inv
from UTCDateTime_funcs import UTCfloor, UTCceil, UTC2dn, UTC2dt64 # custom functions for working with obspy UTCDateTimes

os.environ['TZ'] = 'UTC' # Set the system time to be UTC
time.tzset()

#sys.exit()
#%%
## This next block of code is the Lemon Creek experiment, run on ibest
#network = 'LM'
#station = sys.argv[1]#'BBWL'#TWLV'
#data_dir = '/mnt/lfs2/tbartholomaus/Seis_data/day_vols/LEMON/'

# This next block of code is for the Moscow Mtn test run on my laptop
network = 'XX'
network = 'LM'
#station = 'BBGL'
#data_dir = '/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/LemonCrk_GHT/Seis_analysis/wf_data/Moscow_Mtn/GB/'
station = 'BBEL'#TWLV'
data_dir = '/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/LemonCrk_GHT/Seis_analysis/wf_data/Moscow_Mtn/NM_together/'
data_dir = '/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/LemonCrk_GHT/Seis_analysis/wf_data/testing/'


resp_dir = '../RESP/'
#data_dir = '/mnt/gfs/tbartholomaus/Seis_data/day_vols/data_temp/LEMON/on_ice'

# A set of parameters that define how the script will be run
pp = {'coarse_duration': 3600.0,  # s
      'coarse_overlap' : 0.5,   # Ratio of overlap
      'fine_duration'  : 20.0,   # s
      'fine_overlap'   : 0.5}   # Ratio of overlap

#data_dir = '/Volumes/disk_staff/timb/Seis_data/day_vols/TAKU/day_volumes/' + station + '/'
#data_dir = '../mseed_files/'

# file_names unless we read the individual 10 min geobit files below:
file_names = glob.glob(data_dir + station + '/*HZ*')

file_names.sort()
if station == 'BBGL' and len(file_names)>5:
    file_names = file_names[2:] # For BBGL  "2" when it's day files, 11 when its 10 min files
elif station == 'BBGU':
    file_names = file_names[1:] # For BBGU
elif station == 'BBWU':
    file_names = file_names[1:] # For BBGU
elif station == 'BBEL':
    file_names = file_names[2:] # For BBEL
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

    #data_dir = '/mnt/gfs/tbartholomaus/Seis_data/RAW/LEMON/SV02/msd_10min/'
    #For working with all the Geobit data
    if network == 'LM':
        file_names = []
        for root, dirnames, fnames in os.walk('/mnt/gfs/tbartholomaus/Seis_data/RAW/LEMON/SV02/msd_10min/'+station):
            for fname_cnt in fnmatch.filter(fnames, '*.HHZ.mseed'):
                file_names.append(os.path.join(root, fname_cnt))
    if station =='BBGL':
        file_names = file_names[11:]

elif station[:3] == 'BBW' or station[:3] == 'BBE' or station == 'UI05':
    pre_filt = (0.01, .02, 210, 225.)
    resp_file_nm = resp_dir + "XX.UI02.resp_171112/XX.UI02.HHZ.resp"
    inv = obspy.read_inventory(resp_file_nm)


inv[0][0][0].start_date = UTCDateTime("2017-6-29") # All stations were installed on 6/29, after 8 local
# add another inventory object with the same response, but with the correct station names
#    sta_chan_id = gb[0].get_id()
#    inv_gb = clone_inv(inv_gb, sta_chan_id[:2], sta_chan_id[3:7])
inv = clone_inv(inv, network, station)

    
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
st.merge(fill_value='interpolate')
st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)

#st1 = read('./' + station + '.**..*HZ.2015.248')

# create an np.array of times at which to calculate the median power.
#    The actual, calculated median powers will be calculated for the coarse_duration
#    beginning at this time.
t = np.arange( UTCfloor(st[0].stats.starttime), UTCceil(last_time), 
              pp['coarse_duration'] * pp['coarse_overlap'])#, dtype='datetime64[D]')
#t = t[:-2]
t_dt64 = UTC2dt64(t) # Convert a np.array of obspy UTCDateTimes into datenums for the purpose of plotting

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

print(st[0].stats)
print(st_IC[0].stats)

flag_up = False

for i in range(len(t)): # Loop over all the t's, however, the for loop will never complete
#     the loop ends when file_counter == len(file_names).  Perhaps a while loop would be more elegant 
#%%    
#    tr = st[0]
    tr_trim = st_IC[0].copy() # copy instrument corrected trace

    # the minus small number and False nearest_sample make the tr include the first data point, but stop just shy of the last data point
    tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - 0.00001, nearest_sample=False)
    
#    print()
#    print(i)
#    print(t[i])

#    # If you're at the end of the day volume, and the duration of tr_trim is not coarse_duration long:
#    while t[i] + pp['coarse_duration'] > tr_trim.stats.endtime + 0.01:
#    # If the trimmed trace is shorter than it ought to be, load in a new file
#    while tr_trim.stats.npts/tr_trim.stats.sampling_rate < pp['coarse_duration']:

    # If the trimmed trace ends within pp['coarse_duration'] of the end of the 
    #   instrument corrected trace, then reload the next file.
    #   This keeps away from the end of the st_IC, which is tapered.
    while tr_trim.stats.endtime > st_IC[0].stats.endtime - 2*pp['coarse_duration']:        
        file_counter += 1
        print (t[i])
        
        if file_counter >= len(file_names):
            break # break out of the for loop when there are no more files to load.
        
        print("{:>4.0%}".format(float(i)/len(t)) + ' complete.  Loading file: ' + file_names[file_counter])
         # Read in another day volume as another trace, pre_filter it (which is
         #      critical, otherwise spectra are crap, remove its instrument
         #      response at the same time, and then add that trace after the
         #      existing trace, into the stream "st"
        st += read(file_names[file_counter])
        st.merge(fill_value='interpolate')#method=0) # Merge the new and old day volumes
        print(st[0].stats)
        st.trim(starttime=t[i], endtime=st[0].stats.endtime ) # trim off the part of the merged stream that's already been processed.
        
        
        st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)
        tr_trim = st_IC[0].copy() # copy the trace in st
        tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - 0.00001, nearest_sample=False)

        print(st_IC[0].stats)
        
#        flag_up = True
#    print(tr_trim)
#    print(st)
    
    # If tr_trim does not contain a full coarse_duration of samples, then skip
    #     the present coarse_duration and go onto the next iteration of the for
    #     loop.
    if tr_trim.stats.npts < pp['coarse_duration'] * int(tr_trim.stats.sampling_rate) * 1:
        print('Incomplete coarse_duration at ' + UTCDateTime.strftime(t[i], "%d %b %y %H:%M") + ': Skipping')

#        if flag_up:
#            sys.exit()
            
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

# %% Pickle the output of the big runs

# Saving the objects:
with open('mp' + station + '_test2.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([t, t_dt64, freqs, Pdb_array, pp, data_dir, station], f)

