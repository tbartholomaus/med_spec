#!/opt/anaconda/envs/seisenv/bin python
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


## current version modified by Yoram on oct 28 2022 as a debug version. Par file and terminal inputs are removed to make running in vs code as a notebook easier. 

"""

#%% %matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
from matplotlib import dates as mdates

import obspy
from obspy.core import read
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
#from obspy import signal

import datetime as dt
import pickle

import glob
import os, time
import sys
import configparser
#import imp
#imp.reload(get_med_spectra)

import get_med_spectra_v1
#from clone_inv import clone_inv
from UTCDateTime_funcs import UTCfloor, UTCceil, UTC2dn, UTC2dt64 # custom functions for working with obspy UTCDateTimes

from clone_inv import clone_inv


os.environ['TZ'] = 'UTC' # Set the system time to be UTC
time.tzset()

#sys.exit()

#%% copy of the par file 

# data_source should be 'local' or 'ETH' or 'IRIS'
data_source= 'local'
# if data_source is not local, then the data_dir and resp_file parameters below are irrelevant.  However, they should be retained in the par file none-the-less.

network= 'YG'
station = 'SW7'
channel = 'HHZ'

# include full paths to data, with final forward slash: /
# script assumes miniseed files are organized in data_dir + station/
data_dir = '/data/stor/basic_data/seismic_data/day_vols/TURNER/'


# the latest response file, OK for use with geobit:
resp_file = '/data/stor/basic_data/seismic_data/day_vols/TURNER/resp/Turner_station_20220923.xml'


# full path to the directory in which the output pickle file should be saved
out_dir = '/data/stor/proj/Turner/med_spec_output/run2020_2022/output_files/new_resp_file/'

# durations in seconds, and overlap in ratios
pp_coarse_duration = 3600.0
pp_coarse_overlap = 0.5
pp_fine_duration = 20.0
pp_fine_overlap = 0.5

pre_filt0= 0.02
pre_filt1= .05
pre_filt2= 110
pre_filt3= 125.

# times in form yyyy-mm-dd
t_start_string = '2020-09-10'
t_end_string = '2021-10-20'



# A set of parameters that define how the script will be run
pp = {'coarse_duration': pp_coarse_duration, 
      'coarse_overlap' : pp_coarse_overlap, 
      'fine_duration'  : pp_fine_duration, 
      'fine_overlap'   : pp_fine_overlap } 

pre_filt = (pre_filt0, pre_filt1, pre_filt2, pre_filt3)


#%% data load 

# #print('Loading time: ' + t[0].strftime('%d %b %Y'))
print('Loading first data')
if data_source=='local':
    # READ IN FROM LOCAL FILES
#    inv_file = resp_file + 'TAKU_station.xml'
    inv_file = resp_file # For mseed files
    inv = obspy.read_inventory(inv_file)
    inv = clone_inv(inv, network, station)

    inv = inv.select(channel=channel, station=station) # subset the inventory to just that which is necessary for script

    file_names = glob.glob(data_dir + station + '/*' + channel +'*')
    file_names.sort()
    file_counter = 0
    st = read(file_names[0])


    t_start = UTCDateTime(t_start_string) # for mseed files       Start and end the timeseries according to the dates during which the station was running.
    t_end = UTCDateTime(t_end_string) # for mseed files       

    
    while(st[0].stats.starttime < t_start): #inv[0][0][0].start_date): # If the first stream is before the t_start, then skip and go to the next file
        if file_counter == 0:
            print('File(s) found that pre-date t_start from the par file.')
        print(' Skipping file: ' + file_names[file_counter])
        file_counter += 1 
        st = read(file_names[file_counter])
    st.merge(method=0)
    print('\n')
    print('Skipped ' + str(file_counter) + ' files.')


# retain the initial stream 
st_IC = st.copy()# .remove_response(inv) 
print('First considered stream is' + str(st_IC))

###########################################################################################################################
## THIS IS THE END OF THE DATA LOADING MODIFICICATIONS ###################################################################
###########################################################################################################################

#%%

# create an np.array of times at which to calculate the median power.
t = np.arange( t_start, t_end, pp['coarse_duration'] * pp['coarse_overlap'])

# Convert a np.array of obspy UTCDateTimes into datenums for the purpose of plotting
t_dt64 = UTC2dt64(t) 

# initialize the sampling rate with zero, to ensure proper running in the first iteration
Fs_old = 0 

#  Find the number of frequencies that the FFT will produce
L = int(st[0].stats.sampling_rate) * pp['fine_duration']
freq_nums = int(2**np.ceil(np.log2( L )) / 2) + 1

# initialize the final array of median powers 
Pdb_array = np.ones( (freq_nums, len(t)) ) * np.nan 

#%% Start the big for loop that will go through each time and each miniseed file and calculate the median powers.  
# These are each time of the spectrogram, each x-axis.

print_out = True 

for i in range(len(t)): 


    ## CURRENT LONG WINDOW: ## 
    # copy instrument corrected trace
    tr_trim = st_IC[0].copy() 

    # the minus small number and False nearest_sample make the tr include the first data point, but stop just shy of the last data point
    tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - tr_trim.stats.delta)
    
    ## IF NEEDED, DEAL WITH FILE CHANGE: ##

    # if the trim is within a coarse seample duration of the end, preload a new one 
    while tr_trim.stats.endtime > st_IC[-1].stats.endtime - pp['coarse_duration']:
        file_counter += 1 
        
        # print where we are to output 
        print("{:>4.0%}".format(float(i)/len(t)) + ' complete.  Loading time: ' + t[i].strftime('%d %b %Y, %H:%M'))
        
        try:
            if data_source=='local':

                # Read in new dayvol file:
                print('About to read in new datafile: ' + file_names[file_counter])            
                st += read(file_names[file_counter])

                # Merge the new and old day volumes
                st.merge(fill_value='interpolate')

                # trim off the part of the merged stream that's already been processed.
                st.trim(starttime=t[i] - pp['coarse_duration'], endtime=st[-1].stats.endtime ) 

        # if there is no file of the querried name, break out of while loop    
        except Exception as e:
            print(e)
            print('requested file does not exist for '+ str(t[i]))
            break 

        # remove response from newly loaded stream:            
        st_IC = st.copy()#.remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)                  

        # copy instrument corrected trace 
        tr_trim = st_IC[0].copy() 

        # the minus small number and False nearest_sample make the tr include the first data point, but stop just shy of the last data point
        tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - tr_trim.stats.delta)


    ## DEAL WITH DATA GAPS: 
       
    # if the current stream is longer than 5 days, we are getting to a data gap:
    # this statement skips to the where there is renewed data if that is the case. 
    if st[-1].stats.endtime - st[0].stats.starttime > 5 * 24 * 3600:
        if print_out:
            print('Stream in memory is longer than 5 days at ' + str(t[i]))
            print('Stream: ' + str(st))
            print('skipping ahead...')
            print_out = False 
        # locate the index in t of the end of the stream (which is where there is new data)
        if t[i] < st[-1].stats.endtime: 
            continue 
        else: 
            print('Done skipping at ' + str(t[i]))
            print('Stream: ' + str(st))
            print('Back to normal processing.')
            print_out = True
    
         
       
    ## DEAL WITH EMPTY OR SHORT STREAMS: 
    
    # If tr_trim does not contain a full coarse_duration of samples, then skip the present coarse_duration and go onto the next iteration of the for loop        
    if tr_trim.stats.npts < pp['coarse_duration'] * int(tr_trim.stats.sampling_rate) * 1:
        print('Incomplete coarse_duration at ' + UTCDateTime.strftime(t[i], "%d %b %y %H:%M") + ': Skipping')
        continue
    

    ## ACTUALLY RUN MED SPEC

    print('Got to the med spec function on ' + str(t[i]))   

    # run function           
    #freqs, Pdb, Fs_old = get_med_spectra_v1.med_spec(tr_trim, pp, Fs_old)

    # Save the median spectra into array
    #Pdb_array[:len(freqs),i] = Pdb[:len(freqs)] 
    

# At the end of the big for loop:
print('\n\n' + '===========================================')
print(station + ' run complete: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()))
#print('Elapsed time: ' + str(dt.datetime.now() - run_start_time))
print('===========================================' + '\n\n')

# %% Pickle the output of the big runs

# Saving the objects:
# with open(out_dir + 'mp' + network + '_' + station + '.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
#     pickle.dump([t, t_dt64, freqs, Pdb_array, pp, network, station, run_start_time], f)