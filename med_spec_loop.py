#!/opt/anaconda/envs/seisenv/bin python
# -*- coding: utf-8 -*-
"""
med_spec_loop.py

This wrapper script handles the data loading and manipulation, and passes a coarse_duration length of seismic data to the script get_med_spectra, to
calculate the median spectra of many small samples of seismic data with length fine_duration.  get_med_spectra.py returns the median spectra, which is stored
in an array, and finally plotted.

Created on Fri Feb 19 12:06:45 2016
_v1 Mar 27, 2017: TCB
_v2 Nov 29, 2017: Modified to deconvolve the instrument response from the
    waveforms.
_v3 Nov 15, 2018: Modified to use a parameter file for all parameters.
_v5 on October 28, 2022 (yoram terleth): Modified to accomodate long data gaps, bug fixes, clean up. 
Mar 27, 2024 (TCB): Start parallelizing using parsl

Can run script from commandline and overide par file for network and station using:
    >> nohup python -u ./med_spec_loop.py med_spec.par XF BBWL > BBWL.log &    
"""

#%% load packages

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import dates as mdates

import obspy
from obspy.core import read
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

import datetime as dt
import pickle

import glob
import os, time
import sys
import configparser

from clone_inv import clone_inv

# custom function for calculating median spectra 
sys.path.append('/data/stor/basic_data/seismic_data/git_repos/med_spec')
import get_med_spectra_v1
# custom functions for working with obspy UTCDateTimes
from UTCDateTime_funcs import UTCfloor, UTCceil, UTC2dn, UTC2dt64 

# ----------
# Required packages for parsl
import parsl
from parsl import python_app
from parsl.config import Config
from parsl.executors import HighThroughputExecutor
from parsl.providers import LocalProvider

# Configure the parsl executor
activate_env = 'conda activate seisenv24a'
htex_local = Config(
    executors=[
        HighThroughputExecutor(
            max_workers=5,
            provider=LocalProvider(
                worker_init=activate_env
            )
        )
    ],
)
parsl.clear()
parsl.load(htex_local)

parsl_get_med_spec = python_app(get_med_spectra_v1.med_spec)


# @python_app
# def parsl_get_med_spec(tr_trim, pp, Fs_old):
#     freqs, Pdb, Fs_old = get_med_spectra_v1.med_spec(tr_trim, pp, Fs_old)
#     return freqs, Pdb, Fs_old

#------------



# Set the system time to be UTC
os.environ['TZ'] = 'UTC'
time.tzset()


#%% READ PARAMETERS FROM THE PAR FILE


# # Just overwrite the sys.argv for running from within an IDE
# sys.argv = ['',
#             '~/proj/Turner/med_spec_junk/med_spec_short.par',
#             'YG',
#             'SE14']

config = configparser.ConfigParser()
config.read(sys.argv[1])

# config.read('~/proj/Turner/med_spec_junk/med_spec_short.par')
# print(config['DEFAULT'])

# general info
data_source = config['DEFAULT']['data_source']
network = config['DEFAULT']['network']
station = config['DEFAULT']['station']
channel = config['DEFAULT']['channel']

# directories 
data_dir = config['DEFAULT']['data_dir']
resp_file = config['DEFAULT']['resp_file']
out_dir = config['DEFAULT']['out_dir']

# window lenghts (see example par file for details)
pp = {'coarse_duration': float(config['DEFAULT']['pp_coarse_duration']), 
      'coarse_overlap' : float(config['DEFAULT']['pp_coarse_overlap']), 
      'fine_duration'  : float(config['DEFAULT']['pp_fine_duration']), 
      'fine_overlap'   : float(config['DEFAULT']['pp_fine_overlap']) } 

# frequency filters
pre_filt = ( float(config['DEFAULT']['pre_filt0']), 
             float(config['DEFAULT']['pre_filt1']), 
             float(config['DEFAULT']['pre_filt2']), 
             float(config['DEFAULT']['pre_filt3']) )


# temporal range to analyse 
t_start = UTCDateTime(config['DEFAULT']['t_start']) 
t_end = UTCDateTime(config['DEFAULT']['t_end']) 


#%%  ALLOW FOR OPTIONAL COMMAND LINE OVERRIDES OF THE PAR FILE

# sys.argv consists [0] is the file path, [1] is the par_file path, [2] (if it exists) is the network, [3] is the station, [4] is channel

# one override of the par file
if len(sys.argv) == 3: 
    network = sys.argv[2]

# two overrides of the par file
elif len(sys.argv) == 4: 
    network = sys.argv[2]
    station = sys.argv[3]

# three overrides of the par file
elif len(sys.argv) == 5:  
    network = sys.argv[2]
    station = sys.argv[3]
    channel = sys.argv[4]

#%% PRINT OUTPUTS ABOUT THE RUN, FOR THE PURPOSE OF RECORDING IN THE LOG FILE

run_start_time = dt.datetime.now()
print('\n\n' + '===========================================')
print(station + ' run started: ' + '{:%b %d, %Y, %H:%M}'.format(run_start_time) + 'UTC')
print('===========================================' + '\n')

print('  --  PARSL  --  ')

print('Run executed from "' + os.getcwd() + '/"')
print('Run arguments consist of: ' + str(sys.argv))

print('\n\n' + '-------------------------------------------')
print('Start display of parameter file: ' + sys.argv[1])
par_time_sec = os.path.getmtime(sys.argv[1])
print('Parameter file last saved: ' + time.ctime( par_time_sec )  + 'UTC')

print('-------------------------------------------')
with open(sys.argv[1], 'r') as fin:
    print(fin.read())

print('-------------------------------------------')
print('End display of parameter file: ' + sys.argv[1])
print('-------------------------------------------' + '\n\n')

if len(sys.argv) > 2:
    print('-------------------------------------------')
    print('Parameter file overwritten with: ' + str(sys.argv[2:]) )
    print('-------------------------------------------' + '\n\n')



#%% INITIALIZATION

# this section loads one stream to initialize the output arrays  
print('Initialisation: loading first stream...\n\n')

# read from local day_vol files
if data_source=='local':

    # get the response file 
    inv_file = resp_file 
    inv = obspy.read_inventory(inv_file)
    inv = clone_inv(inv, network, station)

    # # subset the inventory to just that which is necessary for script
    inv = inv.select(channel=channel, station=station) 

    # build an array of all the file names of interest
    file_names = glob.glob(data_dir + station + '/*' + channel +'*')
    file_names.sort()
    file_counter = 0

    # read the first stream from the files
    st = read(file_names[0])

    # if there are files that precede our t_start, skip ahead through them
    while(st[0].stats.starttime < t_start):
        if file_counter == 0:
            print('File(s) found that pre-date t_start from the par file.\n')
        file_counter += 1 
        st = read(file_names[file_counter])

    st.merge(method=0)
    print('\n')
    print('Skipped ' + str(file_counter) + ' files.\n')


# read from a client server 
elif data_source != 'local':

    # defien the client 
    fdsn_client = Client(data_source)

    # read in client response file 
    inv = fdsn_client.get_stations(network=network, station=station, channel=channel, 
                                   location='', level="response")
  
    # Read in client data fron the first day
    st = fdsn_client.get_waveforms(network=network, station=station, location='',
                                   channel=channel, starttime=t_start, endtime=t_start+86400)
    
    # initiate token file counter
    file_counter = 0 

else:
    print('No valid data source given!')


# remove response from initalization data 
st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)
print('First considered stream is ' + str(st_IC))

# build arrays
print('\n\n Initialisation: building arrays...')

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

# intialize printing flag 
print_out = True 

# initialize list that will contain all the futures that get launched within the time for loop
# all_freqs_futures = []
# all_Pdb_futures = []
# all_Fs_old_futures = []
all_futures = []

#%% Start the big for loop that will go through each time and each miniseed file and calculate the median powers. These are each time of the spectrogram, each x-axis.
print('\n\nStarting for loop over time...\n\n')

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
        
        # Block progression of the code until all the parsl app_futures finish running
        if file_counter % 3 == 2: # Run 3 days at a a time
            print('Starting blocking at: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()) )
            done = [one_future.result() for one_future in all_futures]
            print('Finished blocking at: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()) )
        # print(done)

        try:

            # reading in from local files
            if data_source=='local':

                # Read in new dayvol file:
                print('About to read in new datafile: ' + file_names[file_counter])            
                st += read(file_names[file_counter])

                # Merge the new and old day volumes
                st.merge(fill_value='interpolate')

                # trim off the part of the merged stream that's already been processed.
                st.trim(starttime=t[i] - pp['coarse_duration'], endtime=st[-1].stats.endtime ) 

            # reading in from client server 
            if data_source!='local':
                
                # read in the next 24 hours of data 
                st = fdsn_client.get_waveforms(network=network, station=station, location='',
                                   channel=channel, starttime=t[i]-pp['coarse_duration'], endtime=t[i]+86400)


        # if there is no file of the querried name, break out of while loop    
        except Exception as e:
            print(e)
            print('Requested file does not exist for '+ str(t[i]))
            break 

        # remove response from newly loaded stream:            
        st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)                  

        # copy instrument corrected trace 
        tr_trim = st_IC[0].copy() 

        # the minus small number and False nearest_sample make the tr include the first data point, but stop just shy of the last data point
        tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - tr_trim.stats.delta)


    ## DEAL WITH DATA GAPS: 
       
    # if the current stream is longer than 5 days, we are getting to a data gap: this statement skips to the where there is renewed data if that is the case. 
    if st[-1].stats.endtime - st[0].stats.starttime > 5 * 24 * 3600:
        if print_out:
            print('Stream in memory is longer than 5 days at ' + str(t[i]))
            print('Stream: ' + str(st))
            print('skipping ahead...')
            print_out = False 

        ## locate the index in t of the end of the stream (which is where there is new data):

        # as long as we are not to the end_date of the stram in memory, skip ahead
        if t[i] < st[-1].stats.endtime: 
            continue 
        
        # stop skipping ahead if we are at the end of the loaded stream
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
    
    
    ## COMPUTE MEDIAN SPECTRA

    # run function            
    # freqs, Pdb, Fs_old = get_med_spectra_v1.med_spec(tr_trim, pp, Fs_old)
    future_result = parsl_get_med_spec(tr_trim, pp, Fs_old)
    # freqs_future, Pdb_future, Fs_old_future
    # all_freqs_futures.append(freqs_future)
    # all_Pdb_futures.append(Pdb_future)
    # all_Fs_old_futures.append(Fs_old_future)
    all_futures.append(future_result)
    print('    AppFuture launched: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()))


    # Save the median spectra into array
    #   Pdb_array[:len(freqs),i] = Pdb[:len(freqs)] 
    
    # end of the for loop 


all_freqs = []
Pdb_list = []
all_Fs_old_futures = []
for future in all_futures:
    (freqs, Pdb, Fs_old) = future.result()
    Pdb_list.append(Pdb)

Pdb_array = np.array(Pdb_list)


#%% PRINT END MESSAGE
print('\n\n' + '===========================================')
print(station + ' run complete: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()))
print('Elapsed time: ' + str(dt.datetime.now() - run_start_time))
print('===========================================' + '\n\n')



Pdb_array[:len(freqs),i] = Pdb[:len(freqs)] 


#%% SAVE THE OUTPUT TO PICKLE FILE 
with open(out_dir + 'mp' + network + '_' + station + '.pickle', 'wb') as f:  
    pickle.dump([t, t_dt64, freqs, Pdb_array, pp, network, station, run_start_time], f)




#%% PRINT END MESSAGE
print('\n\n' + '===========================================')
print(station + ' run complete: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()))
print('Elapsed time: ' + str(dt.datetime.now() - run_start_time))
print('===========================================' + '\n\n')
