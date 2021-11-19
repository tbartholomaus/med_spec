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

_v2 Nov 29, 2017: Modified to deconvolve the instrument response from the
    waveforms.
_v3 Nov 15, 2018: Modified to use a parameter file for all parameters.
    
Can read input station from commandline using:
    >> nohup python -u ./med_spec_loop_v3.py med_spec.par XF BBWL > BBWL.log &    

"""

#%matplotlib inline

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

#%% READ FROM THE PAR FILE
config = configparser.ConfigParser()
config.read(sys.argv[1])

data_source = config['DEFAULT']['data_source']
network = config['DEFAULT']['network']
station = config['DEFAULT']['station']
channel = config['DEFAULT']['channel']


data_dir = config['DEFAULT']['data_dir']
resp_file = config['DEFAULT']['resp_file']
out_dir = config['DEFAULT']['out_dir']


# A set of parameters that define how the script will be run
pp = {'coarse_duration': float(config['DEFAULT']['pp_coarse_duration']), # 3600.0,  # s
      'coarse_overlap' : float(config['DEFAULT']['pp_coarse_overlap']), # 0.5,   # Ratio of overlap
      'fine_duration'  : float(config['DEFAULT']['pp_fine_duration']), # 20.0,   # s
      'fine_overlap'   : float(config['DEFAULT']['pp_fine_overlap']) } # 0.5}   # Ratio of overlap

pre_filt = ( float(config['DEFAULT']['pre_filt0']), 
             float(config['DEFAULT']['pre_filt1']), 
             float(config['DEFAULT']['pre_filt2']), 
             float(config['DEFAULT']['pre_filt3']) )#0.01, .02, 90, 100.)


t_start = config['DEFAULT']['t_start']
t_end = config['DEFAULT']['t_end']

# ALLOW FOR OPTIONAL ARGV OVERRIDES OF THE PAR FILE
# sys.argv consists [0] is the file path, [1] is the par_file path, 
#   [2] (if it exists) is the network, [3] is the station, [4] is channel
if len(sys.argv) == 3: # i.e., one override of the par file
    network = sys.argv[2]#'7E'
elif len(sys.argv) == 4: # i.e., two overrides of the par file
    network = sys.argv[2]#'7E'
    station = sys.argv[3]#'BBWL'#TWLV'
elif len(sys.argv) == 5: # i.e., two overrides of the par file
    network = sys.argv[2]#'7E'
    station = sys.argv[3]#'BBWL'#TWLV'
    channel = sys.argv[4]


#
#network = 'XX'
#station = 'BBGU'
#resp_file = "/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/Instruments/geobit/180305_Resp/V1/RESP_C100SRi32V1_200sps_MIN.rsp"
#data_dir = '/Users/timb/Documents/syncs/OneDrive - University of Idaho/RESEARCHs/Instruments/geobit/1705_Geobit_testing/from_thumb/msd_files/'#CWU201516/'#/SV03/'
#pre_filt = (0.1, .2, 90, 100.)


#%% PRINT OUTPUTS ABOUT THE RUN, FOR THE PURPOSE OF RECORDING IN THE LOG FILE
run_start_time = dt.datetime.now()
print('\n\n' + '===========================================')
print(station + ' run started: ' + '{:%b %d, %Y, %H:%M}'.format(run_start_time))
print('===========================================' + '\n')

print('Run executed from "' + os.getcwd() + '/"')
print('Run arguments consist of: ' + str(sys.argv))


print('\n\n' + '-------------------------------------------')
print('Start display of parameter file: ' + sys.argv[1])
par_time_sec = os.path.getmtime(sys.argv[1])
print('Parameter file last saved: ' + time.ctime( par_time_sec ) )
print('-------------------------------------------')
with open(sys.argv[1], 'r') as fin:
#with open(par_file[0], 'r') as fin:
    print(fin.read())
print('-------------------------------------------')
print('End display of parameter file: ' + sys.argv[1])
print('-------------------------------------------' + '\n\n')

if len(sys.argv) > 2:
    print('-------------------------------------------')
    print('Parameter file overwritten with: ' + str(sys.argv[2:]) )
    print('-------------------------------------------' + '\n\n')

#print('Loading time: ' + t[0].strftime('%d %b %Y'))
print('Loading first data')
if data_source=='local':
    # READ IN FROM LOCAL FILES
#    inv_file = resp_file + 'TAKU_station.xml'
    inv_file = resp_file # For mseed files
    inv = obspy.read_inventory(inv_file)
    inv = clone_inv(inv, network, station)

    inv = inv.select(channel=channel, station=station) # subset the inventory to just that which is necessary for script

    file_names = glob.glob(data_dir + station + '/*' + channel[-2:] +'*')
    file_names.sort()
    file_counter = 0
    st = read(file_names[0])

#    t_start = inv[0][0].start_date # Start and end the timeseries according to the dates during which the station was running.
#    t_end = inv[0][0].end_date
    t_start = UTCDateTime(config['DEFAULT']['t_start']) # for mseed files       Start and end the timeseries according to the dates during which the station was running.
    t_end = UTCDateTime(config['DEFAULT']['t_end']) # for mseed files       
#    inv[0][0][0].start_date = t_start
#    inv[0][0][0].end_date = t_end
    
    while(st[0].stats.starttime < inv[0][0][0].start_date): # If the first stream is before the t_start, then skip and go to the next file
        if file_counter == 0:
            print('File(s) found that pre-date t_start from the par file.')
        print(' Skipping file: ' + file_names[file_counter])
        file_counter += 1 
        st = read(file_names[file_counter])
    st.merge(method=0)
    print('\n')

if data_source!='local':
    fdsn_client = Client(data_source)#'ETH')#'IRIS')

    # READ IN FROM FDSN CLIENT
    inv = fdsn_client.get_stations(network=network, station=station, channel=channel, 
                                   location='', level="response")
    t_start = inv[0][0].start_date
    t_end = inv[0][0].end_date
    
    # Read in and remove instrument response from first day
    st = fdsn_client.get_waveforms(network=network, station=station, location='',
                                   channel=channel, starttime=t_start, endtime=t_start+86400)


#sys.quit()

#%%

st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)



# create an np.array of times at which to calculate the median power.
#    The actual, calculated median powers will be calculated for the coarse_duration
#    beginning at this time.
t = np.arange( t_start, t_end, pp['coarse_duration'] * pp['coarse_overlap'])#, dtype='datetime64[D]')
#t = t[:-2]
t_dt64 = UTC2dt64(t) # Convert a np.array of obspy UTCDateTimes into datenums for the purpose of plotting

Fs_old = 0 # initialize the sampling rate with zero, to ensure proper running in the first iteration




#%% Find the number of frequencies that the FFT will produce, and initialize the output array
L = int(st[0].stats.sampling_rate) * pp['fine_duration']
freq_nums = int(2**np.ceil(np.log2( L )) / 2) + 1#1025 #2049 # Number of frequencies output.  2049 for 200Hz data, 1025 for 100Hz data.
Pdb_array = np.ones( (freq_nums, len(t)) ) * -500.0 # initialize the final array of median powers with -500 db/Hz
#Pdb_array = np.ones( (2049, len(t)) ) * -500.0 # initialize the final array of median powers with -500 db/Hz

#%% Start the big for loop that will go through each time and each miniseed 
#       file and calculate the median powers.  These are each time of the
#       spectrogram, each x-axis.

flag_up = False

#for i in np.arange(950, 1200):#range(len(t)): # Loop over all the t's, however, the for loop will never complete
for i in range(len(t)): # Loop over all the t's, however, the for loop will never complete
#     the loop ends when file_counter == len(file_names).  Perhaps a while loop would be more elegant 
##%%    
#    tr = st[0]
#    print('>> Starting new loop: ' +str(i) )
#    time_temp = dt.datetime.now()
#    IC_counter += 1
#    
#    if IC_counter == 3: # Note that this threshold IC_counter value is ad-hoc and may not work if different amounts of IC tapering are done, or if the pp['coarse_duration'] is different than 3600
#        # Now that t[i] has moved on IC_counter timesteps away from the 
#        #   beginning of st, remove the instrument response, which will also
#        #   end up with a tapered st_IC.  We don't want to look at the data
#        #   which has been tapered.
#        st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)
#        print('--- Re-initialize st_IC -----')
#        print(t[i])
#        print(st_IC)
        
    tr_trim = st_IC[0].copy() # copy instrument corrected trace
    # the minus small number and False nearest_sample make the tr include the first data point, but stop just shy of the last data point
#    tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - 0.00001, nearest_sample=False)
    tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - tr_trim.stats.delta)

#    print()
#    print(i)
#    print(t[i])    
#    print(tr_trim.stats.npts)
#    print(st)
#    print(st_IC)
#    print(tr_trim)


#    # If you're at the end of the day volume, and the duration of tr_trim is not coarse_duration long:
#    while t[i] + pp['coarse_duration'] > tr_trim.stats.endtime + 0.01:
#    # If the trimmed trace is shorter than it ought to be, load in a new file
#    while tr_trim.stats.npts/tr_trim.stats.sampling_rate < pp['coarse_duration']:

    # If the trimmed trace ends within 2*pp['coarse_duration'] of the end of  
    #   the data stream, then reload the next file.
    #   This keeps away tr_trim away from the end of the st_IC, which is tapered.
    while tr_trim.stats.endtime > st_IC[-1].stats.endtime - pp['coarse_duration']:
        file_counter += 1
#        print('--- Try to load in a new stream ---')
#        print (t[i])
        
#        if t[i]+86400 >= t_end:
#            break # break out of the for loop when there are no more files to load.
        
        print("{:>4.0%}".format(float(i)/len(t)) + ' complete.  Loading time: ' + t[i].strftime('%d %b %Y, %H:%M'))
        
        try:
            if data_source=='local':
                # Read in local file
                        # Read in another day volume as another trace, and merge it 
                        #   into the stream "st".  When the st is instrument corrected, t[i]
                        #   will have moved on, away from the beginning of the st.
#                print('About to read in new datafile: ' + str(dt.datetime.now() - time_temp) )             
                st += read(file_names[file_counter])
                st.merge(fill_value='interpolate')#method=0) # Merge the new and old day volumes

                st.trim(starttime=t[i] - pp['coarse_duration'], endtime=st[0].stats.endtime ) # trim off the part of the merged stream that's already been processed.

            if data_source!='local':
                # Read in from FDSN client
                st = fdsn_client.get_waveforms(network=network, station=station, location='',
                                   channel=channel, starttime=t[i]-pp['coarse_duration'], endtime=t[i]+86400)
            
            
        except Exception as e:
            print(e)
            break # If there is no data to load, break out of the while loop 
                  #    and go to the next time step.

        # This next break loop short-circuits the while loop
#        if st[0].stats.endtime - pp['coarse_duration'] < tr_trim.stats.endtime:
#            break # If we were not able to load enough new data to get us
#                     #   through to the next iteration of the for loop (of t), 
#                     #   then force that next iteration.
        

#        print('about to remove response: ' + str(dt.datetime.now() - time_temp) )            
        st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)
#        print('response now removed: ' + str(dt.datetime.now() - time_temp) )             

#        print(st_IC)
#        IC_counter = 0 # Reset the IC_counter so that new st_IC will be created soon
        tr_trim = st_IC[0].copy() # copy instrument corrected trace
        # the minus small number and False nearest_sample make the tr include the first data point, but stop just shy of the last data point
#        tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - 0.0000001, nearest_sample=False)
        tr_trim.trim(t[i], t[i] + pp['coarse_duration'] - tr_trim.stats.delta)
#        print(tr_trim.stats.npts)

#    print(tr_trim)
#    print(st)
    
    # If tr_trim does not contain a full coarse_duration of samples, then skip
    #     the present coarse_duration and go onto the next iteration of the for
    #     loop.
    if i>46:
        flag_up = True
        
    if tr_trim.stats.npts < pp['coarse_duration'] * int(tr_trim.stats.sampling_rate) * 1:
        print('Incomplete coarse_duration at ' + UTCDateTime.strftime(t[i], "%d %b %y %H:%M") + ': Skipping')

#        if flag_up:
#            sys.exit()
            
#        print(tr_trim.stats.npts)    
        continue
    
#    print('Calculating median spectra for ' + UTCDateTime.strftime(t[i], "%d/%m/%y %H:%M"))
    # pass the seismic data with coarse_duration to the helper function for
    #     calculation of the median spectra.
#    print('Sending for get_med_spectra_v1: ' + str(dt.datetime.now() - time_temp) )             
    freqs, Pdb, Fs_old = get_med_spectra_v1.med_spec(tr_trim, pp, Fs_old)
    Pdb_array[:len(freqs),i] = Pdb[:len(freqs)] # Save the median spectra into an array
    

# At the end of the big for loop:
print('\n\n' + '===========================================')
print(station + ' run complete: ' + '{:%b %d, %Y, %H:%M}'.format(dt.datetime.now()))
print('Elapsed time: ' + str(dt.datetime.now() - run_start_time))
print('===========================================' + '\n\n')

# %% Pickle the output of the big runs

# Saving the objects:
with open(out_dir + 'mp' + network + '_' + station + '.pickle', 'wb') as f:  # Python 3: open(..., 'wb')
    pickle.dump([t, t_dt64, freqs, Pdb_array, pp, network, station, run_start_time], f)

