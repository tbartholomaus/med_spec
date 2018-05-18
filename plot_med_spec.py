#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_med_spec.py

Script reads in the saved (pickled) data from a run that calculates median
spectra, and then plots them.  The core of the script is a function spec_plt,
which includes the core plotting routines.  After this function is defined, a
list of stations (i.e., saved data sets) is defined.  Each station is loaded
into memory, and a series of plots are created according to specified frequency
and date limits.

Created on Fri Dec  1 13:06:40 2017

@author: timb
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.dates as mdates

power_lims = [-200, -120] # Lower and upper bounds to use for color scale in
    # plotting the spectrograms

# First, define a function call that will do the plotting.
def spec_plt(plt_type, freq_lims, date_lims):
    mask_val1 = Pdb_array<=-300
    mask_val2 = np.isinf(Pdb_array)
    
    # Create another mask array that will rely on unusually low power in the microseism band
    mask_val3 = mask_val1.copy()
    freq_ind = np.where(freqs>.2)[0][0] # find the index of the first frequency greater than 0.2 Hz
    mask_ind3 = np.where(Pdb_array[freq_ind,:]<-170)[0] # If the power within the dominant microseism band (~0.2 Hz) is less than -170 dB, then be suspicious of the data.
    mask_val3[:, mask_ind3] = True
    
    Pdb_array_mask = np.ma.masked_where(np.any([mask_val1, mask_val2, mask_val3], axis=0), Pdb_array)
    
    fig, ax = plt.subplots()#figsize=(8, 4))
    qm = ax.pcolormesh(t_dt64, freqs[1:], Pdb_array_mask[1:,:], cmap='YlOrRd')#, extent = [0, len(file_names), freqs[1], freqs[freq_nums-1]])
    ax.set_yscale('log')
#    ax.set_ylim()
    ax.set_ylim(freq_lims)
    ax.set_ylabel('Frequency (Hz)')
    
    # Set the date limits for the plot, to make each station's plots consistent
    ax.set_xlim(date_lims)
    
    # Format the xaxis of the pcolormesh to be dates
    #ax.xaxis.set_major_locator(mdates.AutoDateLocator())
#    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
#    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b, %H:%M'))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%d %b'))
    fig.autofmt_xdate()
    
    qm.set_clim(vmin=power_lims[0], vmax=power_lims[1])
    cb = plt.colorbar(qm, ticks=np.arange(power_lims[0], power_lims[1], 10))  # typically -200 to -150
    cb.set_label('Power (dB, rel. 1 (m/s)^2/Hz)')
    plt.title(station)
    
    # #%%
    plt.savefig('output_figs/Spec_' + plt_type + '_' + station + '', dpi=150) # _ght, _fld
    plt.show()

# %% Read in and plot each station.  In this case, output data is saved
#       according to the network code name, and the station name:
#stations = ['BBWU', 'BBEU', 'BBGU', 'BBWL', 'BBEL', 'BBGL']
stations = ['XX_BBGL']#['XH_FX01', 'XH_FX03', 'XH_FX06', 'XH_FX10', 'XH_FX11', 'XH_FX12']
stations = ['XH_FX11']#['XH_FX01', 'XH_FX03', 'XH_FX06', 'XH_FX10', 'XH_FX11', 'XH_FX12']
#['7E_DL1', '7E_S1', '7E_S2', '7E_S4', '7E_S5', '7E_S6']#, 'BBGL']
#['XF_BOOM', 'XF_DOST', 'XF_GRAP']#
stations = [#'ZQ_ETIP', 'ZQ_TWLV', 'ZQ_RTBD', 'ZQ_TAKN', 'ZQ_TAKE', 
            'ZQ_TAKW', 
            'ZQ_TAKC', 'ZQ_HITW', 'ZQ_GAGA', 'ZQ_GUGU', 'ZQ_GIW1', 'ZQ_GIW2']

stations = ['XX_BBGL', 'XX_BBGU']
stations = ['ZQ_GAGA', 'ZQ_GUGU']

#
for station in stations:
#    station = 'BBEU'
    #filename = 'mpBBEL_test2.pickle'
    filename = 'output_results/mp' + station + '.pickle'
    with open(filename, 'rb') as f:  # Python 3: open(..., 'rb')
        t, t_dt64, freqs, Pdb_array, pp, network, station, run_start_time = pickle.load(f, encoding='latin1')
#        t, t_dt64, freqs, Pdb_array, pp, data_dir, station 
    t_start = t_dt64[0]
    t_end = t_dt64[-1]

#    freqs = np.linspace(0,100, 2049)
    #%% PLOT COMMANDS:
    # Each block of code takes a plt_type (which is a name for that style of
    #   plot), the frequency limits for the spectrogram, and the date limits
    #   for each type of spectrogram
    
    #date_lims = np.array(['2017-10-02', '2017-10-07'], dtype='datetime64' )
    
    plt_type = 'comp'
    freq_lims = [0.1, 90] # Hz  [0.5, 80] 
    date_lims = np.array([t_start, t_end], dtype='datetime64' ) # 
#    date_lims = np.array(['2016-04-08', '2016-04-13'], dtype='datetime64' )
    spec_plt(plt_type, freq_lims, date_lims)
    
#    plt_type = 'early'
#    freq_lims = [0.1, 90] # [0.5, 80]
#    date_lims = np.array(['2016-03-23', '2016-04-15'], dtype='datetime64' )
##    date_lims = np.array(['2016-04-08', '2016-04-13'], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)
#    
#    plt_type = 'ght'
#    freq_lims = [.8, 15] # [0.5, 80]
#    date_lims = np.array([t_start, t_end], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)
    
    
#
#    
#    plt_type = 'fld'
#    freq_lims = [0.5, 220]
#    date_lims = np.array(['2017-07-03', '2017-07-09'], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)
#    
#    plt_type = 'ght'
#    freq_lims = [1.5, 30] # [0.5, 80]
#    date_lims = np.array(['2017-06-27', '2017-09-27'], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)
#
#    plt_type = 'spd'
#    freq_lims = [0.1, 100] # [0.5, 80]
#    date_lims = np.array(['2017-07-22', '2017-07-31'], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)
    
#    plt_type = 'glac_comp'
#    freq_lims = [0.2, 80] # [0.5, 80]
#    date_lims = np.array(['2017-06-30', '2017-08-15'], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)

#    plt_type = 'glac_fld'
#    freq_lims = [0.2, 80] # [0.5, 80]
#    date_lims = np.array(['2017-07-05', '2017-07-09'], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)
#    
#    plt_type = 'glac_spd'
#    freq_lims = [0.2, 80] # [0.5, 80]
#    date_lims = np.array(['2017-07-26', '2017-07-30'], dtype='datetime64' )
#    spec_plt(plt_type, freq_lims, date_lims)    