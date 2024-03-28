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

Created on Fri Dec  1 13:06:40 2017 by @author: timb

Updated most recently in February 2024 by Yoram Terleth. - cleaned up the script, added storage of spectrogram as a .mat file. 
"""


import numpy as np
import matplotlib.pyplot as plt
import pickle
import matplotlib.dates as mdates
from scipy.io import savemat

# Define a function for plotting spectrograms
def spec_plt(plt_type, freq_lims, date_lims, power_lims,path_to_pickle_files, save_as_mat):
    mask_val1 = Pdb_array <= -300
    mask_val2 = np.isinf(Pdb_array)

    # Create a mask array based on unusually low power in the microseism band
    mask_val3 = mask_val1.copy()
    freq_ind = np.where(freqs > 0.2)[0][0]  # Find the index of the first frequency greater than 0.2 Hz
    mask_ind3 = np.where(Pdb_array[freq_ind, :] < -170)[0]  # If power within the microseism band is less than -170 dB, be suspicious of the data
    mask_val3[:, mask_ind3] = True

    Pdb_array_mask = Pdb_array  # Use a mask array to filter out suspicious data

    fig, ax = plt.subplots(figsize=(20, 5))
    qm = ax.pcolormesh(t_dt64, freqs[1:], Pdb_array_mask[1:, :], cmap='YlOrRd')  # Plot the spectrogram
    ax.set_yscale('log')
    ax.set_ylim(freq_lims)
    ax.set_ylabel('Frequency (Hz)')
    ax.set_xlim(date_lims)  # Set date limits for the plot
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))  # Format x-axis as dates
    qm.set_clim(vmin=power_lims[0], vmax=power_lims[1])  # Set color scale limits
    cb = plt.colorbar(qm, ticks=np.arange(power_lims[0], power_lims[1], 10))  # Add color bar
    cb.set_label('Power (dB, rel. 1 $m^{2}s^{-1}$)')
    plt.title(station)
    plt.savefig(path_to_pickle_files + 'output_figs/Spec_' + plt_type + '_' + station + '', dpi=700)  # Save plot as image
    plt.show()

    # Save spectrogram data as .mat file
    if save_as_mat:
        savemat(path_to_pickle_files + '/spectrogram' + station + '.mat', {'a': [t_dt64, freqs[1:], Pdb_array_mask[1:, :]]}, appendmat=False)
    



# define path to the stored pickle files, produced by med_spec_loop.py
path_to_pickle_files = '/data/stor/proj/Turner/med_spec_2024/output_files_pickles/'

# do you want to save the spectrogram as a matlab file? 
save_as_mat = True 

# define stations, can be as many as you want, will produce individual images for each station.  
stations = ['YG_SE7', 'YG_SW7']  # List of stations to process

# Loop over each station
for station in stations:
    # Load data
    filename = path_to_pickle_files + 'mp' + station + '.pickle'
    with open(filename, 'rb') as f:
        t, t_dt64, freqs, Pdb_array, pp, network, station, run_start_time = pickle.load(f, encoding='latin1')
    t_start = t_dt64[0]
    t_end = t_dt64[-1]

    # define spectrogram characteristics
    plt_type = 'comp'
    freq_lims = [0.5, 80]
    power_lims = [-185, -135]
    date_lims = np.array(['2020-09-15', '2023-09-01'], dtype='datetime64')

    # execute the function defined above
    spec_plt(plt_type, freq_lims, date_lims, power_lims, path_to_pickle_files, save_as_mat)


