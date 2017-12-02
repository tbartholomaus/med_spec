#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 13:06:40 2017

@author: timb
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle

#%% Plot the output of the big runs as median spectrograms

def spec_plt(plttype, freqlims, datelims):
    mask_val1 = Pdb_array<=-300
    mask_val2 = np.isinf(Pdb_array)
    Pdb_array_mask = np.ma.masked_where(np.any([mask_val1, mask_val2], axis=0), Pdb_array)
    
#    t_datenum = UTC2dn(t) # Convert a np.array of obspy UTCDateTimes into datenums for the purpose of plotting
    
#    plt.imshow(np.log10(Pxx_vals[0:-2,]), extent = [0, len(file_names), freqs[1], freqs[freq_nums-1]])
    fig, ax = plt.subplots()#figsize=(8, 4))
    qm = ax.pcolormesh(t_dt64, freqs, Pdb_array_mask, cmap='YlOrRd')#, extent = [0, len(file_names), freqs[1], freqs[freq_nums-1]])
    ax.set_yscale('log')
    ax.set_ylim()
    ax.set_ylim(freqlims)
    ax.set_ylabel('Frequency (Hz)')
    
    # Set the date limits for the plot, to make each station's plots consistent
    ax.set_xlim(datelims)
    #ax.set_xlim(mdates.date2num([dt.date(2017, 7, 3), dt.date(2017, 7, 8)]))
    #ax.set_xlim(mdates.date2num([dt.datetime(2017, 10, 4, 23,0,0), dt.datetime(2017, 10, 7, 1, 0, 0)]))
    
    # Format the xaxis of the pcolormesh to be dates
    #ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    fig.autofmt_xdate()
    
    qm.set_clim(vmin=-200, vmax=-150)
    cb = plt.colorbar(qm, ticks=np.arange(-200,-150, 10))
    cb.set_label('Power (dB, rel. 1 (m/s)^2/Hz)')
    plt.title(station)
    
    # #%%
#    plt.savefig('output_figs/Spec_' + plttype + '_' + station + '', dpi=150) # _ght, _fld
    plt.show()

# %% Getting back the objects:
station = 'BBEL'
filename = 'mpBBEL_test2.pickle'
#filename = 'output_results/mp' + station + '.pickle'
with open(filename, 'rb') as f:  # Python 3: open(..., 'rb')
    t, t_dt64, freqs, Pdb_array, pp, data_dir, station = pickle.load(f, encoding='latin1')

#%%
#datelims = np.array(['2017-10-02', '2017-10-07'], dtype='datetime64' )


plttype = 'comp'
freqlims = [0.1, 250] # [0.5, 80]
datelims = np.array(['2017-06-27', '2017-09-27'], dtype='datetime64' )
spec_plt(plttype, freqlims, datelims)

plttype = 'fld'
freqlims = [0.5, 80]
datelims = np.array(['2017-07-02', '2017-07-08'], dtype='datetime64' )
spec_plt(plttype, freqlims, datelims)

plttype = 'ght'
freqlims = [1.5, 20] # [0.5, 80]
datelims = np.array(['2017-06-27', '2017-09-27'], dtype='datetime64' )
spec_plt(plttype, freqlims, datelims)

