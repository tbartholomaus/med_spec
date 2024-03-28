# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 19:04:07 2016

@author: tbartholomaus

#==============================================================================
# % This is the core of the Calc_periodogram scripts.  It reads in the data,
# % calculates the DFTs and creates PSDs
# %
# % Modified/taken out of the core scripts on Feb. 4, 2014
# % Tim Bartholomaus
#==============================================================================

Ported to Python/obspy, from matlab, starting on Feb 24, 2016
Completed March 27, 2017

This script is called by med_spec_loop_v1.py
It is passed a sample of seismic data with coarse_duration, calculates the
median spectra of many subsamples of seismic data with fine_duration, and then
returns the median result to med_spec_loop for further analysis and plotting.

This implements the methods described in Bartholomaus et al., 2015, in GRL.

_v0.py had been unchanged at least since October 3, 2017.

_v1.py Version updated Nov. 29, 2017 to assume that waveform data has already
    been corrected for its instrument response.
    
"""


import numpy as np
from scipy import signal


def med_spec(tr, pp, Fs_old):
    # global WIND, NFFT, idx # on Mar 28, 2024, TCB comments these- don't use globals!

#%%   
#        tr = irisFetch.Traces(netw,station,'*',chan, ...
#            datestr(t(i), 'yyyy-mm-dd HH:MM:SS'), datestr(t(i)+coarse_duration/86400, 'yyyy-mm-dd HH:MM:SS'));
#%         [datestr([tr.startTime]), datestr([tr.endTime])]
#        if length(tr)~=1
#            continue
#        end
    Fs = int(tr.stats.sampling_rate)
#        if rem(Fs,1) ~= 0 % Check to see if the sample rate is an integer
#            fprintf(1,'\n');
#            disp('Sampling frequency (Fs) is not an integer     ')
#            disp('Progress:      ');
#            continue
#        end
#                if ~strcmp(tr.sensitivityUnits, 'M/S')
#            error('Sensitivity units are not in m/s.  Modify function appropriately and re-run.')
#        end
#        data = tr.data/tr.sensitivity;
    data = tr.data# * 10**-9
#        coarse_dur_retrieved(i) = tr.endTime - tr.startTime;
#    end

#%%
# START PROCESSING THE DATA
#if Fs*pp['course_duration'] > length(data)+Fs % Skip incomplete hours
#    continue
#end
    
    if Fs != Fs_old: # % Only recalculate these metrics when necessary.  Most of the time they can stay the same.
        print('>>> Recalulating idx and WIND <<<')        
        L = int(Fs * pp['fine_duration']) # Length of fine_window data vector
        NFFT = int(2**np.ceil(np.log2(L))) # Next power of 2 from length of fine data

        id1 = np.arange(0, L) # [0:L-1]';
        start_offset = np.arange(0, Fs*pp['coarse_duration']-L, L*(1-pp['fine_overlap']))
        start_offset = start_offset.astype('int')        
        idx = id1 + start_offset[:, None]

        Fs_old = Fs
    
    DATA = data[idx] # % turn the coarse-window data into many fine-window snippets of length L
#    print('DATA complete')  


#    print('DATA ready for periodogram')
    freqs, Pxx = signal.periodogram(DATA, fs=Fs, window='hann', nfft=NFFT,
                                    detrend='linear', return_onesided=True,
                                    scaling='density', axis=1)
    # With simple sin wave, TCB confirms Mar 1, 2018 
    #    that np.sum(DATA**2)/L = np.sum(Pxx) * delta_freq, for unwindowed data
    # Tim further confirms that when the data is windowed in the signal.periodogram
    #    call, that the powers are re-scaled so that Parseval's theorem holds.
    #    but that if the data is hanning windowed first, then sig.per has no 
    #    idea (of course), and the powers are reduced by a factor of 2.67
    #    On Mar 1, 2018, Tim changes the script so that windowing is done in
    #    the call to sig.per, and not externally.
    Pdb = 10 * np.log10(np.median(Pxx,0))

#%%
#==============================================================================
# plt.plot(data[idx[100:103,:]].T)
# plt.show()
# plt.plot(DATA[100:103,:].T)
# plt.show()
# #%%
# plt.plot(freqs, 10* np.log10(Pxx[100:103,:].T))
# plt.plot(freqs, Pdb, 'k')
# plt.xscale('log')
# plt.xlim([0.1, 100])
# plt.ylim([-220, -140])
#==============================================================================
#%%
#    print('Ready to return results from get_med_spectra')
    return freqs, Pdb, Fs
# %%


#==============================================================================
# fig, ax = plt.subplots()
# 
# 
# ax.plot(freqs, 10*np.log10(Pxx.T) )
# ax.plot(freqs, Pdb, 'k')
# 
# ax.set_ylim([-190, -150])
# ax.set_xlim([0, 50])
#==============================================================================


#==============================================================================
#                 
#     freqs, Pxx = mlab.psd(DATA, NFFT=L, Fs=Fs, detrend='none', 
#         noverlap=0, pad_to=NFFT, 
#         sides='default', scale_by_freq=True)    
#     
#     
#     DATA = detrend(DATA); % demeans AND detrends DATA with this one command
#     DATA = WIND.*DATA; % Apply window to DATA, prior to FFT, to reduce sidelobe leakage.
#     
#     MAG = abs(fft(DATA,NFFT));
#     MAG = MAG(1 : NFFT/2 +1, :);
#     POWR = (1/(L*NFFT)) * MAG.^2;
#     POWR(2:end-1, :) = 2 * POWR(2:end-1, :);  % This is the PSD.  Matrix of PSDs for all the short window noise.
#     
#==============================================================================
