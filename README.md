# med_spec
Calculates the median spectrograms of waveforms from a seismic deployment.  Method follows that described in *Bartholomaus et al., 2015, GRL*.

Several key functions and a couple of helpers.
The core is **med_spec_loop_v2.py**:
    This file loops over a long duration seismic deployment, calling **get_med_spectra_v1.py** every time there's a snippet of data to calculate the median spectra of.  The output, which is the ingredients of a median spectrogram, includes time, frequency and power.
    This file corrects the data using the instrument response and requires use of an appropriate RESP file located in resp_dir

**plot_med_spec.py**:
    Read in the saved output from med_spec_loop_v2.py and create spectrograms of the processed data.
    
The next file is **load_all_GHT.py**:
    This script reads in a median spectrogram and then plots a time series of GHT amplitude by integrating the PSD over a specified range of frequencies.
