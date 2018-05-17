# med_spec
Calculates the median spectrograms of waveforms from a seismic deployment.  Method follows that described in *Bartholomaus et al., 2015, GRL*.

Several key functions and a couple of helpers.
The core is **med_spec_loop_v2.py**:
    This file loops over a long duration seismic deployment, calling **get_med_spectra_v1.py** every time there's a snippet of data to calculate the median spectra of.  The output, which is the ingredients of a median spectrogram, includes time, frequency and power.  As of March 1, 2018, get_med_spectra_v1.py satisfies Parseval's theorem.  Previous versions of data run with this function reported power that was too low due to windowing that reduced amplitudes by an unaccounted for scalar.
    The file med_spec_loop_v2.py corrects the data using the instrument response and requires use of an appropriate RESP or station xml file, which is assumed located in resp_dir.

**plot_med_spec.py**:
    Read in the saved output from med_spec_loop_v2.py and create spectrograms of the processed data.
    
The next file is **load_all_GHT.py**:
    This script reads in a median spectrogram and then plots a time series of GHT amplitude by integrating the PSD over a specified range of frequencies.
