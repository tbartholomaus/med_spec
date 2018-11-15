# med_spec
Calculates the median spectrograms of waveforms from a seismic deployment.  Method follows that described in *Bartholomaus et al., 2015, GRL*.

This repository consists of several key functions and a couple of helpers to calculate the median power at a range of frequencies for seismic timeseries at single stations.  These median spectra are less sensitive to brief, high power transient events such as earthquakes and icequakes, and thus are appropriate for use in noisy seismic environments, like glaciers.

The core of these scripts begins with **med_spec_loop_v3.py**. This file loops over a long duration seismic deployment and calls **get_med_spectra_v1.py** every time there's a snippet of data to calculate the median spectra of.  The output of get_med_spectra_v1.py, which contains the ingredients of a median spectrogram, includes time, frequency and power.  As of March 1, 2018, get_med_spectra_v1.py satisfies Parseval's theorem.  Previous versions of data run with this function reported power that was too low due to windowing that reduced amplitudes by an unaccounted-for scalar.
    The file med_spec_loop_v3.py corrects the data using the instrument response and requires use of an appropriate RESP or station xml file.
    
As of version 3 of med_spec_loop (15 November 2018), neither med_spec_loop_v3.py nor get_med_spectra_v1.py should require modification to run in a variety of configurations and with a variety of datasets. Key parameters and paths are defined in a parameter file, which includes all of these.  An example parameter file is included in this repository as **med_spec.par**. 

---

At the command prompt, the successful execution of the med_spec algorithm looks like:
```
python ./med_spec_loop_v3.py med_spec.par  
```

A more feature-rich execution of this code could be
```
nohup python -u ./med_spec_loop_v2.py med_spec.par LM BBWL > BBWL.log &  
```
which is expected to run on a server with nohup preventing killing the process if the server is disconnected from.  The `-u` flag forces the text output of the python script to stdout, so that it can then be redirected to a log file (`> BBWL.log`).  Here network and station codes from the par file are also overridden by `LM` and `BBWL`.  If present, these network, station, and channel codes (in that order) take precedence over whatever is in the par file.

---

Ancillary scripts include the following: 
**plot_med_spec.py**:
    Read in the saved output from med_spec_loop_v2.py and create spectrograms of the processed data.
    
The next file is **load_all_GHT.py**:
    This script reads in a median spectrogram and then plots a time series of GHT amplitude by integrating the PSD over a specified range of frequencies.
