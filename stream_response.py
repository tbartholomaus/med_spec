#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 16:58:54 2018

@author: tbartholomaus
"""

import obspy
from clone_inv import clone_inv

import numpy as np

network= 'XX'
station= 'SE47'
channel= 'HHZ'
pre_filt = (0.1, .2, 90, 100.)

#%%

# include full paths to data, with final forward slash: /
# script assumes miniseed files are organized in data_dir + station/
data_dir = '/Users/timb/Desktop/'

# include full path to resp file
resp_file = '/Users/timb/Desktop/Resp/XX.UI05.resp/XX.UI05.HHZ.resp'

#%timeit 
st = obspy.read(data_dir + station + '/SE47.XX..HHZ.2018.130')

#assdf
#%timeit 
inv = obspy.read_inventory(resp_file)
inv = clone_inv(inv, network, station)
inv = inv.select(channel=channel, station=station) # subset the inventory to just that which is necessary for script


#%timeit st_IC = st.copy().remove_response(inventory=inv, output="VEL", pre_filt=pre_filt)
%timeit fft_out = np.fft.fft(st[0].data)
