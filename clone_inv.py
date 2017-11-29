#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  4 11:06:15 2017

@author: timb
"""


import obspy
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
#from obspy.clients.nrl import NRL

def clone_inv(inv, net_name, sta_name):
    
    net = Network(
        # This is the network code according to the SEED standard.
        code=net_name,
        # A list of stations. We'll add one later.
        stations=[],
#        description="A test stations.",
        # Start-and end dates are optional.
#        start_date=obspy.UTCDateTime(2016, 1, 2))
        )
    
    sta = Station(
        # This is the station code according to the SEED standard.
        code=sta_name,
        latitude=  inv[0][0].latitude,
        longitude= inv[0][0].longitude,
        elevation= inv[0][0].elevation,
        creation_date=obspy.UTCDateTime(2016, 1, 2),
        site=Site(name="First station"))
    
    cha = Channel(
        # This is the channel code according to the SEED standard.
        code="HHZ",
        # This is the location code according to the SEED standard.
        location_code="",
        # Note that these coordinates can differ from the station coordinates.
        latitude=  inv[0][0][0].latitude,
        longitude= inv[0][0][0].longitude,
        elevation= inv[0][0][0].elevation,
        depth=     inv[0][0][0].depth,
#        azimuth=0.0,
#        dip=-90.0,
        sample_rate= inv[0][0][0].sample_rate)
    
    
    
    # Now tie it all together.
    cha.response = inv[0][0][0].response #response
    sta.channels.append(cha)
    net.stations.append(sta)
    inv.networks.append(net)
    
    return inv