#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 09:33:36 2020

@author: mario
"""

import numpy as np
from get_gtis import get_gtis
from get_level2data import get_level2data
from EnLxLa_mask import EnLxLa_mask
from stingray.events import EventList

def import_data_and_make_lc(energy_range, dt, gti_fits='usergti.fits', 
                            level2_fits='level2.event.fits', gti_index=None):
    
    # Import level2data (and gtis) and split up columns into their respective vars
    gtis, gtis_min = get_gtis(gti_fits)
    level2data, ncounts, mjdref = get_level2data(level2_fits)

    time     = level2data[:,0]
    channel  = level2data[:,1]
    layer    = level2data[:,2]
    laxpc_no = level2data[:,3]
    energy   = level2data[:,4]
    
    # Define a function to make lightcurve using the Stingray "evts.to_lc" function.
    # You can pass a "gti_index", which essentially uses only a single GTI. So, for 
    # example, if you wrote "gti_index = 1", then you would use only the second GTI 
    # to create the lightcurve object. Otherwise it uses all GTIs. 

    print('Creating mask.')
    mask = EnLxLa_mask(energy, laxpc_no, layer, energy_range=energy_range)
    
    print('Creating EventList object.')
    if gti_index != None:
        evts = EventList(time=time[mask], ncounts=ncounts, mjdref=mjdref, 
                         gti=np.array([gtis[gti_index]]), energy=energy[mask], 
                         pi=channel[mask])
    else:
        evts = EventList(time=time[mask], ncounts=ncounts, mjdref=mjdref, 
                         gti=gtis, energy=energy[mask], pi=channel[mask])
                         
    print('Creating Lightcurve object.')
    lc = evts.to_lc(dt)
    
    return lc
