#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 10:59:58 2019

@author: mario
"""

import matplotlib.pyplot as plt
from get_gtis import get_gtis
from get_level2data import get_level2data
from EnLxLa_mask import EnLxLa_mask
from stingray.events import EventList

print('Importing level2.event.fits and usergti.fits')
gtis, gtis_min = get_gtis('usergti.fits')
level2data, ncounts, mjdref = get_level2data('level2.event.fits')

time     = level2data[:,0]
channel  = level2data[:,1]
layer    = level2data[:,2]
laxpc_no = level2data[:,3]
energy   = level2data[:,4]

print('Mask created with default values: 3.0-80.0 keV, laxpc 10 & 20, all layers')

mask = EnLxLa_mask(energy, laxpc_no, layer)
evts = EventList(time=time[mask], ncounts=ncounts, mjdref=mjdref, 
                 gti=gtis, energy=energy[mask], pi=channel[mask])

dt = float(input('Input the time resolution you\'d like for the lightcurve: (default = 0.1s) ')) or 0.1

lc = evts.to_lc(dt)

plt.figure(1)
plt.clf()
plt.errorbar(lc.time,lc.countrate,yerr=lc.countrate_err,fmt='k.',ms=2,
             elinewidth=0.2,label='source+bkg')
plt.xlabel('Time (s)',size=14)
plt.ylabel(r'Count rate (s$^{-1}$)',size=14)
plt.legend(loc='best',fontsize=16)
plt.tight_layout()
plt.show()
