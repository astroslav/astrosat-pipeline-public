#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 12:54:19 2020

@author: mario
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from stingray_import import import_data_and_make_lc
from stingray import Powerspectrum, AveragedPowerspectrum
from marioplots import marioplots

# Stylize the plots with my settings
marioplots()

# Make a lightcurve to make sure the lightcurve object was created properly. Here 
# I use a time resolution of 100 ms
lc = import_data_and_make_lc(np.array([3.0,80.0]), 0.1)

# Using the lc._truncate_by_time function to check how it works. Plot result
plt.plot(lc.time, lc.counts, 'k.')
lc2 = lc._truncate_by_time(lc.gti[2][0], lc.gti[2][1])
plt.plot(lc2.time, lc2.counts, 'ro')
plt.show()

# Test Powerspectrum function from Stingray
ps_list = []
lc_list = []
for i in range(len(lc.gti)):
    lc_split = lc._truncate_by_time(lc.gti[i][0], lc.gti[i][1])
    lc_split.gti = np.array([[lc.gti[i][0], lc.gti[i][1]]])
    lc_list.append(lc_split)
    ps_list.append(Powerspectrum(lc_split))

colors = plt.cm.jet(np.linspace(0,1,len(ps_list)))
fig, ax1 = plt.subplots(1,1,figsize=(9,6))#, sharex=True)
i=0
for ps in ps_list:
    ax1.plot(ps.freq, ps.power, lw=1, color=colors[i])
    i+=1
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Power (raw)")
ax1.set_yscale('log')
plt.show()

# Test AveragedPowerspectrum function from Stingray
avg_ps = AveragedPowerspectrum(lc, 100)
fig, ax1 = plt.subplots(1,1,figsize=(9,6))
ax1.plot(avg_ps.freq, avg_ps.power, lw=1, color='b')
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Power (raw)")
ax1.set_yscale('log')
plt.show()

#%%
avg_ps_list = []
for i in range(len(lc.gti)):
    lc_split = lc._truncate_by_time(lc.gti[i][0], lc.gti[i][1])
    lc_split.gti = np.array([[lc.gti[i][0], lc.gti[i][1]]])
    avg_ps_list.append(AveragedPowerspectrum(lc_split, 10))

colors = plt.cm.jet(np.linspace(0,1,len(ps_list)))
fig, ax1 = plt.subplots(1,1,figsize=(9,6))#, sharex=True)
i=0
for avg_ps in avg_ps_list:
    ax1.plot(avg_ps.freq, avg_ps.power, lw=1, color=colors[i])
    i+=1
ax1.set_xlabel("Frequency (Hz)")
ax1.set_ylabel("Power (raw)")
ax1.set_yscale('log')
plt.show()
