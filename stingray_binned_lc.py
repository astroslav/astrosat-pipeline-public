#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 09:33:36 2020

@author: mario
"""

import numpy as np
import matplotlib.pyplot as plt
from stingray_import import import_data_and_make_lc
from marioplots import marioplots

# Stylize the plots with my settings
marioplots()

# Make a lightcurve to make sure the lightcurve object was created properly. Here 
# I use a time resolution of 100 ms
lc = import_data_and_make_lc(np.array([3.0,80.0]), 0.1)

# Making a plot of the lightcurve 
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9,9), sharex=True, sharey=True)
fig.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.97, hspace=0)
ax1.plot(lc.time, lc.counts, 'k.')
ax1.set_ylabel("Time (s)")
ax1.set_ylabel("Counts (#)")


# Now let's create bins of 1000 counts each.
counts_per_bin = int(input('How many counts per bin would you like? (Default = 100000) ') or 100000)
print('Creating bins of %.0f counts each.' % (counts_per_bin))

# Create a "bin_size" with respect to how long (in seconds) it takes to reach X 
# counts on average. Then create a "bin_arr" so that we can cut the lightcurve 
# into chunks. 
bin_size = int(round(counts_per_bin/np.mean(lc.counts)))
bin_arr = np.arange(0, len(lc.time), bin_size, dtype='int')

lc_bins = []
for i in range(len(bin_arr)-1):
    lc_bins.append(lc[bin_arr[i]:bin_arr[i+1]])

# Plot the entire binned matrix to compare to the full lightcurve from earlier.
# (Mental check)
for lc_bin in lc_bins:
    ax2.plot(lc_bin.time, lc_bin.counts, '.')

ax1.set_ylabel("Counts (#)")
plt.show()
