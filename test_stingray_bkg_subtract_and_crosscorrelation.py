#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 13:38:51 2020

@author: mario
"""

import numpy as np
import glob
import sys
import os
import matplotlib.pyplot as plt
from subprocess import call
from get_gtis import get_gtis
from get_level2data import get_level2data
from EnLxLa_mask import EnLxLa_mask
from stingray.events import EventList
from stingray.events import Lightcurve
from stingray.crosscorrelation import CrossCorrelation
from scipy.interpolate import interp1d


#%% ===========================================================================
# Import data (level2data (and gtis) and backlc)
# =============================================================================
gtis, gtis_min = get_gtis('usergti.fits')
level2data, ncounts, mjdref = get_level2data('level2.event.fits')

time     = level2data[:,0]
channel  = level2data[:,1]
layer    = level2data[:,2]
laxpc_no = level2data[:,3]
energy   = level2data[:,4]

def make_lc(energy_range, dt):
    print('Creating mask.')
    mask = EnLxLa_mask(energy, laxpc_no, layer, energy_range=energy_range)
    print('Creating EventList object.')
    evts = EventList(time=time[mask], ncounts=ncounts, mjdref=mjdref, gti=gtis,
                     energy=energy[mask], pi=channel[mask])
    print('Creating Lightcurve object.')
    lc = evts.to_lc(dt)
    return lc

dt = 0.1

lc_3_10  = make_lc(np.array([3.0,10.0]),dt)
lc_10_30 = make_lc(np.array([10.0,30.0]),dt)
lc_30_55 = make_lc(np.array([30.0,55.0]),dt)
lc_55_80 = make_lc(np.array([55.0,80.0]),dt)


#%% ===========================================================================
# FUNCTIONS
# =============================================================================
# Find GTI's from backlc
def create_gti_from_time_arr(time_array,dt,addpad=False):
    gti_l_edges = [time_array[0]]
    gti_r_edges = []
    for i in range(len(time_array)-1):
        if (time_array[i+1] - time_array[i]) > dt*1.001:
            gti_l_edges.append(time_array[i+1])
            gti_r_edges.append(time_array[i])

    gti_r_edges.append(time_array[-1])
    gti = np.array([gti_l_edges,gti_r_edges]).T

    if addpad:
        gti[:,0] = gti[:,0] - 0.5*dt
        gti[:,1] = gti[:,1] + 0.5*dt

    return gti


# Function that interpolates in between given array
def self_interpolate(old_array,old_dt,new_dt):
    new_array = []
    for i in range(len(old_array)-1):
        # Multiply old_dt by 1.001 to get around floating point errors (skips)
        if (old_array[i+1] - old_array[i]) <= old_dt*1.001:
            new_array.append(np.linspace(old_array[i],old_array[i+1],
                                         int(round((old_dt/new_dt)+1))))

    return np.unique(np.asarray(new_array).flatten())


# Function that adds poisson noise to a Lightcurve object
def add_poisson_to_lc(lc):
    poisson_counts = np.random.poisson(lc.counts)
    poisson_countrate = poisson_counts/lc.dt
    return Lightcurve(time=lc.time, counts=poisson_countrate, gti=lc.gti, mjdref=mjdref,
                      input_counts=False)


# Generate interpolation function and interpolate time and countrates
def interp_time_series(lc,new_dt,mental_check=False):
    interp_func = interp1d(lc.time, lc.countrate, kind='cubic')#, bounds_error=False,
                           #fill_value='extrapolate')

    # Self-interpolate time
    print('Step 1/2: Self-interpolating time.')
    lc_interp_time = self_interpolate(lc.time,lc.dt,new_dt)

    # Interpolate count rate
    print('Step 2/2: Interpolating countrate.')
    lc_interp_countrate = interp_func(lc_interp_time)

    if mental_check:
        print('Mental check.')
        lc_interp_orig_time = self_interpolate(lc.time,lc.dt,lc.dt)

        if np.equal(lc_interp_orig_time,lc.time).all() == False:
            sys.exit('self_interpolate function is doing something wacky.')
        else:
            print('Returning Lightcurve object.')
            return Lightcurve(time=lc_interp_time, counts=lc_interp_countrate, gti=lc.gti,
                              mjdref=mjdref, input_counts=False)
    else:
        return Lightcurve(time=lc_interp_time, counts=lc_interp_countrate, gti=lc.gti,
                          mjdref=mjdref, input_counts=False)


# Subtracts background from Lightcurve object and returns it
def subtract_background(lc,backlc):
    interp_func = interp1d(backlc.time, backlc.counts, kind='linear')

    # if lc time is outside background time, ignore
    time_bounds = (backlc.time[0] <= lc.time) * (lc.time <= backlc.time[-1])

    interp_counts = interp_func(lc.time[time_bounds])
    subtracted_counts = lc.counts[time_bounds] - interp_counts

    return Lightcurve(time=lc.time[time_bounds], counts=subtracted_counts, mjdref=mjdref, gti=lc.gti)


#%% ===========================================================================
# Import background data
# =============================================================================

fnames = sorted(glob.glob('combined_backlc*.txt'))

backlc = []
for i in np.arange(0,len(glob.glob('combined_backlc*.txt')),1):
    print('data',i)
    backlc.append(np.loadtxt(fnames[i]))
backlc = np.asarray(backlc)

backlc_gti = create_gti_from_time_arr(backlc[0,:,0],0.1,addpad=True)

backlc_3_10 = Lightcurve(time=backlc[2,:,0], counts=backlc[2,:,1], err=backlc[2,:,2],
                         gti=backlc_gti, mjdref=mjdref, input_counts=False)
backlc_10_30 = Lightcurve(time=backlc[0,:,0], counts=backlc[0,:,1], err=backlc[0,:,2],
                          gti=backlc_gti, mjdref=mjdref, input_counts=False)
backlc_30_55 = Lightcurve(time=backlc[1,:,0], counts=backlc[1,:,1], err=backlc[1,:,2],
                          gti=backlc_gti, mjdref=mjdref, input_counts=False)
backlc_55_80 = Lightcurve(time=backlc[3,:,0], counts=backlc[3,:,1], err=backlc[3,:,2],
                          gti=backlc_gti, mjdref=mjdref, input_counts=False)


#%% ===========================================================================
# Working with Stingray on background data
# =============================================================================

def In_Po_Sub(lc,backlc,new_dt):
    lc_interp = interp_time_series(lc,new_dt)
    backlc_interp = interp_time_series(backlc,new_dt)
    backlc_poisson = add_poisson_to_lc(backlc_interp)
    lc_sub = subtract_background(lc_interp,backlc_poisson)
    return lc_sub

lc_3_10_sub100ms = In_Po_Sub(lc_3_10, backlc_3_10, 0.1)
lc_55_80_sub100ms = In_Po_Sub(lc_55_80, backlc_55_80, 0.1)

lc_3_10_sub10ms = In_Po_Sub(lc_3_10, backlc_3_10, 0.01)
lc_55_80_sub10ms = In_Po_Sub(lc_55_80, backlc_55_80, 0.01)

lc_3_10_sub = In_Po_Sub(lc_3_10, backlc_3_10, 0.001)
lc_55_80_sub = In_Po_Sub(lc_55_80, backlc_55_80, 0.001)

# lc_3_10_interp = interp_time_series(lc_3_10, 0.01)
# backlc_3_10_interp = interp_time_series(backlc_3_10, 0.01)
# backlc_3_10_poisson = add_poisson_to_lc(backlc_3_10_interp)
# lc_3_10_sub = subtract_background(lc_3_10_interp, backlc_3_10_poisson)

# lc_55_80_interp = interp_time_series(lc_55_80, 0.01)
# backlc_55_80_interp = interp_time_series(backlc_55_80, 0.01)
# backlc_55_80_poisson = add_poisson_to_lc(backlc_55_80_interp)
# lc_55_80_sub = subtract_background(lc_55_80_interp, backlc_55_80_poisson)

#%% diff poisson
# poisson_countrate = np.random.poisson(backlc_3_10_interp.countrate)
# backlc_3_10_poisson2 = Lightcurve(backlc_3_10_interp.time, poisson_countrate,
#                                   input_counts=False, gti=backlc_3_10_interp.gti)

#%%

def plot_Lightcurve(lc):
    plt.figure()
    plt.plot(lc.time, lc.counts, 'ko', ms=5)
    plt.xlabel('Time [s]', size=14)
    plt.ylabel('Counts', size=14)
    plt.show()

plot_Lightcurve(backlc_3_10_interp)
plot_Lightcurve(backlc_3_10_poisson)
plot_Lightcurve(lc_3_10_interp)
plot_Lightcurve(lc_3_10_sub)

#%%

def plot_CrossCorrelation(cr):
    time_bounds = (cr.time_lags > -100) * (cr.time_lags < 100)

    plt.figure()
    plt.plot(cr.time_lags[time_bounds], cr.corr[time_bounds], 'k.', ms=0.5)
    plt.xlabel('Time Lag (seconds)', size=14)
    plt.ylabel('Correlation', size=14)
    plt.show()



cr_lo_ex = CrossCorrelation(lc_3_10_sub,lc_55_80_sub)
cr_10ms = CrossCorrelation(lc_3_10_sub10ms,lc_55_80_sub10ms)
cr_100ms = CrossCorrelation(lc_3_10_sub100ms,lc_55_80_sub100ms)




plot_CrossCorrelation(cr_lo_ex)
plot_CrossCorrelation(cr_10ms)
plot_CrossCorrelation(cr_100ms)



#%% ===========================================================================
# Working with Stingray
# =============================================================================

cr_lo_me = CrossCorrelation(lc_sub_3_10,lc_sub_10_30)
cr_lo_hi = CrossCorrelation(lc_sub_3_10,lc_sub_30_55)
cr_lo_ex = CrossCorrelation(lc_sub_3_10,lc_sub_55_80)

cr_lo_me.plot(labels = ['Time Lag (seconds)', 'Correlation (3-10 keV vs. 10-30 keV)'])
cr_lo_hi.plot(labels = ['Time Lag (seconds)', 'Correlation (3-10 keV vs. 30-55 keV)'])
cr_lo_ex.plot(labels = ['Time Lag (seconds)', 'Correlation (3-10 keV vs. 55-80 keV)'])

#%% ----

#%%

lc_interp_3_10 = interp_time_series(lc_sub_3_10,0.001)
lc_interp_55_80 = interp_time_series(lc_sub_55_80,0.001)

#%%

cr_lo_ex = CrossCorrelation(lc_sub_3_10,lc_sub_55_80)
cr_lo_ex.plot(labels = ['Time Lag (seconds)', 'Correlation (3-10 keV vs. 55-80 keV)'])


# interp_func = interp1d(backlc_time, backlc_countrate, kind='cubic', bounds_error=False,
#                         fill_value='extrapolate')

# # Self-interpolate time and interpolate countrate (i.e. rebin lc) to 10 ms and 1 ms
# backlc_0_01_interp_time = self_interpolate(backlc_time,0.1,0.01)
# backlc_0_01_interp_countrate = interp_func(backlc_0_01_interp_time)

# backlc_0_001_interp_time = self_interpolate(backlc_time,0.1,0.001)
# backlc_0_001_interp_countrate = interp_func(backlc_0_001_interp_time)

# # Mental check (if we rebin at 0.1s, arrays should be equal)
# backlc_0_1_rebin_time = self_interpolate(backlc_time,0.1,0.1)
# if (backlc_0_1_rebin_time == backlc_time).all():
#     print('Time arrays match.')
# else:
#     quit('Something is wrong. Time arrays didn\'t match.')


#%% ===========================================================================
# Generate Lightcurve objects and add Poisson noise
# =============================================================================
# backlc_gti = create_gti_from_time_arr(backlc_time,0.1,addpad=True)

# backlc_0_1 = Lightcurve(time=backlc_time, counts=backlc_countrate, err=backlc_countrate_err,
#                         gti=backlc_gti, mjdref=mjdref, input_counts=False)

# backlc_0_01 = Lightcurve(time=backlc_0_01_interp_time, counts=backlc_0_01_interp_countrate,
#                          mjdref=mjdref, input_counts=False)

# backlc_0_001 = Lightcurve(time=backlc_0_001_interp_time, counts=backlc_0_001_interp_countrate,
#                           mjdref=mjdref, input_counts=False)

# # Add poisson noise to Lightcurves (wP i.e. with Poisson)
# add_poisson_to_lc(backlc_0_1)
# add_poisson_to_lc(backlc_0_01)
# add_poisson_to_lc(backlc_0_001)

# # Print some stats
# print('')
# print('Original mean and std:   %.4f +/- %.4f' % (np.mean(backlc_countrate),np.std(backlc_countrate)))
# print('Interp 10 ms:            %.4f +/- %.4f' % (np.mean(backlc_0_01_interp_countrate),np.std(backlc_0_01_interp_countrate)))
# print('Interp 1 ms:             %.4f +/- %.4f' % (np.mean(backlc_0_001_interp_countrate),np.std(backlc_0_001_interp_countrate)))
# print('Original w/ Poisson:     %.4f +/- %.4f' % (np.mean(backlc_0_1.countrate),np.std(backlc_0_1.countrate)))
# print('Interp 10 ms w/ Poisson: %.4f +/- %.4f' % (np.mean(backlc_0_01.countrate),np.std(backlc_0_01.countrate)))
# print('Interp 1 ms w/ Poisson:  %.4f +/- %.4f' % (np.mean(backlc_0_001.countrate),np.std(backlc_0_001.countrate)))

##### Plot
# start = 0
# end = 10000
# end = end+1

# plt.figure(1)
# plt.clf()
# plt.plot(backlc_time[start:end],backlc_countrate[start:end],'ko',ms=8,label='original 100 ms')
# plt.plot(backlc_0_01_interp_time[start*10:end*10],backlc_0_01_interp_countrate[start*10:end*10],
#          'ro',ms=5,label='interp 10 ms')
# plt.plot(backlc_0_001_interp_time[start*100:end*100],backlc_0_001_interp_countrate[start*100:end*100],
#          'bo',ms=5,label='interp 1 ms')
# plt.plot(backlc_lc.time[start:end],backlc_countrate_wP[start:end],'kd',ms=8,label='original w/ Poisson')
# plt.plot(backlc_lc_0_01.time[start*10:end*10],backlc_countrate_0_01_wP[start*10:end*10],
#          'rd',ms=5,label='10 ms w/ Poisson')
# plt.plot(backlc_lc_0_001.time[start*100:end*100],backlc_countrate_0_001_wP[start*100:end*100],
#          'bd',ms=5,label='1 ms w/ Poisson')
# plt.legend(loc='best')
# plt.show()


#%% ===========================================================================
# Subtract backlc from lc (for entire range of 3.0 - 80.0 keV)
# =============================================================================
lc_0_1 = make_lc(np.array([3.0,80.0]),0.1)
lc_0_01 = make_lc(np.array([3.0,80.0]),0.01)
lc_0_001 = make_lc(np.array([3.0,80.0]),0.001)

def subtract_background(lc,backlc): # Returns a Lightcurve object
    interp_func = interp1d(backlc.time, backlc.counts, kind='cubic', bounds_error=False,
                           fill_value='extrapolate')
    interp_counts = interp_func(lc.time)
    subtracted_counts = lc.counts - interp_counts
    return Lightcurve(time=lc.time, counts=subtracted_counts, mjdref=mjdref)

lc_sub_0_1 = subtract_background(lc_0_1,backlc_0_1)
lc_sub_0_01 = subtract_background(lc_0_01,backlc_0_01)
lc_sub_0_001 = subtract_background(lc_0_001,backlc_0_001)

from stingray.crosscorrelation import AutoCorrelation

ac_0_1 = AutoCorrelation(lc_sub_0_1)
ac_0_01 = AutoCorrelation(lc_sub_0_01)
ac_0_001 = AutoCorrelation(lc_sub_0_001)

ac_0_1.plot()
ac_0_01.plot()
ac_0_001.plot()

