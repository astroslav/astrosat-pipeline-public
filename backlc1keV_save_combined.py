#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:11:07 2020

@author: mario
"""

import numpy as np
import glob
from astropy.io import fits
from find_nearest_idx import find_nearest_idx



def get_combined_countrate(eneinput):
    '''

    Parameters
    ----------
    eneinput : str
        A text file with two columns with minimum and maximum energies in the
        first and second column respectively

    Returns
    -------
    time : array-like
        Array of times from background lightcurves
    countrate : array-like
        Array of count rates from background lightcurves

    '''
    filelist = []
    for emin, emax in eneinput:
        filelist.append(glob.glob('Back_lightcurve_%.1f_%.1fkeV.lc' % (emin,emax)))

    back_lc_countrate_matrix = []
    for file in filelist:
        back_lc_file = fits.open(file[0])  # the [0] is because glob.glob adds another set of [] brackets
        back_lc_data = np.array(back_lc_file['lightcurve'].data.tolist())
        back_lc_countrate_matrix.append(back_lc_data[:,1])  # the 1st column is countrate

    # get the times
    back_lc_times = back_lc_data[:,0]

    # make the matrix a proper numpy array
    back_lc_countrate_matrix = np.asarray(back_lc_countrate_matrix)

    # sum over all columns (i.e. sum all countrates at the same time)
    back_lc_countrate_sum = np.sum(back_lc_countrate_matrix,axis=0)

    return back_lc_times, back_lc_countrate_sum



def combined_backlc(emin,emax,eneinput_file):
    '''

    Parameters
    ----------
    emin : float
        Minimum energy used (in keV)
    emax : float
        Maxumum energy used (in keV)
    eneinput_file : str
        A text file with two columns with minimum and maximum energies in the
        first and second column respectively

    Returns
    -------
    time : array-like
        Array of times from background lightcurve
    countrate : array-like
        Array of count rates from background lightcurve
    eneinput : array-like
        New energy range array that is specific to the emin and emax used here

    '''
    eneinput_python = np.loadtxt(eneinput_file)

    emin_idx = find_nearest_idx(eneinput_python[:,0],emin)
    emax_idx = find_nearest_idx(eneinput_python[:,1],emax)

    print('The actual range read is %.2f to %.2f keV.'
          % (eneinput_python[emin_idx,0],eneinput_python[emax_idx,1]))

    eneinput = eneinput_python[emin_idx:emax_idx+1]
    time, countrate = get_combined_countrate(eneinput)

    return time, countrate, eneinput



print('')
print('This python script will save a "combined_backlc_emin_emax.txt" in the current folder')
print('')

emin = float(input('Input the lower energy bound for the lightcurve: '))
emax = float(input('Input the upper energy bound for the lightcurve: '))
eneinput = str(input('What is the name of the eneinput file? (Default = eneinput_1keV_intervals) ')) or "eneinput_1keV_intervals"

back_time, back_countrate, eneinput = combined_backlc(emin,emax,eneinput)
backlc_textsave = np.array([back_time,back_countrate,np.sqrt(back_countrate)]).T

np.savetxt('combined_backlc_%.1f_%.1f.txt' % (emin,emax),backlc_textsave,fmt='%.8f')
print('"combined_backlc_%.1f_%.1f.txt" saved.' % (emin,emax))
print('')
