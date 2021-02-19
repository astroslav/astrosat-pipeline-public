#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:37:00 2020

@author: mario
"""

import sys
import glob
import numpy as np
from astropy.io import fits
from find_nearest_idx import find_nearest_idx

# =============================================================================
# Find and read rmf files
# =============================================================================

def get_ebounds_from_rmf(rmf_string):
    rmf_file = fits.open(rmf_string)
    rmf_data = np.array(rmf_file['EBOUNDS'].data.tolist())
    emin = rmf_data[:,1]
    emax = rmf_data[:,2]
    return np.array([emin,emax]).T

ebounds10 = get_ebounds_from_rmf(glob.glob('lx10cshp*')[0])
ebounds20 = get_ebounds_from_rmf(glob.glob('lx20cshp*')[0])
ebounds30 = get_ebounds_from_rmf(glob.glob('lx30cshp*')[0])

# =============================================================================
# Start with getting all the indices for ebounds20 using one of two methods:
#    EXCLUDING extremes: If any values are contained within a bin with a range
#                        outside the min/max, that bin is ignored
#    INCLUDING extremes: If any values are contained within a bin with a range
#                        outside the min/max, that bin is included
# =============================================================================

extreme_method = str(input('Do you want to include the extremes of the energy range? Default: n (y/n) > ')) or 'n'
extreme_method = extreme_method.lower()
if extreme_method == 'n':
    print('EXCLUDING extremes of energy ranges')
elif extreme_method == 'y':
    print('INCLUDING extremes of energy ranges')
else:
    sys.exit('Invalid method choice. Use "y" or "n".')


def get_ebounds20idx():
    # Explicitly defined as 3-80 keV. This is the working range for LAXPC
    vals3to80 = np.arange(3,81,1) 
    idx = []
    for val in vals3to80:
        nearest20idx = find_nearest_idx(ebounds20[:,0],val)

        # =============================================================================
        # Dealing with extremes (e.g. including 3.0 for a 3.0-80.0 range)
        # =============================================================================
        if val == vals3to80[0]:
            if extreme_method == 'n' and ebounds20[:,0][nearest20idx] < val:
                nearest20idx += 1
            elif extreme_method == 'y' and ebounds20[:,0][nearest20idx] > val:
                nearest20idx -= 1

        if val == vals3to80[-1]:
            if extreme_method == 'n' and ebounds20[:,0][nearest20idx] > val:
                nearest20idx -= 1
            elif extreme_method == 'y' and ebounds20[:,0][nearest20idx] < val:
                nearest20idx += 1

        nearest20 = ebounds20[:,0][nearest20idx]
        idx.append(np.where(ebounds20[:,0] == nearest20)[0][0])

    return np.asarray(idx)

ebounds20lims = ebounds20[get_ebounds20idx(),0]

# =============================================================================
# Round up so that the values chosen can be included in the range we'll define
# =============================================================================
ebounds20ceil = np.ceil(ebounds20lims*100)/100

# =============================================================================
# Now, get the indices for the ebounds10 and ebounds30 where we have the the
# last values that are smaller than the specified value of ebounds20. This
# will be used to make sure that all of them are included the defined range
# =============================================================================

def get_eboundsidx(ebounds,ebounds20ceil):
    idx = []
    for val in ebounds20ceil:
        idx.append(np.where(ebounds[:,0] < val)[0][-1])
    return np.asarray(idx)

ebounds10lims = ebounds10[get_eboundsidx(ebounds10,ebounds20ceil),0]
ebounds30lims = ebounds30[get_eboundsidx(ebounds30,ebounds20ceil),0]

# =============================================================================
# Vertically stack the arrays so we can now find the lower bound of what needs
# to be included to make sure the energy bounds are included properly. We'll
# find the minimum of all three and then round down to make sure we include
# everything. I'll have to check to make sure there are no overlaps afterward.
# =============================================================================

def get_ebounds(ebounds10lims,ebounds20lims,ebounds30lims,decimals=2):
    eboundslims = np.vstack((ebounds10lims,ebounds20lims,ebounds30lims))

    # Get the specified decimal places by using the operations below
    eboundslimsMIN = np.floor(np.min(eboundslims,axis=0)*(10**decimals))/(10**decimals)
    eboundslimsMAX = np.ceil(np.max(eboundslims,axis=0)*(10**decimals))/(10**decimals)

    # Mental check. eboundslimsMAX should be equal to ebounds20ceil
    if (ebounds20ceil == eboundslimsMAX).all() == True:
        print('Everything is gucci. (ebounds seem correctly set-up)')
    else:
        print('Something\'s not gucci. (ebounds are incorrectly written)')

    return np.array([eboundslimsMIN[0:-1],eboundslimsMAX[1:]]).T

ebounds = get_ebounds(ebounds10lims,ebounds20lims,ebounds30lims)

# =============================================================================
# Check the counts of each row in the given ebounds ranges. All should be used
# only once so that we don't double count counts. (tongue twister?)
# =============================================================================

def check_arr_counts(ebounds_rmf_ranges,ebounds):
    ebounds_inside_lim = []
    for row in ebounds_rmf_ranges:
        if min(ebounds[:,0]) < row[0] and row[1] < max(ebounds[:,1]):
            ebounds_inside_lim.append(row)

    ebounds_inside_lim = np.asarray(ebounds_inside_lim)
    ebounds_arr_counts = np.zeros(len(ebounds_inside_lim))

    for row in ebounds:
        for i in range(len(ebounds_inside_lim)):
            if row[0] < ebounds_inside_lim[i,0] and ebounds_inside_lim[i,1] < row[1]:
                ebounds_arr_counts[i] += 1

    if (ebounds_arr_counts == 1).all() == True:
        return True  # all counts are used only once
    else:
        return False  # some counts are counted more than once or not at all

if check_arr_counts(ebounds10,ebounds) *\
   check_arr_counts(ebounds20,ebounds) *\
   check_arr_counts(ebounds30,ebounds) == True:
    print('All counts counted once.')
else:
    sys.exit('I counted the counts wrong :(')
