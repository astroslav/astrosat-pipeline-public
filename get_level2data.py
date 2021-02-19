#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:09:21 2020

@author: mario
"""

import numpy as np
from astropy.io import fits

def get_level2data(evtfits_file):
    '''

    Parameters
    ----------
    evtfits_file : str, optional
        Call the filename of the event file, preferably the one that is ASBary
        corrected. The default is 'level2_ASBary.fits'.

    Returns
    -------
    level2data : array-like
        Returns five columns of data. They are time, channel, layer, laxpc_no,
        and energy, respectively.
    ncounts : int
        The length of each of the five arrays in level2data.
    mjdref : int
        A reference to the Modified Julian Date. Needed when making EventList.

    '''
    level2fits = fits.open(evtfits_file)
    ncounts = level2fits['event file'].header['NAXIS2']
    mjdref = level2fits['event file'].header['MJDREFI']
    level2data = np.array(level2fits['event file'].data.tolist())
    return level2data, ncounts, mjdref
