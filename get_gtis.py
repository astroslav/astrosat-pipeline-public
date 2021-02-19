#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 11:54:37 2020

@author: mario
"""

import numpy as np
from astropy.io import fits

def get_gtis(gtifits_file):
    '''

    Parameters
    ----------
    gtifits_file : string, optional
        Call the filename of the GTI, preferably the one that is ASBary
        corrected. The default is 'gti_ASBary.fits'.

    Returns
    -------
    gtis : array-like
        All of the GTIs from the laxpc_make_gti.
    gtis_min : array-like
        The minimum of the first GTI and the maximum of the last GTI. Used for
        the purposes of Stingray subtracting the lightcurves and applying all
        the GTIs properly.

    '''
    gtifits = fits.open(gtifits_file)

    try:  # There's only two different names LAXPC uses for GTI as far as I know.
        gtis = np.array(gtifits['GTI'].data.tolist())
    except KeyError:
        gtis = np.array(gtifits['USER_GTI'].data.tolist())

    gtis_min = np.array([[gtis[0,0],gtis[-1,-1]]])
    return gtis, gtis_min
