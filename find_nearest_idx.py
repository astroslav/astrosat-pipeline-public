#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 15:41:59 2020

@author: mario
"""

import numpy as np

def find_nearest_idx(array, value):
    '''

    Parameters
    ----------
    array : array-like
        The array to be checked
    value : float
        The value in the array you're searching for

    Returns
    -------
    idx : int
        The index in the array closest to the given value

    '''
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
