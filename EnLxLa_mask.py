#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:46:57 2020

@author: mario
"""

import numpy as np

def EnLxLa_mask(energy_arr, laxpc_no_arr, layer_arr,
                energy_range=np.array([3.0,80.0]),
                laxpc_no_range=np.array([1,2]),
                layer_range=np.array([1,2,3,4,5])):
    '''

    Parameters
    ----------
    energy_arr : array-like
        Input laxpc array for energy
    laxpc_no_arr : array-like
        Input laxpc array for the laxpc instrument number. 1 = 10, 2 = 20, 3 = 30.
    layer_arr : array-like
        Input laxpc array for the laxpc layer number. 1 = top layer, 5 = bottom layer.
    energy_range : array-like, optional
        Which energies to include in the mask. The default is 3.0 keV to 80 keV.
    laxpc_no_range : array-like, optional
        Which laxpc instruments to include in the mask. The default is laxpc10 and laxpc20.
    layer_range : array-like, optional
        Which laxpc layers to include in the mask. The default is all five.

    Returns
    -------
    array-like, boolean
        A mask that only includes only where all elements are True.

    '''
    # Energy mask
    en_mask_low = energy_range[0] <= energy_arr
    en_mask_high = energy_arr <= energy_range[1]
    en_mask = en_mask_low * en_mask_high

    # Laxpc instrument number
    lx_mask = np.isin(laxpc_no_arr, laxpc_no_range)

    # Laxpc layer number
    la_mask = np.isin(layer_arr, layer_range)

    return en_mask*lx_mask*la_mask
