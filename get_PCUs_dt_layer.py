#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 12:02:26 2020

@author: mario
"""

import sys

def get_PCUs_dt_layer():
    # =============================================================================
    # PCUs used
    # =============================================================================
    print('')
    PCUs = str(input('Which PCU units are you using? (Default = 12 for LAXPC10 and LAXPC20) > ')) or '12'

    if (PCUs == '1' or PCUs == '2' or PCUs == '3'):
        print('Using only LAXPC{}0.'.format(PCUs))
    elif (PCUs == '12' or PCUs == '13' or PCUs == '23'):
        print('Using LAXPC{}0 and LAXPC{}0.'.format(PCUs[0], PCUs[1]))
    elif (PCUs == '123'):
        PCUs = 'all'
        print('Using all instruments.')
    else:
        sys.exit('The input ' + PCUs + ' for the PCUs is not a valid input.')

    # =============================================================================
    # Time resolution
    # =============================================================================
    print('')
    dt = str(input('What time resolution are you using? (Default = 0.1 for dt = 0.1s) > ') or 0.1)
    print('Using a time resolution of {}s'.format(dt))

    # =============================================================================
    # Layer used in each instrument
    # =============================================================================
    print('')
    layer = str(input('How many layers are you using? (Default = all (1 is the other option)) > ')) or '0'

    if layer == '0':
        print('Using all layers.')
    elif layer == '1':
        print('Using the first layer.')
    else:
        sys.exit('The input ' + layer + ' for the layers is not a valid input. Only "0" (or all) or "1" (only first layer) allowed.')

    return PCUs, dt, layer
