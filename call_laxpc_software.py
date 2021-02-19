#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 13:45:17 2020

@author: mario
"""

import os
import sys
from subprocess import call

def call_spectra(usergti,layer,level2event):
    if not os.path.exists('spectrum_10.pha' or 'spectrum_10_L1.pha'):
        print('')
        print('============ SPECTA ============')
        call(['laxpc_make_spectra -u '+usergti+' -l '+layer+' '+level2event],shell=True)
        print('Spectra made.')
    else:
        print('Spectra already exists. Delete \'spectra_10*.pha\' if you want to overwrite them.')


def call_backspectra(usergti,layer):
    if not os.path.exists('backlxp10.pha'):
        print('')
        print('============ BACK SPECTRA ============')
        call(['laxpc_make_backspectra -u '+usergti+' -l '+layer+' filterfiles'],shell=True)
        print('Back spectra made.')
    else:
        print('Back spectra already exists. Delete \'backlxp10.pha\' if you want to overwrite them.')
        

def call_backlc(PCUs,dt,usergti,layer,backlc_choice):
    if not os.path.exists('Back_lightcurve.txt'):
        print('')
        print('============ BACK LC ============')
        if backlc_choice == '1':
            call(['laxpc_make_backlightcurve -p '+PCUs+' -t '+dt+' -u '+usergti
                  + ' -e eneinput_3-80 -l '+layer+' filterfiles'],shell=True)
        elif backlc_choice == '2':
            call(['laxpc_make_backlightcurve -p '+PCUs+' -t '+dt+' -u '+usergti
                  + ' -e eneinput_1keV_intervals -l '+layer+' filterfiles'],shell=True)
        else:
            sys.exit('Invalid backlc_choice = '+backlc_choice+'.')
    
        if os.path.exists('Back_lightcurve.txt'):
            print('Back lightcurves made.')
        else:
            sys.exit('Back lightcurves not made. LAXPC didn\'t like something.')
    else:
        print('Back lightcurves already exists. Check files.')
