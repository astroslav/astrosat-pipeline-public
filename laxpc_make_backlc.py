#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:46:50 2020

@author: mario
"""

import numpy as np
import os
from subprocess import call
import glob
from context_manager import cd
from symlink_startfiles import symlink
from get_PCUs_dt_layer import get_PCUs_dt_layer
from split_gti import split_gti
from call_laxpc_software import call_spectra, call_backspectra, call_backlc
from recombine_backlc import recombine_backlc

print('')
print('Warning: HEASoft must be initiated prior to running this script!')
print('------------------------------------------------------')
print('This python script will use the LAXPC software to generate back lightcurves.')
print('Making spectra and background spectra is necessary to make backlc.')
print('It will also divide the usergti at higher time resolutions.')

symlink()
usergti = str(input('GTI fits file name? (Default "usergti.fits") > ')) or 'usergti.fits'
level2f = str(input('Level 2 event fits file name? (Default "level2.event.fits") > ')) or 'level2.event.fits'

PCUs, dt, layer = get_PCUs_dt_layer()

print('')
split_gti(float(dt),usergti)

print('')
print('How do you want to generate the BACK lightcurve files?')
print('   1. Minimum GTI interval (e.g. 3.0 to 80.0 keV)')
print('   2. At 1 keV intervals (3-4, 4-5, etc.) WARNING: Takes a long time.')
backlc_choice = str(input('(Default = 1) > ') or 1)
print('')

if os.path.exists('usergti_0000.fits'):
    
    for i in np.arange(0,len(glob.glob('usergti_*.fits')),1):
        print('')
        print('============ gti_'+str(i).zfill(4)+' ============')
        print('')

        if not os.path.exists('gti_'+str(i).zfill(4)):
            os.mkdir('gti_'+str(i).zfill(4))
            print('Made folder for GTI '+str(i).zfill(4))
        else:
            print('Folder "gti_'+str(i).zfill(4)+'" already exists')

        with cd('gti_'+str(i).zfill(4)):
            symlink()
            call('cp ../usergti_'+str(i).zfill(4)+'.fits .', shell=True)
            current_usergti = 'usergti_'+str(i).zfill(4)+'.fits'

    # =============================================================================
    # Make spectra and backspectra
    # =============================================================================
            call_spectra(current_usergti,layer,level2f)
            call_backspectra(current_usergti,layer)

    # =============================================================================
    # Import ebounds from rmf files (to be used for backlc file creation)
    # =============================================================================
            print('')
            print('Importing energy bounds from rmf files.')
            from get_ebounds_from_rmf import ebounds
            np.savetxt('eneinput_1keV_intervals', ebounds, fmt='%.2f')

            ebounds3_80 = np.array([ebounds[0,0],ebounds[-1,-1]])  # The 3-80 keV bounds
            np.savetxt('eneinput_3-80', ebounds3_80, fmt='%.2f')

    # =============================================================================
    # Make back lightcurves
    # =============================================================================
            call_backlc(PCUs,dt,current_usergti,layer,backlc_choice)
    
    # =============================================================================
    # Recombining all backlc's into one file in original working directory   
    # =============================================================================
    recombine_backlc()

else:
    # =============================================================================
    # Make spectra and backspectra
    # =============================================================================
    call_spectra(usergti,layer,level2f)
    call_backspectra(usergti,layer)

    # =============================================================================
    # Import ebounds from rmf files (to be used for backlc file creation)
    # =============================================================================
    print('')
    print('Importing energy bounds from rmf files.')
    from get_ebounds_from_rmf import ebounds
    np.savetxt('eneinput_1keV_intervals', ebounds, fmt='%.2f')
    
    eboundsmin = np.array([ebounds[0,0],ebounds[-1,-1]])  # The 3-80 keV bounds
    np.savetxt('eneinput_3-80', eboundsmin, fmt='%.2f')
    
    # =============================================================================
    # Make back lightcurves
    # =============================================================================
    call_backlc(PCUs,dt,usergti,layer,backlc_choice)
