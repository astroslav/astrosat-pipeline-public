#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 12:53:22 2020

@author: mario
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from get_gtis import get_gtis
from stingray.lightcurve import Lightcurve

# This assumes that the usergti.fits file is in the working directory
gtis, gtis_min = get_gtis()

fname = str(input('What is the name of the background lc file? (Default = Back_lightcurve_3.0_80.0keV.lc) ')) or "Back_lightcurve_3.0_80.0keV.lc"

backLC = fits.open(fname)
backLCdata = np.array(backLC['lightcurve'].data.tolist())
backLCtime      = backLCdata[:,0]
backLCcountrate = backLCdata[:,1]
backLCerr       = backLCdata[:,2]

back_lc = Lightcurve(backLCtime, backLCcountrate, err=backLCerr, input_counts=False, gti=gtis_min)

plt.figure(1, clear=True)
plt.errorbar(back_lc.time, back_lc.countrate, yerr=back_lc.countrate_err,
             ls='none', marker='.', color='r', markersize=2,
             elinewidth=0.2)
plt.xlabel('Times')
plt.ylabel('Counts')
plt.show()
