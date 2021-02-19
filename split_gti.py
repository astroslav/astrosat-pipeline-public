#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 10:42:31 2020

@author: mario
"""

import numpy as np
import sys, os
from astropy.table import Table

def split_gti(dt,usergti_file):
    # Quick check for already existing files
    if os.path.exists('usergti_0000.fits'):
        return print('Custom "usergti_*.fits" files exist. Delete "usergti_0000.fits" to overwrite them.')

    # Index "j" is used for iterating through writing fits files
    j = 0
    table = Table.read(usergti_file)

    # The laxpc_make_backlightcurve software has a maximum limit of 1000000
    # points that can be ran through the software at once. If we want higher
    # time resolution, we have to make sure all our GTI's with our chosen
    # dt are less than that. I've arbitrarily chosen 90% of the max to be the
    # limit to avoid any boundary problems. -shrug-
    laxpc_max = 1000000 * 0.9

    # Total time is calculated as the (max of GTI - min of GTI)
    total_time = table[-1][-1] - table[0][0]

    if total_time/dt/laxpc_max > 1:
        print('LAXPC would raise a "too many points" error with this dt.')
        print('Splitting usergti into individual rows.')

        for row in table:
            row_diff = row[1] - row[0]  # row_diff = STOP - START

            factor = row_diff/dt/laxpc_max

            if factor > 1:
                # The following two lines are used for making intermediary rows.
                row0 = row[0]
                step = row_diff/np.ceil(factor)

                for i in range(int(np.ceil(factor))):
                    # The following for loop modifies the table to a small enough
                    # GTI so we can bypass that 1000000 point limit at the given dt.
                    row[0] = row0 + step*(i)
                    row[1] = row0 + step*(i+1)

                    # Writing all the rows using the index "j" starting from 0.
                    table[[row.index]].write('usergti_'+str(j).zfill(4)+'.fits',overwrite=True)
                    j += 1
            else:
                table[[row.index]].write('usergti_'+str(j).zfill(4)+'.fits',overwrite=True)
                j += 1

        print(usergti_file+' split into manageable chunks for LAXPC.')

    else:
        print('No need to split up GTI into smaller chunks.')
