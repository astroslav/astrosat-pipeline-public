#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 13:40:51 2020

@author: mario
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from context_manager import cd
import glob

def recombine_backlc(bkg_lc_fname='Back_lightcurve.txt'):
    # Check if combined_backlc.txt already exists
    if os.path.exists('combined_backlc.txt'):
        return print('"combined_backlc.txt" found. Skipping recombination.')

    time = []
    count_rate = []
    number_of_rows = len(glob.glob('usergti_*.fits'))
    for i in np.arange(0,number_of_rows,1):
        if not os.path.exists('gti_'+str(i).zfill(4)):   
            sys.exit('No gti_'+str(i).zfill(4)+' folder exists in this dir.')
        with cd('gti_'+str(i).zfill(4)):
            if i == 0:
                print('Currently reading and appending gti '+str(i)+' out of '+str(number_of_rows-1))
            else:
                print('... gti '+str(i)+' out of '+str(number_of_rows-1))

            data = np.loadtxt(bkg_lc_fname, comments='!')
            time.append(data[:,3])
            count_rate.append(data[:,1])
    
    time = np.concatenate(np.asarray(time))
    count_rate = np.concatenate(np.asarray(count_rate))
    np.savetxt('combined_backlc.txt',np.array([time,count_rate]).T,fmt='%.8f %.14f')
    print('"combined_backlc.txt" saved in this directory.')
    
    plot_userinput = str(input('Would you like to see the plot? (Default: n) > ')) or 'n'
    
    if plot_userinput.lower() == 'y':
        plt.figure(1)
        plt.clf()
        plt.plot(time,count_rate,'k.')
        plt.show()
    else:
        print('Skipping plot.')
