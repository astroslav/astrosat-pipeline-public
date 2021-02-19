#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 13:22:31 2020

@author: mario
"""

import numpy as np
import matplotlib.pyplot as plt

fname = str(input('What is the name of the backlc txt file? (Default = Back_lightcurve.txt) ')) or "Back_lightcurve.txt"

data = np.loadtxt(fname, comments='!')

time = data[:,3]
count_rate = data[:,1]
# count_rate_err = data[:,2]

plt.figure(1, clear=True)
plt.plot(time,count_rate,'k.')
plt.show()
