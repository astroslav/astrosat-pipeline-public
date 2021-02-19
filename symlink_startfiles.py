#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 15:05:54 2020

@author: mario
"""

import os
from subprocess import call

def symlink():
    print('')
    if not os.path.exists('level2.event.fits'):
        print('Symbolically linking files to this directory.')
        call(['ln -s /home/mario/Documents/repos/astrosat-pipeline/start_files/* .'], shell=True)
    else:
        print('There are linked files in this directory already. Delete \'level2.event.fits\' if you want to overwrite.')
