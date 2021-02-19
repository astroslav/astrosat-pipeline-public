#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 16:02:54 2020

@author: mario
"""

import matplotlib as mpl

# Customizing options here: 
# https://matplotlib.org/3.3.2/tutorials/introductory/customizing.html

def marioplots():
    mpl.rcParams['axes.labelsize'] = 16
    mpl.rcParams['errorbar.capsize'] = 4
    mpl.rcParams['font.family'] = 'Times New Roman'
    
    mpl.rcParams['lines.markersize'] = 8
    mpl.rcParams['lines.markeredgecolor'] = 'k'
    mpl.rcParams['lines.markeredgewidth'] = 0.5
    
    mpl.rcParams['xtick.minor.visible'] = True
    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['xtick.labelsize'] = 'large'
    
    mpl.rcParams['ytick.minor.visible'] = True
    mpl.rcParams['ytick.right'] = True
    mpl.rcParams['ytick.direction'] = 'in'
    mpl.rcParams['ytick.labelsize'] = 'large'
    
    return print('Set matplotlib rcParams to my faves.')