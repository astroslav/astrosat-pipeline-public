#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 13:32:51 2021

@author: mario
"""

import os
import sys
import glob
from context_manager import cd
import numpy as np
from xspec import *
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
import time
import corner
from astropy.table import Table
from astropy.table import vstack
# import bxa.xspec as bxa
print('imported all modules')

if not os.path.exists('usergti_0000.fits'):
    sys.exit('usergti_0000.fits file doesn\'t exist in this dir')


def define_model_params(temp, model):
    
    model.TBabs.nH = 0.6
    model.TBabs.nH.frozen = True
    model.nthComp.Gamma.values = "1.626, 0.01, 0.8, 1.2, 2.0, 2.2"
    model.nthComp.kT_e = 400.0
    model.nthComp.kT_e.frozen = True
    model.nthComp.kT_bb.link = "11"
    model.nthComp.inp_type = 1.0
    model.nthComp.norm = 0.02049
    model.gaussian.LineE = "7.189, 0.1, 5.0, 5.5, 8.5, 9.0"
    model.gaussian.Sigma = 0.01
    model.gaussian.Sigma.frozen = True
    model.gaussian.norm = 0.000154
    model.diskbb.Tin = str(temp) + ", 0.1, 0.05, 0.2, 5, 6"
    # model.diskbb.Tin = "0.4266, 0.1, 0.1, 0.2, 1.8, 2.0"
    model.diskbb.norm = "179.6, 0.1, 0, 0, 1e20, 1e24"
    model.gaussian_5.LineE = "51.56, 1, 30, 35, 70, 75"
    model.gaussian_5.Sigma = "4.400, 0.05, 1, 2, 20, 30"
    model.gaussian_5.norm = 0.000677
    
    model.nthComp.Gamma.frozen = True
    model.gaussian.LineE.frozen = True
    model.gaussian.norm.frozen = True 
    model.diskbb.Tin.frozen = True  
    model.gaussian_5.LineE.frozen = True 
    model.gaussian_5.Sigma.frozen = True
    model.gaussian_5.norm.frozen = True
    
    Fit.perform()
    
    model.nthComp.Gamma.frozen = False
    model.gaussian.LineE.frozen = False
    model.gaussian.norm.frozen = False 
    model.diskbb.Tin.frozen = False  
    model.gaussian_5.LineE.frozen = False 
    model.gaussian_5.Sigma.frozen = False
    model.gaussian_5.norm.frozen = False
    
    


firstPass = True # var used for creating DataFrames

# Delete xspec_data.csv and xspec_time.csv if you want to remake the data
# with different parameters
if not os.path.exists('xspec_data0000.csv'):
    print('PyXspec data not found. Time to make it.')
    
    # Delete previous chain*.fits files to avoid having to confirm rewriting
    # files from the ChainManager using MCMC
    filelist = glob.glob('chain*.fits')
    for filepath in filelist:
        try:
            os.remove(filepath)
        except:
            sys.exit("Error while deleting file: ", filepath)
    
    # Defining Fit statistics. cstat used for MCMC (doesn't like chisq?)
    Fit.nIterations = 100
    Fit.query = 'yes'
    Fit.statMethod = "cstat"
    
    # Defining MCMC statistics
    AllChains.defAlgorithm = 'gw'
    AllChains.defBurn = 100000
    AllChains.defLength = 100000
    AllChains.defWalkers = 50
    AllChains.defProposal = 'gaussian fit'
    AllChains.defRand = True
    
    
    # Defining Plotting parameters
    Plot.device = "/xw"
    Plot.xAxis = "keV"
    # Plot.yLog = True
    Plot.add = True

    # n is the number of usergti_0*.fits files in the current working dir. 
    n = len(glob.glob('usergti_0*.fits')) 
    
    for i in range(n):
    # for i in [0]:
        # Clearing all the data, models, and chains from previous loops
        AllData.clear()
        AllChains.clear()
        
        # I couldn't figure out a better way to call in the data
        usergti_fits = fits.open('usergti_'+str(i).zfill(4)+'.fits')
        AllData += 'gti_'+str(i).zfill(4)+'/spectrum_grp_20.pha'.format(i)
        
        # if i < 10:
        #     usergti_fits = fits.open('usergti_000{}.fits'.format(i))
        #     AllData += "gti_000{}/spectrum_grp_20.pha".format(i)
        # if i >= 10:
        #     usergti_fits = fits.open('usergti_00{}.fits'.format(i))
        #     AllData += "gti_00{}/spectrum_grp_20.pha".format(i)
        
        # We're ignoring all data in the spectrum less than 3.0 keV and higher 
        # than 80.0 keV. This also does not include the extremes (3.0 and 80.0)
        # are not part of the data. If you want to include them, use the 
        # function "notice"
        AllData.ignore("**-3.0, 80.0-**")
        
        # Create models and chains for different values of Tin
        c_list = []

        temps = np.arange(0.2,1.7,0.2)
        # temps = [3.0]
        for j in range(len(temps)):
            AllModels.clear()
            mod = Model("tbabs*(nthComp+gauss+diskbb+gauss)")
            define_model_params(temps[j], mod)
            # Fit.renorm()
            Fit.perform()
            # input('wait...')
            c_list.append(Chain("chain_temp{}.fits".format(j)))
        
        # Write all chains to a single FITS file
        base_table = Table.read("chain_temp0.fits", format='fits')
        for k in range(1,len(temps)):
            table = Table.read("chain_temp{}.fits".format(k), format='fits')
            base_table = vstack([base_table,table])
        base_table.write('chain{}.fits'.format(i), format='fits')
        
        # Remove all temporary chain fits files
        [os.remove(filepath) for filepath in glob.glob('chain_temp*.fits')]
        
        # Record all the time data
        t = {'Time_start': usergti_fits['USER_GTI'].data['START'][0],
             'Time_end': usergti_fits['USER_GTI'].data['STOP'][0]
             }
        
        # Open the chain{}.fits data to write to Pandas DataFrame
        chain_data = fits.open("chain{}.fits".format(i))
        d = {'nthComp_Gamma': chain_data['CHAIN'].data['Gamma__2'],
             'nthComp_norm': chain_data['CHAIN'].data['norm__7'],
             'gaussian_LineE': chain_data['CHAIN'].data['LineE__8'],
             'gaussian_norm': chain_data['CHAIN'].data['norm__10'],
             'diskbb_Tin': chain_data['CHAIN'].data['Tin__11'], 
             'diskbb_norm': chain_data['CHAIN'].data['norm__12'],
             'gaussian2_LineE': chain_data['CHAIN'].data['LineE__13'],
             # 'gaussian2_Sigma': chain_data['CHAIN'].data['Sigma__14'],
             'gaussian2_norm': chain_data['CHAIN'].data['norm__15'],
             'FIT_STATISTIC': chain_data['CHAIN'].data['FIT_STATISTIC']
             }    
        
        # The first index creates the DataFrame table, the rest of the indices
        # appends to the table
        if firstPass:  
            t_df = pd.DataFrame(t, index=[0])
            # d_df = pd.DataFrame(d)
            firstPass = False
        else:
            t_df = pd.concat([t_df, pd.DataFrame(t, index=[i])], ignore_index=True)
            # d_df = pd.concat([d_df, pd.DataFrame(d)], ignore_index=True)
    
        # Within each loop, write all the chain data to *.csv files
        d_df = pd.DataFrame(d)
        d_df.to_csv('xspec_data'+str(i).zfill(4)+'.csv')
        
    t_df.to_csv('xspec_time.csv')
    print('xspec_time.csv and xspec_data*.csv created.')
    
else:
    n = len(glob.glob('xspec_data0*.csv')) 
    d_df_list = []
    for i in range(n):
        d_df_list.append(pd.read_csv('xspec_data'+str(i).zfill(4)+'.csv'))
        
    t_df = pd.read_csv('xspec_time.csv')
    print('xspec_time.csv and xspec_data.csv loaded.')


    
print('')
print('The data columns are as follows')
print(d_df_list[0].columns)

# Defining the one, two, and three sigma percentiles
onesigma = [0.1587, 0.5, 0.8414]
twosigma = [0.0228, 0.5, 0.9773]
threesigma = [0.0014, 0.5, 0.9987]

len_chain = 100000

# Function to call a specific parameter name and plot it
def plot_stat(varname, time_df=t_df, sigma=onesigma, n_fits=1, lenChain=1000, color='k'):
    varname_stats = []
    for i in range(n_fits):
        varname_stats.append(d_df[varname][i*lenChain:(i+1)*lenChain].describe(percentiles=sigma))

        plt.vlines(time_df.Time_start[i], varname_stats[i][6], varname_stats[i][4], colors=color)
        plt.plot(time_df.Time_start[i], varname_stats[i][1], color+'o')


temps = np.arange(0.2,1.7,0.2)
# temps = [3.0]
n_Tins = len(temps)
n_length = len_chain




# plt.figure(1, clear=True)
# 
# # # plot_stat('nthComp_Gamma',sigma=threesigma, color='g')
# # # plot_stat('nthComp_Gamma',sigma=twosigma, color='r')
# # # plot_stat('nthComp_Gamma',sigma=onesigma, color='k')
# # 
# plot_stat('diskbb_Tin',sigma=threesigma, color='g')
# plot_stat('diskbb_Tin',sigma=twosigma, color='r')
# plot_stat('diskbb_Tin',sigma=onesigma, color='k')
# 
# plt.xlabel('Time (s)', size=14)
# plt.show()


for d_df in d_df_list:
    data_Tin = d_df.gaussian2_norm.values.reshape(n_Tins,n_length)
    data_FIT = d_df.FIT_STATISTIC.values.reshape(n_Tins,n_length)

    fig, axs = plt.subplots(n_Tins, 1, sharex=True, sharey=True, figsize=(10,8))
    fig.subplots_adjust(left=0.05, bottom=0.1, right=0.97, top=0.97, hspace=0)
    minFIT = data_FIT.min()
    print('minFIT = ', minFIT)

    print(len(np.where(d_df.FIT_STATISTIC.values < minFIT+1)[0]))

    for i in range(n_Tins):
    #     axs[i].plot(np.arange(0,n_length,1),data_Tin[i,:],'k.',ms=0.5)
    #     axs[i].set_ylabel('Tin = %.1f' % (temps[i]))
    #     axs[i].set_ylim(40,70)
    # 
    # plt.xlim(0,n_length)
    # plt.xlabel('Steps',size=14)
    # plt.tight_layout()
    # plt.show()

        axs[i].plot(data_Tin[i,:],data_FIT[i,:],'k.',ms=0.5)
        # axs.plot(np.arange(0,n_length,1),data_FIT[i,:],'k.',ms=0.5)
        # axs[i].plot(np.arange(0,n_length,1),data_Tin[i,:],'k.',ms=0.5)

        # axs[i].set_ylabel('%.1f' % (temps[i]), size=10)



    plt.ylim(minFIT,minFIT+2.7055)
    # plt.xlim(0,2)
    # plt.xlabel('Steps',size=14)
    # plt.tight_layout()
    plt.show()
    
    input('wait...')




# data_Tin = d_df.diskbb_Tin.values.reshape(16,10000)
# 
# figure = corner.corner(data_Tin.T, quantiles=[0.16, 0.5, 0.84])
# plt.show()
