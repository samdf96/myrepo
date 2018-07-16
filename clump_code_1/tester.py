#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:12:57 2018

@author: sfielder
"""
import glob
import os
from definitions import jComparisonPlotter
from definitions import FiducialPlotter

data_dir = '/Users/sfielder/Documents/Astro_Data/'

# Comparison Plots by specific timestep happens here:
flist_plots = glob.glob(data_dir + '**/*.fits', recursive=True)
flist_plots.sort()

for i in range(0,len(flist_plots)):
    current_file = flist_plots[i]
    jComparisonPlotter(current_file)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
print("j Comparison Plots by Timestep (Specific) completed.")   
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# =======================================================

# Comparison Plots over Fiducial Runs (by timestep) happens here:
data_check_list = ['0060','0070','0080','0090','0100']
overwrite_fiducial_plots = True
carry_on_fiducial_plots = False #Initialize value - this is the default
flist_config = glob.glob(data_dir+'config_*')
flist_config.sort()

for k in range(0,len(flist_config)): #Adding Loop for config directories
    #Write in os function to create appropiate directory for Fiducial Plots
    fid_dir = flist_config[k]+'/Fiducial_Plots/'
    if os.path.isdir(fid_dir) == True:
        print('Fiducial Directory Detected.\n')
        if overwrite_fiducial_plots==False:
            print('Overwrite set to False. Plots will not be made.\n')
        else:
            carry_on_fiducial_plots = True
            print('Overwrite set to True. Continuing with plot making.\n')
    else:
        os.mkdir(fid_dir)
        carry_on_fiducial_plots = True

    #If Carry On Value is True then continue with plotting
    if carry_on_fiducial_plots == True:
        for i in range(0,len(data_check_list)): #len(data_check_list)
            flist = glob.glob(flist_config[k]+'/**/*'+data_check_list[i]+'*.fits',
                          recursive=True)
            flist.sort()
            data_check_list_print = data_check_list[i]
            print("Current Timestep being worked on: {}".format(data_check_list_print))
            #Calling Main Function Here
            FiducialPlotter(flist,fid_dir,data_check_list_print)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
print("j Comparison Plots by Timestep (for all Fiducials) completed.")   
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  