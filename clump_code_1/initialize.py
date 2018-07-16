#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 15:26:55 2018

@author: sfielder
"""

"""
This piece of code is run to data analyze simulation data

Manual Inputs:
    - Overwrite protection: boolean
        - True or False for if data is overwritten when code runs
    - flist argument: string
        - This argument points to where the data is stored
    - data_check_list: list of strings
        - Decides what simulation time checks to use
    - config_dir: string
        - Points to where the config files are stored
    - data_dir: string
        - Points to where the data should be saved

Returns:
    Tree like directory structure matching different simulation input parameters
    All input parameters are stored via directory naming in the Fiducialxx convention
    
    Directory Tree Output Structure:
        
        tree_top / data_dir / config_x / Fiducialxx / Time Stamp / Output Files /
    
    tree_top - this is where the config files are stored for input
    data_dir - this is set at the highest level of the tree for output
    config_x - this is to seperate different input parameters for the analyzer
    Fiducialxx - this is for different simulation data sets
    Time Stamp - this is for all the checkpoints in the simulation
    Output Files: Generated Currently:
        - FITS File
            - Contains all the necessary clump information for the simulation
            - Projection Plots for axial directions with overlayed center
                of mass coordinates for all the clumps found
        - Header_info.txt
            - Contains information about all the clumps found in each config
            directory found by glob.glob
        - PDF Plots
            - A regular j/j plot
            - A colormapped j/j plot with reference to log_10(mass)
            Notes:
                - Both of these are auto saved within the directory of where
                the fits files are located: /time stamp/
            - j Comparison Plots by Fiducial run per time step
                - These are saved in a created directory in the / Fiducialxx/
                level.

Notes:
    - SEE OVERWRITE PROTECTION BELOW.
"""

import clump_code_1 as cc
import header_printer as hp
import glob
import yaml
import io
import os
from definitinos import jComparisonPlotter
from definitions import FiducialPlotter

"""
Overwrite Protection Here:
    - header and timestep_plots have absolute status and will not even
      start looking for files if value below is set to FALSE.
    
    - analyzer and fiducial_plots have conditional status, and will look for
      for already created directories and will skip those that are found
      if the value below is set to FALSE.
    
    - Setting any of the following values to TRUE will overwrite files
      even if files and directories are found by the code to exist.
"""
overwrite_analyzer = False
overwrite_header = False
overwrite_timestep_plots = True
overwrite_fiducial_plots = True
# =============================================================================
#INPUTS HERE

#Creates a list of directories with the appropriate files for analysis
flist = glob.glob('/mnt/bigdata/erosolow/Orion2/*/data.*.hdf5')

#This is to filter out the timestamps that we want to analyze over
data_check_list = ['0060','0070','0080','0090','0100']

#This is where the config files are
tree_top_dir = '/home/sfielder/Documents/Clumps/'
data_dir = '/home/sfielder/Documents/Clumps/Output/'

# =============================================================================
#Creating empty list for data sorting
flist_data = []
for i in range(0,len(flist)):
    main_string = flist[i].split("/")
    out_string = main_string[-1].split(".")
    time_stamp = out_string[1]
    #This checks if timestamp
    if any(x in time_stamp for x in data_check_list):
        flist_data.append(flist[i])
        print('Directory: ' + flist[i] + ' was detected for analysis.')

#Sorting the Data by filename
flist_data.sort()

#Make a list of all the yaml files found in data_dir
flist_config_yaml = glob.glob(tree_top_dir+'*.yaml')
flist_config_yaml.sort()

carry_on_analyzer = False #Initialization Value - this is the default
for i in range(0,len(flist_config_yaml)):
    #Grabs the config file name here
    naming_string = flist_config_yaml[i].split("/")[-1].split(".")[0]
    #Creating String for Directory
    save_dir = data_dir + naming_string + '/'
    
    #Checking if Directory Exists
    if os.path.isdir(save_dir) == True:
        if overwrite_analyzer == True:
            print("Warning!!! Overwrite has been set to True and Directory: " +
                  save_dir +
                  "is detected as a valid directory. Proceeding with Analysis.")
            carry_on_analyzer = True
            
        else:
           print("Warning!!! Overwrite has been set to False and Directory: " +
                  save_dir +
                  "is detected as a valid directory. Skipping Analysis.")
    else:
        os.mkdir(save_dir)
        carry_on_analyzer = True
        
    #If Carry_on_Analyzer has been set to true, then run the analysis.    
    if carry_on_analyzer == True:
        #Importing Config File settings here
        with io.open(flist_config_yaml[i], 'r') as stream:
            data_loaded = yaml.load(stream)
        
        #Call main code here
        #Testing first file here
        for j in range(0,len(flist_data)):
            cc.analyzer(flist_data[j],
                    data_loaded['l'],
                    data_loaded['cmin'],
                    data_loaded['step'],
                    data_loaded['beta'],
                    data_loaded['clump_sizing'],
                    save_dir)
print('~~~~~~~~~~~~~~~~~~~~~~~~~')            
print('Analysis Section Complete')
print('~~~~~~~~~~~~~~~~~~~~~~~~~')  
# =============================================================================
if overwrite_header == True:

    # Call Header Printer Script to compute the .txt file needed for
    # summary of analysis
    # Grab the config_x directories, and puts them in a list
    flist_config_dir = glob.glob(data_dir + 'config_*')
    
    #Run Loop over all the config_x directories found that are in the list
    
    for i in range(0,len(flist_config_dir)):
        print("Processing Summary File for directory: ", flist_config_dir[i])
        hp.Header_Printer(flist_config_dir[i])

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
    print("Initialize.py file completed.")   
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
# =============================================================================
if overwrite_timestep_plots == True:

    # Comparison Plots by specific timestep happens here:
    flist_plots = glob.glob(data_dir + '**/*.fits', recursive=True)
    flist_plots.sort()
    
    for i in range(0,len(flist_plots)):
        current_file = flist_plots[i]
        jComparisonPlotter(current_file)

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
    print("j Comparison Plots by Timestep (Specific) completed.")    
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
# =============================================================================
# Comparison Plots over Fiducial Runs (by timestep) happens here:
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
            #Calling Main Function Here
            print("Current Timestep being worked on: {}".format(data_check_list_print))
            FiducialPlotter(flist,fid_dir,data_check_list_print)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
print("j Comparison Plots by Timestep (for all Fiducials) completed.")  
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')    
# =============================================================================



        