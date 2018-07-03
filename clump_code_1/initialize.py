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

Notes:
    - Both functions being called below are seperated out thus;
    If overwrite for data_analysis is turned off, then script prints a statement
    out and then moves forward onto the header_printer function
"""


import clump_code_1 as cc
import FITS_Header_Printer as hp
import glob
import yaml
import io
import os

# Overwrite Protection Here
# Set to True - Will Overwrite Data already on disk
# Set to False - Will NOT Overwrite Data already on disk
# This is specific to the data analysis portion
overwrite = True



#Creates a list of directories with the appropriate files for analysis
flist = glob.glob('/mnt/bigdata/erosolow/Orion2/*/data.*.hdf5')
#Creating empty list for data sorting
flist_data = []
#This is to filter out the timestamps that we want to analyze over
data_check_list = ['0060','0070','0080','0090','0100']

for i in range(0,len(flist)):
    main_string = flist[i].split("/")
    out_string = main_string[-1].split(".")
    time_stamp = out_string[1]
    #This checks if timestamp
    if any(x in time_stamp for x in data_check_list):
        flist_data.append(flist[i])
        print('Directory: ' + flist[i] + ' was grabbed for analysis.')

#Sorting the Data by filename
flist_data.sort()

#This is where the config files are
tree_top_dir = '/home/sfielder/Documents/Clumps/'
data_dir = '/home/sfielder/Documents/Clumps/Output/'

#Make a list of all the config files found in  data_dir
flist_config_files = glob.glob(tree_top_dir+'*.yaml')
flist_config_files.sort()

for i in range(0,len(flist_config_files)):
    #Grabs the config file name here
    naming_string = flist_config_files[i].split("/")[-1].split(".")[0]
    #Creating String for Directory
    save_dir = data_dir + naming_string + '/'
    
    #Checking if Directory Exists, if so overwrite is checked and either analysis
    #is skipped or overwritten
    if os.path.isdir(save_dir) == True:
        if overwrite == True:
            print("Warning!!! Overwrite has been set to True and Directory: " +
                  save_dir +
                  "is detected as a valid directory. Proceeding with Analysis.")
            #Importing Config File settings here
            with io.open(flist_config_files[i], 'r') as stream:
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
        else:
           print("Warning!!! Overwrite has been set to False and Directory: " +
                  save_dir +
                  "is detected as a valid directory. Skipping Analysis.")
    else:
        #If Directory is not present, then it gets created here and the code
        #is ran through the analyzer
        os.mkdir(save_dir)

        #Importing Config File settings here
        with io.open(flist_config_files[i], 'r') as stream:
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Call Header Printer Script to compute the .txt file needed for summary of analysis
#Grab the config_x directories, and puts them in a list
flist_config_dir = glob.glob(data_dir + 'config_*')
print("Debug: flist_config_dir set as: ", flist_config_dir)

#Run Loop over all the config_x directories found that are in the list
for i in range(0,len(flist_config_dir)):
    print("Processing Summary File for directory: ", flist_config_dir[i])
    hp.Header_Printer(flist_config_dir[i])

print("Initialize.py file completed.")    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

        