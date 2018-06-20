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
        
        / data_dir / Config_x / Fiducialxx / Time Stamp / Output Files /
        
    data_dir - this is set at the highest level of the tree
    config_x - this is to seperate different input parameters for the analyzer
    Fiducialxx - this is for different simulation data sets
    Time Stamp - this is for all the checkpoints in the simulation
    Output Files: Generated Currently:
        - FITS File
            - Contains all the necessary clump information for the simulation
            - Projection Plots for axial directions with overlayed center
                of mass coordinates for all the clumps found
"""


import clump_code_1 as cc
import glob
import yaml
import io
import os

# Overwrite Protection Here
# Set to True - Will Overwrite Data already on disk
# Set to False - Will NOT Overwrite Data already on disk
overwrite = True


#Creates a list of directories with the appropriate files for analysis
flist = glob.glob('/mnt/bigdata/erosolow/Orion2/*/data.*.hdf5')
#Creating empty list for data sorting
flist_data = []
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
config_dir = '/home/sfielder/Documents/Clumps/'
data_dir = '/home/sfielder/Documents/Clumps/Output/'

#Make a list of all the config files found in  data_dir
flist_config = glob.glob(config_dir+'*.yaml')

for i in range(0,len(flist_config)):
    #Grabs the config file name here
    naming_string = flist_config[i].split("/")[-1].split(".")[0]
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
            with io.open(flist_config[i], 'r') as stream:
                data_loaded = yaml.load(stream)
        
            #Call main code here
            #Testing first file here
            for j in range(0,len(flist)):
                cc.analyzer(flist[j],
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
        with io.open(flist_config[i], 'r') as stream:
            data_loaded = yaml.load(stream)
        
        #Call main code here
        #Testing first file here
        for j in range(0,len(flist)):
            cc.analyzer(flist[j],
                    data_loaded['l'],
                    data_loaded['cmin'],
                    data_loaded['step'],
                    data_loaded['beta'],
                    data_loaded['clump_sizing'],
                    save_dir)
            