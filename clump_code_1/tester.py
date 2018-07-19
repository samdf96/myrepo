#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:12:57 2018

@author: sfielder
"""

import clump_code_1 as cc
import header_printer as hp
import glob
import yaml
import io
import os
from definitions import jComparisonPlotter
from definitions import FiducialPlotter

# Logging Info Here:

import logging
import logging.config

logging.config.fileConfig('logging.conf')

# create logger
logger = logging.getLogger('initialize')

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
overwrite_analyzer = True
overwrite_header = True
overwrite_timestep_plots = True
overwrite_fiducial_plots = True

logger.info("Overwrite Protection for Analyzer function has been set to: ",
            overwrite_analyzer)
logger.info("Overwrite Protection for Header Printer function has been set to: ",
            overwrite_header)
logger.info("Overwrite Protection for Timestep Plots function has been set to: ",
            overwrite_timestep_plots)
logger.info("Overwrite Protection for Fiducial Plots function has been set to: ",
            overwrite_fiducial_plots)

# =============================================================================
#INPUTS HERE

#Creates a list of directories with the appropriate files for analysis
flist = glob.glob('/Users/sfielder/Documents/Astro_Data/**/data.*.hdf5', recursive=True)

logger.info("Glob function has found the following to be sifted: ", flist)

#This is to filter out the timestamps that we want to analyze over
data_check_list = ['0060','0070','0080','0090','0100']

logger.info("Data Check List has been set to: ", data_check_list)

#This is where the config files are
tree_top_dir = '/Users/sfielder/Documents/Astro_Data/Clumps'
data_dir = '/Users/sfielder/Documents/Astro_Data/Clumps/Output/'

logger.info("tree_top_dir has been set to: ", tree_top_dir)
logger.info("data_dir has been set to: ", data_dir)

# =============================================================================
logger.info("Analysis Section Started.")
#Creating empty list for data sorting
flist_data = []
logger.debug("Filtering data by data_check_list_entry.")
for i in range(0,len(flist)):
    main_string = flist[i].split("/")
    out_string = main_string[-1].split(".")
    time_stamp = out_string[1]
    #This checks if timestamp
    if any(x in time_stamp for x in data_check_list):
        flist_data.append(flist[i])

#Sorting the Data by filename
flist_data.sort()

logger.info("Glob function has found the following to be analyzed: ",
            flist_data)
#Make a list of all the yaml files found in data_dir
logger.info("Finding All Config Files.")
flist_config_yaml = glob.glob(tree_top_dir+'**/*.yaml', recursive=True)
flist_config_yaml.sort()

logger.info("The following files will be inputs to the analyzer: ",
            flist_config_yaml)

carry_on_analyzer = False #Initialization Value - this is the default
logger.debug("Setting Initialization value for carry_on to: ",
             carry_on_analyzer)
for i in range(0,len(flist_config_yaml)):
    logger.debug("Currentl Working on Config File: ", flist_config_yaml[i])
    #Grabs the config file name here
    naming_string = flist_config_yaml[i].split("/")[-1].split(".")[0]
    #Creating String for Directory
    save_dir = data_dir + naming_string + '/'
    
    #Checking if Directory Exists
    if os.path.isdir(save_dir) == True:
        logger.info("Save Directory: ", save_dir, " has been detected as existing.")
        if overwrite_analyzer == True:
            logger.info("Overwrite for Analyzer has been set to TRUE.")
            logger.info("Carry On Value is being set to TRUE.")
            carry_on_analyzer = True
        else:
            logger.info("Overwrite for Analyzer has been set to FALSE.")
            logger.info("Carry On Value will remain as FALSE.")
    else:
        logger.info("Save Directory: ", save_dir, " has been detected as non-existent.")
        logger.info("Save Directory: ", save_dir, " is being created.")
        os.mkdir(save_dir)
        logger.info("Setting Carry On Value to TRUE.")
        carry_on_analyzer = True
        
    #If Carry_on_Analyzer has been set to true, then run the analysis.    
    if carry_on_analyzer == True:
        logger.info("Carry On Value has been detected as TRUE.")
        logger.info("Proceeding with Analysis.")
        #Importing Config File settings here
        with io.open(flist_config_yaml[i], 'r') as stream:
            data_loaded = yaml.load(stream)
            logger.info("Config File has been opened and settings extracted.")
        
        #Call main code here
        #Testing first file here
        for j in range(0,len(flist_data)):
            logger.info("Currently working on file: ", flist_data[j])
            logger.info("Invoking Analyzer function (cc).")
            cc.Analyzer(flist_data[j],
                    data_loaded['l'],
                    data_loaded['cmin'],
                    data_loaded['step'],
                    data_loaded['beta'],
                    data_loaded['clump_sizing'],
                    save_dir)
    else:
        logger.info("Carry On Value has been detected as FALSE.")
        logger.info("Skipping Analysis for current config file.")

logger.info("Analysis Section Complete.")
# =============================================================================
logger.info("Header Printer Section Started.")
if overwrite_header == True:
    logger.info("Overwrite for Header Printer has been set to TRUE.")
    # Call Header Printer Script to compute the .txt file needed for
    # summary of analysis
    # Grab the config_x directories, and puts them in a list
    flist_config_dir = glob.glob(data_dir + 'config_*')
    
    #Run Loop over all the config_x directories found that are in the list
    
    for i in range(0,len(flist_config_dir)):
        logger.info("Processing Summary File for directory: ", flist_config_dir[i])
        logger.info("Invoking HeaderPrinter function from header_printer.py")
        hp.HeaderPrinter(flist_config_dir[i])
else:
    logger.info("Overwrite for Header Printer has been set to FALSE.")
    logger.info("Skipping Header Printer Function.")

logger.info("Header Printer Section Completed.")
# =============================================================================
logger.info("Timestep Plots Section Started.")
if overwrite_timestep_plots == True:
    logger.info("Overwrite for Timestep Plots has been set to TRUE.")

    # Comparison Plots by specific timestep happens here:
    flist_plots = glob.glob(data_dir + '**/*.fits', recursive=True)
    flist_plots.sort()
    logger.info("Files to loop over: ", flist_plots)
    
    for i in range(0,len(flist_plots)):
        current_file = flist_plots[i]
        logger.info("Current File being worked on: ", current_file)
        logger.info("Invoking jComparisonPlotter function.")
        jComparisonPlotter(current_file)
else:
    logger.info("Overwrite for Timestep Plots has been set to FALSE.")
    logger.info("Skipping Timestep Plot Creation.")

logger.info("Timestep Plots Section Completed.")
# =============================================================================
logger.info("Fiducial Plots Section Started.")
# Comparison Plots over Fiducial Runs (by timestep) happens here:
carry_on_fiducial_plots = False #Initialize value - this is the default
logger.info("Carry On Value for Fiducial Plots is initialized as: ",
            carry_on_fiducial_plots)
flist_config = glob.glob(data_dir+'config_*')
flist_config.sort()
logger.info("Output config_x subdirectory list found to be: ", flist_config)
for k in range(0,len(flist_config)): #Adding Loop for config directories
    logger.info("Currently working on config: ", flist_config[k])
    #Write in os function to create appropiate directory for Fiducial Plots
    fid_dir = flist_config[k]+'/Fiducial_Plots/'
    if os.path.isdir(fid_dir) == True:
        logger.info("Save Directory: ", fid_dir, " has been detected to exist.")
        if overwrite_fiducial_plots==True:
            logger.info("Overwrite for Fiducial Plots has been set to TRUE.")
            logger.info("Carry On Value is being set to TRUE.")
            carry_on_fiducial_plots = True
        else:
            logger.info("Overwrite for Fiducial Plots has been set to FALSE.")
            logger.info("Carry On Value will remain as FALSE.")
    else:
        logger.info("Save Directory: ", fid_dir, " has been detected as non-existent.")
        logger.info("Creating Save Directory.")
        os.mkdir(fid_dir)
        logger.info("Setting Carry On Value to TRUE.")
        carry_on_fiducial_plots = True

    #If Carry On Value is True then continue with plotting
    if carry_on_fiducial_plots == True:
        logger.info("Carry On Value has been detected as TRUE.")
        logger.info("Proceeding with Analysis.")
        logger.info("Looping timesteps inputted by data_check_list.")
        for i in range(0,len(data_check_list)):
            flist = glob.glob(flist_config[k]+'/**/*'+data_check_list[i]+'*.fits',
                          recursive=True)
            flist.sort()
            data_check_list_print = data_check_list[i]
            #Calling Main Function Here
            logger.info("Current Timestep being worked on: ", data_check_list_print)
            logger.info("Invoking FiducialPlotter function.")
            FiducialPlotter(flist,fid_dir,data_check_list_print)

logger.info("Fiducial Plots Section Completed.") 
# =============================================================================