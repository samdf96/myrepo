#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:12:57 2018

@author: sfielder
"""


import analyzer as analyzer
import header_printer as hp
import glob
import yaml
import io
import os
from definitions import jComparisonPlotter
from definitions import TimestepPlotter

# Logging Info Here:

import logging
import logging.config

#Twillo Account and SMS integration.
messaging = False
from twilio.rest import Client

with io.open('/home/sfielder/Documents/TwilloAccountInfo.yaml', 'r') as TwillAcnt:
    TwillAcntDict = yaml.load(TwillAcnt)

account_sid = TwillAcntDict['Acc_SID']
auth_token = TwillAcntDict['Acc_TKN']
twilio_phone_number = TwillAcntDict['T_N']
my_phone_number = TwillAcntDict['T_M']

client = Client(account_sid, auth_token)

if messaging == True:
    message_start = client.messages.create(to=my_phone_number, from_=twilio_phone_number,
                                           body="Python Script has Started.")


#INPUTS HERE

#Creates a list of directories with the appropriate files for analysis
# THIS WILL NEED TO BE CHANGED FOR THE NEWER DESIGN SIMULATIONS

#This batches input filters the search criteria to only look for `batches` simulation directories
batches = 'Design'
#flist = glob.glob('/mnt/bigdata/erosolow/Orion2/*'+batches+'*/data.*.hdf5')

#Call ALL Found Files in input directory
flist = glob.glob('/mnt/bigdata/erosolow/Orion2/**/data.*.hdf5')

#This is to filter out the timestamps that we want to analyze over
data_check_list = ['0060','0070','0080','0090','0100']


#This is where the config files are
tree_top_dir = '/Users/sfielder/Documents/Astro_Data/'
data_dir = '/Users/sfielder/Documents/Astro_Data/Fiducial00'

#Load CONFIG FILE HERE
logging.config.fileConfig('logging.conf', defaults={'logfilename': data_dir+'output.log'})

#Set the Batch to do here:
#Batches will be by Simulation Directory General Name (ex. Fiducial or Design...)


# create logger
logger = logging.getLogger('initialize')

logger.info("Batches search value set to: %s", batches)
logger.info("Glob function has found the following to be sifted: %s", flist)
logger.info("Data Check List has been set to: %s", data_check_list)
logger.info("tree_top_dir has been set to: %s", tree_top_dir)
logger.info("data_dir has been set to: %s", data_dir)

"""
Overwrite Protection Here:
    - header and timestep_plots have absolute status and will not even
      start looking for files if value below is set to FALSE.
    
    - analyzer and simulation_plots have conditional status, and will look for
      for already created directories and will skip those that are found
      if the value below is set to FALSE.
    
    - Setting any of the following values to TRUE will overwrite files
      even if files and directories are found by the code to exist.
"""

overwrite_analyzer = False
overwrite_header = False
overwrite_timestep_plots = False
overwrite_simulation_plots = False

logger.info("Overwrite Protection for Analyzer function has been set to: %s",
            overwrite_analyzer)
logger.info("Overwrite Protection for Header Printer function has been set to: %s",
            overwrite_header)
logger.info("Overwrite Protection for Timestep Plots function has been set to: %s",
            overwrite_timestep_plots)
logger.info("Overwrite Protection for Simulation Plots function has been set to: %s",
            overwrite_simulation_plots)

# =============================================================================
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

logger.info("Glob function has found the following to be analyzed: %s",
            flist_data)
#Make a list of all the yaml files found in data_dir
logger.info("Finding All Config Files.")
flist_config_yaml = glob.glob(tree_top_dir+'*.yaml')
flist_config_yaml.sort()

logger.info("The following files will be inputs to the analyzer: %s",
            flist_config_yaml)

carry_on_analyzer = False #Initialization Value - this is the default
logger.debug("Setting Initialization value for carry_on to: %s",
             carry_on_analyzer)
for i in range(0,len(flist_config_yaml)):
    logger.debug("Currentl Working on Config File: %s", flist_config_yaml[i])
    #Grabs the config file name here
    naming_string = flist_config_yaml[i].split("/")[-1].split(".")[0]
    #Creating String for Directory
    save_dir = data_dir + naming_string + '/'
    
    #Checking if Directory Exists
    if os.path.isdir(save_dir) == True:
        logger.info("Save Directory: %s has been detected as existing.", save_dir)
        if overwrite_analyzer == True:
            logger.info("Overwrite for Analyzer has been set to TRUE.")
            logger.info("Carry On Value is being set to TRUE.")
            carry_on_analyzer = True
        else:
            logger.info("Overwrite for Analyzer has been set to FALSE.")
            logger.info("Carry On Value will remain as FALSE.")
    else:
        logger.info("Save Directory: %s , has been detected as non-existent.", save_dir)
        logger.info("Save Directory: %s , is being created.", save_dir)
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
            logger.info("Currently working on file: %s", flist_data[j])
            logger.info("Invoking Analyzer function (cc).")
            analyzer.Analyzer(flist_data[j],
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

