#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 15:26:55 2018

@author: sfielder
"""

import clump_code_1 as cc
import glob
import yaml
import io
import os




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
    
    #Checking if Directory Exists, if so, it skips analysis completely
    #to avoid overwritting data.
    if os.path.isdir(save_dir) == True:
            print("Warning!!! Directory: " +
                  save_dir +
                  "is detected as a valid directory. Skipping Analysis.")
    else:
        #If Directory is not present, then it gets created here and the code
        #is ran through the analyzer
        os.mkdir(save_dir)
        """
        #Importing Config File settings here
        with io.open(flist_config[i], 'r') as stream:
            data_loaded = yaml.load(stream)
        
        #Call main code here
        for j in range(0,len(flist)):
            cc.analyzer(flist[j],
                    data_loaded['l'],
                    data_loaded['cmin'],
                    data_loaded['step'],
                    data_loaded['beta'],
                    data_loaded['clump_sizing'],
                    save_dir)
        """
"""
for i in range(0,len(flist)):
    for j in range(0,len(flist_config)):
        naming_string = str(j)
        
        #Text String for Output Directory
        save_dir = data_dir+'Config_'+ naming_string + '/'
        
        if os.path.isdir(save_dir) == True:
            print("Warning!!! Directory: " +
                  save_dir +
                  "is detected as a valid directory. Skipping Analysis.")
        else:
            #If Directory is not present, then it gets created here and the code
            #is ran through the analyzer
            os.mkdir(save_dir)
            
            #Importing Config File settings here
            with io.open(flist_config[j], 'r') as stream:
                data_loaded = yaml.load(stream)
            
            #Call main code here
            cc.analyzer(flist[i],
                    data_loaded['l'],
                    data_loaded['cmin'],
                    data_loaded['step'],
                    data_loaded['beta'],
                    data_loaded['clump_sizing'],
                    save_dir)
"""