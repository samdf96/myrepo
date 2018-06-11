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

# =============================================================================
# OLD CODE FOR LOCAL MACBOOK AIR DIRECTORIES
# #Directory to save all data and config files
# data_dir = '/Users/sfielder/Documents/Astro_Data/'
# flist = glob.glob(data_dir+'*.hdf5')
# save_dir_fits = data_dir+'FITS_Files/'
# config_dir = data_dir+'configs/'
# =============================================================================


#Finds all files with the .hdf5 extension in the bigdata directory
flist = glob.glob('/home/sfielder/bigdata/**/*.hdf5')
data_dir = '/home/sfielder/Documents/Clumps/'
save_dir = data_dir+'FITS_YAML/'

with io.open(data_dir+"config_1.yaml", 'r') as stream:
    data_loaded = yaml.load(stream)

for i in range(0,len(flist)):
    cc.analyzer(flist[i],
                data_loaded['l'],
                data_loaded['cmin'],
                data_loaded['step'],
                data_loaded['beta'],
                data_loaded['clump_sizing'],
                save_dir)
