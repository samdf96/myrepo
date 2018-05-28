#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 15:26:55 2018

@author: sfielder
"""

import clump_code_1 as cc
import glob


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# To be used with YAML file which we will create later on
#import sys
#filename = sys.argv[1]


#Defining Global Variables Here - Move to File After
l=10	#Length of Original Data Set
cmin = 5e-21	#Minimum Density Threshold for Clump Finding
#cmax = 5e-20	#Maximum Density Threshold for Clump Finding
step = 100	 #Step-size multiplier for Clump Finding
beta = 1        # For Implied Angular Momentum Calculation.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

flist = glob.glob('/Users/sfielder/Documents/Astro_AngMom/*.hdf5')

#%%
for i in range(0,len(flist)):
    cc.analyzer(flist[i],l,cmin,step,beta)