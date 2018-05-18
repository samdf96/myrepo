#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 18 14:49:42 2018

@author: sfielder


This script will take incoming np.array archive file exported by the
main script and attach units as well as table them, using astropy
"""
#%%
from astropy import units as u
import numpy as np


data = np.load('clump_data.npz')

#Load data into a dictionary so it can be loaded arbitrarily.
print(data.files)
data_dict = dict(data)

#Attach Units to certain data in the dictionary

clump_number = data_dict['clump_number']
clump_number *= u.dimensionless_unscaled

com = data_dict['center_of_mass']
center_of_mass_x = []
center_of_mass_y = []
center_of_mass_z = []
for i in range(0,len(com)):
    center_of_mass_x.append(com[i,1,0])
    center_of_mass_y.append(com[i,1,1])
    center_of_mass_z.append(com[i,1,2])






