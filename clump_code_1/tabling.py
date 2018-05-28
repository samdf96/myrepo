#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 28 12:11:25 2018

@author: sfielder
"""

from astropy.io import fits
from astropy.table import Table

main_table = Table.read('data_clump.fits')


#Import Data file here:
hdu_main = fits.open('data_clump.fits')

hdu_table_data = hdu_main[1].data
cols = hdu_main[1].columns

master_names = []
for i in range(0,len(cols.names)):
    master_names.append(cols.names[i])

#Creation of quantities:
master_dict_data = {}
master_dict_units = {}

for i in range(0,len(cols.names)):
    master_dict_data[str(master_names[i])] = hdu_table_data.field(str(master_names[i]))

