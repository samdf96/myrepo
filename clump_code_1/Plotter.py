#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 10:57:36 2018

@author: sfielder

This will take an input fits file, and compute the appropriate j/j plots
"""


from astropy.io import fits
import glob
import matplotlib.pyplot as plt

## Import Definitions from definitions.py here:
from definitions import data_grabber
from definitions import j_comp_plotter


#Main Tree Start
input_dir = '/Users/sfielder/Documents/Astro_Data/'
config_file_string = 'config_1/' #Add slash after file name



#Here is where we define the specific config folder which houses the config files
flist = glob.glob(input_dir + config_file_string + '**/**/*.fits')
flist.sort() #Sorts the Config files by name for easier readability

## Now we have all the fits files in a list: flist


for i in range(0,1):
    
    current_file = flist[i]
    
    
    hdu = fits.open(current_file)
        
    #Grabbing BinHDUTable object here, should always be second object
    hdu_table = hdu[1]
    data = hdu_table.data
    
    #Grabbing All the Data from the hdu_table here, makes dictionary    
    data_dict = data_grabber(data)
    
   
    ## Making First Figure Here for x-LOS
    
    x = data_dict['imp_x_los']
    y1 = data_dict['actual_tot']
    y2 = data_dict['actual_par_yz']
    
    x_str = 'Implied Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    y_str = 'Actual Specific Angular Momentum [$kg \ m^2 \ s^{-1}$]'
    title_str = 'Comparison for X Line-of-Sight'
    
    #Calling Function Here for Plotting
    fig = j_comp_plotter(x,
                         y1,
                         y2,
                         x_str,
                         y_str,
                         title_str) 
    plt.tight_layout()
    fig.savefig('Test_Fig.pdf', bbox_inches='tight')
    
    

    