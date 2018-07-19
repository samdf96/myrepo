#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 09:25:24 2018

@author: sfielder



FITS File print header: Creates .txt with all the info for the .fits files
computed by the analyzer function. Only works with one config directory input

Since only one config file at a time will be processed this way, then we can
simply put all the config file info at the very start, things to do with the 
input parameters for the analyzer, and then grab all the comments for each fits
file and output them into the txt file that this file generates

Inputs:
    - config_list: string
        - path to config_x file found in initilize.py (see file for details)

Returns:
    - txt file with all the comments about the fits files that are in the 
        config_x directory (see initialize.py file for tree structure)
"""

from astropy.io import fits
import glob
import logging

module_logger = logging.getLogger("initialize.header_printer")

def HeaderPrinter(config_input):
    """
    See Above for details about this function.
    """
    logger = logging.getLogger("initialize.header_printer.HeaderPrinter")

    #Making input string into valid directory with slash ending
    config_dir = config_input + '/'
    
    #Here is where we define the specific config folder which houses the config files
    flist = glob.glob(config_dir + '**/**/*.fits')
    flist.sort() #Sorts the Config files by name for easier readability
    
    logger.info("Fits files found for function: ", flist)
    #Creating text file name here
    text_filename = config_dir + 'Header_Info.txt' 
    logger.info("Save Directory set as: ", text_filename)
    
    #Opens the file to write - Will always overwrite data due to '+' argument
    with open(text_filename, 'w+') as txt:
        logger.info("text_filename has been opened.")
        #Generic Print Statements for TXT File
        print('Below presents the Simulation Input parameters ' +
              'and COMMENT(S) for the FITS Files from the Directory: '+
              config_dir, file=txt)
        print('', file=txt)
    
        for i in range(0,len(flist)):    #Loops over all fits files here
            logger.info("Currently working on FITS file: ", flist[i])
            current_file = flist[i]
            
            #Opening and closing a first fits file of the config_dir to get input
            #parameters, these write to the txt file opened above
            if i == 0:
                logger.ingo("First flist object detected, writing extra info.")
                main_hdu = fits.open(current_file)
                param_hdu = main_hdu[1]
                param_header = param_hdu.header
                print('##############################', file=txt)
                print('The following are input parameters for the data in directory: ' +
                      config_dir, file=txt)
                print('Minimum Clump Number = ',
                      param_header['MINCLMP'],
                      file=txt)
                print('Step Size = ',
                      param_header['STEP'],
                      file=txt)
                print('Beta = ',
                      param_header['BETA'],
                      file=txt)
                print('Length (pc) = ',
                      param_header['LENGTH'],
                      file=txt)
                print('CMIN (g/cm^3) = ',
                      param_header['CMIN'],
                      file=txt)
                print('##############################', file=txt)
                print('', file=txt)
                #Close the file here, just being safe, output did not change when added
                main_hdu.close()
    
            with fits.open(current_file) as hdu:
                hdu_table = hdu[1]
                #Extracting header here                
                header = hdu_table.header
                print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', file=txt)
                print('Comments for FITS File: ', current_file, file=txt)
                while True:
                    try:
                        print(header['COMMENT'], file=txt)
                    except KeyError:
                        print('No Comments for this FITS File', file=txt)
                        pass
                    print('', file=txt)
                    break
                
                #Print Ending Statements for File Analysis
                clump_tot = hdu_table.columns['Clump Number'].array[-1]
                print('Total Number of Clumps Found: ', clump_tot, file=txt)
                print('', file=txt)
    o
    logger.info("HeaderPrinter has been run successfully.")     
    return()
    
