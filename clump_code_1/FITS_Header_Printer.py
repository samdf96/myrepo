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
    - input_dir : string
        - Top level of the tree where the config files are
    - config_file_string : string
        - Config file name

Returns:
    - txt file with all the comments about the fits files that are in the 
        config_dir directory
"""

from astropy.io import fits
import glob


def Header_Printer(input_dir, config_file_string):
    
    '''
    FOR LOCAL TESTING ONLY: DELETE AFTER TESTING
    #Main Tree Start
    input_dir = '/Users/sfielder/Documents/Astro_Data/'
    config_file_string = 'config_1/' #Add slash after file name
    '''
    
    
    config_dir = input_dir + config_file_string[:-1]
    
    #Here is where we define the specific config folder which houses the config files
    flist = glob.glob(input_dir + config_file_string + '**/**/*.fits')
    flist.sort() #Sorts the Config files by name for easier readability
    
    #File Writing Here - Not using conventional ending with slash for next line
    
    #Creating text file name here
    config_string = config_dir.split("/")[-1]
    text_filename = config_dir + '/' + config_string + '_Header_Info.txt'
    
    with open(text_filename, 'w') as txt:   #Opens the file to write
        
        #Generic Print Statements for TXT File
        print('Below presents the Simulation Input parameters ' +
              'and COMMENT(S) for the FITS Files from the Directory: '+
              config_string, file=txt)
        print('', file=txt)
    
        for i in range(0,len(flist)):    #Loops over all fits files here
        
            current_file = flist[i]
            
            #Opening and closing a first fits file of the config_dir to get input
            #parameters, these write to the txt file opened above
            if i == 0:
                main_hdu = fits.open(current_file)
                param_hdu = main_hdu[1]
                param_header = param_hdu.header
                print('##############################', file=txt)
                print('The following are input parameters for the data in directory: ' +
                      config_string, file=txt)
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
                
#                file_string = current_file.split()
#                fits_string_true = ['fits' in k for k in file_string]
#                fits_string_id = [j for j, x in enumerate(fits_string_true) if x]
#                fits_string = file_string[fits_string_id[0]]
                
                header = hdu_table.header
                #header_keys = list(header.keys())
                #print(header_keys)
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
                
    return()
    
