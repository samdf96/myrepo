#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 15:40:39 2018

@author: sfielder
"""

"""
Fiducial Plotter Function

Input:
    - data_dir: string
        - Will be the otuput directory of the data.
    - data_check_list: list of strings
        - Will be the timestamps that we are looking for, are going to be the 
        same ones as the data is sorted through.
"""

#Setting Local inputs for now:
data_dir = 'C:/Users/Samuel Fielder/Documents/Astro/'
data_check_list = ['0060','0070','0080','0090','0100']

import glob
from astropy.io import fits
import matplotlib.pylab as plt
import numpy as np

#From Definitions Here:
from definitions import data_grabber


#Find all data sets with one time stamp here - will loop
for i in range(0,1): #len(data_check_list)
    flist = glob.glob(data_dir+'**/*'+data_check_list[i]+'.fits', recursive=True)
    flist.sort()

def FiducialPlotter(flist, equal_axis, percentage):
    
    for i in range(0,len(flist)):
        current_file = flist[i]
        hdu = fits.open(current_file)
        #Print Statements for File being worked on:
        filename_printing = current_file.split("/")[-2]
        print("Current TimeStamp being worked on: {}".format(filename_printing))
        
        
        #Grabbing BinHDUTable object here, should always be second object
        hdu_table = hdu[1]
        data = hdu_table.data #Grabs data stored in the table -> FITS REC   
        d = data_grabber(data) 
        hdu.close()
        
        #Making the Output Directory for current file:
        save_dir_list = current_file.split("/")[:-3]
        save_dir = '/'.join(save_dir_list)
        
        #Make Masked arrays for data that is nan valued
        """
        Will use the act_tot_mask array for master mask condition
        Below creates masked data from the dictionary entries using the condition 
        stated above to parse through the data and extract non-nan values
        """
        
        imp_x_los_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                            d['imp_x_los'])
        imp_y_los_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                            d['imp_y_los'])
        imp_z_los_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                            d['imp_z_los'])
        act_tot_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                          d['actual_tot'])
        act_par_xy_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                             d['actual_par_xy'])
        act_par_xz_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                             d['actual_par_xz'])
        act_par_yz_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                             d['actual_par_yz'])
        mass_mask = np.ma.masked_where(np.isnan(d['actual_tot']) == True,
                                       d['mass'])
        
        #Compress all the masked array to array with correct data
        imp_x_los_mask_compressed = np.ma.compressed(imp_x_los_mask)
        imp_y_los_mask_compressed = np.ma.compressed(imp_y_los_mask)
        imp_z_los_mask_compressed = np.ma.compressed(imp_z_los_mask)
        act_tot_mask_compressed = np.ma.compressed(act_tot_mask)
        act_par_xy_mask_compressed = np.ma.compressed(act_par_xy_mask)
        act_par_xz_mask_compressed = np.ma.compressed(act_par_xz_mask)
        act_par_yz_mask_compressed = np.ma.compressed(act_par_yz_mask)
        mass_mask_compressed = np.ma.compressed(mass_mask)
        
        
        # =========================================================================
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    X LOS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        axis_str = 'X'
        if equal_axis == False:
            filename = 'fiducial_comp_X_los.pdf'
        else:
            filename = 'fiducial_comp_X_los_eq.pdf'
        
        savefig_dir = save_dir + '/' + filename
        
        #Calling Function Here for Plotting
        fig = j_comp_plotter(imp_x_los_mask_compressed,
                             act_tot_mask_compressed,
                             act_par_yz_mask_compressed,
                             axis_str,
                             equal_axis,
                             percentage)
        plt.tight_layout()
        fig.savefig(savefig_dir, bbox_inches='tight')
        plt.close(fig)
        
        #Colormap Variant Here
        if equal_axis == False:
            filename = 'j_comp_X_los_colormapped.pdf'
        else:
            filename = 'j_comp_X_los_colormapped_eq.pdf'
    
        savefig_dir = save_dir + '/' + filename
        
        #Calling Function Here for Plotting
        fig = j_comp_plotter_colormap(imp_x_los_mask_compressed,
                             act_tot_mask_compressed,
                             act_par_yz_mask_compressed,
                             mass_mask_compressed,
                             axis_str,
                             equal_axis,
                             percentage)
        plt.tight_layout()
        fig.canvas.draw()
        fig.savefig(savefig_dir, bbox_inches='tight')
        plt.close(fig)
        # =========================================================================
        
        # =========================================================================
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Y LOS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        axis_str = 'Y'
        if equal_axis == False:
            filename = 'j_comp_Y_los.pdf'
        else:
            filename = 'j_comp_Y_los_eq.pdf'
    
        savefig_dir = save_dir + '/' + filename
        
        #Calling Function Here for Plotting
        fig = j_comp_plotter(imp_y_los_mask_compressed,
                             act_tot_mask_compressed,
                             act_par_xz_mask_compressed,
                             axis_str,
                             equal_axis,
                             percentage)
        plt.tight_layout()
        fig.savefig(savefig_dir, bbox_inches='tight')
        plt.close(fig)
        
        #Colormap Variant Here
        if equal_axis == False:
            filename = 'j_comp_Y_los_colormapped.pdf'
        else:
            filename = 'j_comp_Y_los_colormapped_eq.pdf'
    
        savefig_dir = save_dir + '/' + filename
        
        #Calling Function Here for Plotting
        fig = j_comp_plotter_colormap(imp_y_los_mask_compressed,
                             act_tot_mask_compressed,
                             act_par_xz_mask_compressed,
                             mass_mask_compressed,
                             axis_str,
                             equal_axis,
                             percentage)
        plt.tight_layout()
        fig.canvas.draw()
        fig.savefig(savefig_dir, bbox_inches='tight')
        plt.close(fig)
        # =========================================================================
        
        # =========================================================================
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Z LOS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        axis_str = 'Z'
        if equal_axis == False:
            filename = 'j_comp_Z_los.pdf'
        else:
            filename = 'j_comp_Z_los_eq.pdf'
    
        savefig_dir = save_dir + '/' + filename
        
        #Calling Function Here for Plotting
        fig = j_comp_plotter(imp_z_los_mask_compressed,
                             act_tot_mask_compressed,
                             act_par_xy_mask_compressed,
                             axis_str,
                             equal_axis,
                             percentage)
        plt.tight_layout()
        fig.savefig(savefig_dir, bbox_inches='tight')
        plt.close(fig)
        
        #Colormap Variant Here
        if equal_axis == False:
            filename = 'j_comp_Z_los_colormapped.pdf'
        else:
            filename = 'j_comp_Z_los_colormapped_eq.pdf'
    
        savefig_dir = save_dir + '/' + filename
        
        #Calling Function Here for Plotting
        fig = j_comp_plotter_colormap(imp_z_los_mask_compressed,
                             act_tot_mask_compressed,
                             act_par_xy_mask_compressed,
                             mass_mask_compressed,
                             axis_str,
                             equal_axis,
                             percentage)
        plt.tight_layout()
        fig.canvas.draw()
        fig.savefig(savefig_dir, bbox_inches='tight')
        plt.close(fig)
        # =========================================================================
    
    
        return()


