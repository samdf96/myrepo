#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 19:53:52 2018

@author: sfielder
"""


#Importing packages here
import yt
import numpy as np
import astropy.units as u
from astropy.io import fits
import os

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Section for importing all the definitions from defitions.py file
from definitions import octant_split
from definitions import master_clump_maker
from definitions import clump_finder
from definitions import center_of_mass
from definitions import bounding_box
from definitions import clump_box
from definitions import velocity_array
from definitions import velocity_array_reducer
from definitions import array_flattener
#Add in if needed
#from definitions import myplane
#from definitions import plane_fit_visualization
from definitions import plane_fit
from definitions import gradient
from definitions import angular_momentum_implied
from definitions import angular_momentum_actual
from definitions import proj_creator

#Section for Importing Exceptions classes
from exceptions import YTErrorValue, YTErrorReshape
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

def analyzer(filename,l,cmin,step,beta,clump_sizing,save_dir_fits):
    '''
    Takes in certain global parameters, and a simulation datafile and computes
    all the relevant values for angular momentum, implied and actual.
    
    Arguments:
    ---------- 
    filename: string
        must point to a data-simulation file (mostly .hdf5 files)
        needs to be the full directory string of the filename
    l: float
        length of data simulation in pc along one side of cube
    cmin: float
        Minimum Density Threshold for Clump Finding
    step: float
        Step-size multiplier for Clump Finding
    beta: float
        For Implied Angular Momentum Calculation
    clump_sizing: float
        Minimum Value for Clump Finding
    save_dir: string
        Directory imput to save all FITS files
        
    Returns:
    -------- 
    FITS File contatining all the data computed.
        - Has all clump data
        - COMMENTS in the header for any clumps that returned errors
    
    
    '''
    
    print('Simluation File Currently Working On: ',filename)
    
    #Loads data into File
    ds = yt.load(filename)
    
    #Creatingt Text Strings for Output here
    
    #Splitting the Whole Directory String by / 
    main_string = filename.split("/")
    #Finding where the list matches Fiducial
    fid_string_true = ['Fiducial' in k for k in main_string]
    #Grabbing that indice
    fid_string_id = [i for i, x in enumerate(fid_string_true) if x]
    #Making a string that is called the Directory found above ex. Fiducial00
    fid_str = main_string[fid_string_id[0]] #Should only have one entry
    #Grabs last component of filename
    out_string = main_string[-1].split(".") #Splits around periods
    time_stamp = out_string[1] #Ex. 0060
    
    #Creating Sub Directory for individual Data Simulation File Timestamp
    
    save_dir_specific = save_dir_fits + time_stamp + "/"
    
    
    master_dist_data = int(ds.domain_dimensions[0])
    
    #Creating Directory to store picture files and FITS files in.
    #First detects if directory is there, if so, it deletes and rebuilds that
    #directory
    if os.path.isdir(save_dir_specific) == True:
        print("Warning!!! Directory: " +
              save_dir_specific +
              "is detected as a valid directory." +
              "FITS Files will be overwritten.")
    else:
        os.mkdir(save_dir_specific)
    
    #Creates a Data Object containing all the Simulation Data
    ad = ds.all_data()
    
    #Splits the Data into Octants
    octant = octant_split(ad,ds,l)
    
    #Grabs each of the octants and runs them through the Clump Finding algorithm
    clumps = [] #Defining Empty List for loop
    for i in range(0,len(octant)):
        master_clump_main = master_clump_maker(octant[i])
        cmax = octant[i]["gas", "density"].max()
        print("Now Finding clumps for Octant:",i)
        lc = clump_finder(master_clump_main,clump_sizing,cmin,cmax,step)
        for j in range(0,len(lc)):
            clumps.append(lc[j])
    
    # =========================================================================
    # Making empty arrays to store the data, variable to allow for different
    # size of clumps list. Actual array (1,3) and (3,2) will not change as
    # that is how they are written with no difference between clump objects.
    # =========================================================================
    
    #Creating Arrays for later export (see bottom of script for this step)
    com = np.zeros((len(clumps),1,3))
    bregion = np.zeros((len(clumps),3,2))
    
    #Defining Unit for Angular Momentum
    ang_mom_unit = u.kg * (u.m**2) / u.s
    per_sec_unit = u.s**(-1)
    
    # Creates lists for center_of_mass of clump list -> stores in com
    # Creates list for bounding boxes of clump list -> stores in bregion
    for i in range(0,len(clumps)):
    
        com[i] = center_of_mass(clumps[i])
        bregion[i] = bounding_box(clumps[i])

    # =========================================================================
    # CHECK IN POINT: What we have so far here in the script
    #
    #     com: arrays that contain x,y,z center of mass values for clumps
    #     bregion: x,y,z min/max value to build the boxes for clumps
    # ========================================================================        
    
    data_object_clump = [] #Setting Empty list for loop
    mass_clump_g = []       #Setting Empty list for loop
    dist_span_cm = np.zeros((len(bregion),3,1))
    #   Setting Empty Arrays to store individual coordinates for COM of Clumps
    #   and distance spanned over each axis for each clump
    com_x = np.zeros((len(clumps),1))
    com_y = np.zeros((len(clumps),1))
    com_z = np.zeros((len(clumps),1))
    dist_span_cm_x = np.zeros((len(clumps),1))
    dist_span_cm_y = np.zeros((len(clumps),1))
    dist_span_cm_z = np.zeros((len(clumps),1))
        
    for i in range(0,len(bregion)):
        print('Creating Clump Number:',i)
        data_object_clump.append(clump_box(ds,bregion[i]))
        
        mass_clump_g.append(np.array(data_object_clump[i].quantities.total_mass())[0])
        
        dist_x = bregion[i,0,1] - bregion[i,0,0]
        dist_y = bregion[i,1,1] - bregion[i,1,0]
        dist_z = bregion[i,2,1] - bregion[i,2,0]
        
        dist_tot = np.array([[dist_x],[dist_y],[dist_z]])
        
        # Already in array format for export later
        dist_span_cm[i] = dist_tot
        
        #Store individual clump COM coordinates
        com_x[i,0] = com[i,0,0]
        com_y[i,0] = com[i,0,1]
        com_z[i,0] = com[i,0,2]
        
        #Store individual clump spanning distances for each axis
        dist_span_cm_x[i,0] = dist_tot[0,0]
        dist_span_cm_y[i,0] = dist_tot[1,0]
        dist_span_cm_z[i,0] = dist_tot[2,0]
    
    
    
    # =========================================================================
    #Creating Tuples for Center of Mass along different LOS directions
    com_plotting = np.concatenate((com_x,com_y,com_z), axis=1)
    
    x_los_com = np.concatenate((com_y,com_z), axis=1)
    y_los_com = np.concatenate((com_x,com_z), axis=1)
    z_los_com = np.concatenate((com_x,com_y), axis=1)
    # Creatin of projection plots with markers for clump center of masses
    proj_creator(ds,
                 ad,
                 com_plotting,
                 x_los_com,
                 y_los_com,
                 z_los_com,
                 save_dir_specific,
                 fid_str,
                 time_stamp)
    
    print('PLOT SAVING COMPLETED')
    #Turn into Array for later export (see bottom of script for this step)
    mass_clump_g = np.array(mass_clump_g)
    
    
    
    #Creation of Columns for Astropy FITS for all quantities computed above
    q_com_x = fits.Column(name='Center of Mass Coordinate (x)',
                          format='D',unit='cm',array=com_x)
    q_com_y = fits.Column(name='Center of Mass Coordinate (y)',
                          format='D',unit='cm',array=com_y)
    q_com_z = fits.Column(name='Center of Mass Coordinate (z)',
                          format='D',unit='cm',array=com_z)
    
    q_mass_clump = fits.Column(name='Mass',format='D',unit='g',
                               array=mass_clump_g)
    
    q_dist_span_x = fits.Column(name='Length (x axis)',
                                format='D',
                                unit='cm',
                                array=dist_span_cm_x)
    q_dist_span_y = fits.Column(name='Length (y axis)',
                                format='D',
                                unit='cm',
                                array=dist_span_cm_y)
    q_dist_span_z = fits.Column(name='Length (z axis)',
                                format='D',
                                unit='cm',
                                array=dist_span_cm_z)
    
    
    
    
    # =========================================================================
    # CHECK IN POINT: What we have so far here in the script
    #
    #   data_object_clump: list
    #       a list of YTRegions that represent the clumps found
    #   mass_clump_g: array
    #       an array representing the total mass of each clump
    #   dist_span_cm: array
    #       an array with each superset array having all axis distances
    #   dist_span_cm_(x,y,z): array
    #       an array with just (x,y,z)-spanning distance for all clumps
    #   com_(x,y,z): array
    #       an array with just (x,y,z)-coordinates for all clumps
    # =========================================================================
    
    #Defining Empty Lists to store data into, will convert into arrays later on
    grad_x = []
    grad_y = []
    grad_z = []
    am_implied_x = []
    am_implied_y = []
    am_implied_z = []
    am_actual_total = []
    am_actual_partial_xy = []
    am_actual_partial_xz = []
    am_actual_partial_yz = []
    err_string = [] #Used for tracking errors in functions
    
    print('All Clumps Created Successfully')
    
    for i in range(0,len(bregion)):
        clump = data_object_clump[i]
        
        print('Working on Analysis for Clump Number:',i)

        # =====================================================================
        # Computing Integrated Velocity Arrays
        arr_z, vz = velocity_array(clump,'velocity_z','z',master_dist_data,l)            
        arr_x, vx = velocity_array(clump,'velocity_x','x',master_dist_data,l)
        arr_y, vy = velocity_array(clump,'velocity_y','y',master_dist_data,l)
        # =====================================================================
        
        # =====================================================================
        # Computing Reduced Velocity Arrays
        
        # This while loop raises errors if broken values from output of 
        # definitions file is not zero.
        while True:
            try:
                arr_x_red, vx_py, vx_pz, broken = velocity_array_reducer(arr_x,
                                                         vx,
                                                         'x',
                                                         master_dist_data)
                #Check if broken statement is triggered
                if broken == 1:
                    raise YTErrorValue
                if broken == 2:
                    raise YTErrorReshape
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                arr_y_red, vy_px, vy_pz, broken = velocity_array_reducer(arr_y,
                                                         vy,
                                                         'y',
                                                         master_dist_data)
                #Check if broken statement is triggered
                if broken == 1:
                    raise YTErrorValue
                if broken == 2:
                    raise YTErrorReshape
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                arr_z_red, vz_px, vz_py, broken = velocity_array_reducer(arr_z,
                                                         vz,
                                                         'z',
                                                         master_dist_data)
                #Check if broken statement is triggered
                if broken == 1:
                    raise YTErrorValue
                if broken == 2:
                    raise YTErrorReshape
                
                # =============================================================
                # Flattening Arrays for Plane Fitting Process
                arr_x_red_flat = array_flattener(arr_x_red)
                arr_y_red_flat = array_flattener(arr_y_red)
                arr_z_red_flat = array_flattener(arr_z_red)
                
                vx_py_flat = array_flattener(vx_py)
                vx_pz_flat = array_flattener(vx_pz)
                vy_px_flat = array_flattener(vy_px)
                vy_pz_flat = array_flattener(vy_pz)
                vz_px_flat = array_flattener(vz_px)
                vz_py_flat = array_flattener(vz_py)
                # =============================================================
                
                #Plane Fitting Process
                result_x = plane_fit(vx_py_flat,vx_pz_flat,arr_x_red_flat)
                result_y = plane_fit(vy_px_flat,vy_pz_flat,arr_y_red_flat)
                result_z = plane_fit(vz_px_flat,vz_py_flat,arr_z_red_flat)
                
                #Computing Gradients from Results
                gradient_x = gradient(result_x)
                gradient_y = gradient(result_y)
                gradient_z = gradient(result_z)
            
                #Computing Specific Angular Momentum
                # BETA VALUE INPUT TO DEFINITION WILL DEFAULT TO ONE
                ang_mom_implied_x = angular_momentum_implied(gradient_x,
                                                              dist_span_cm[i,1,0],
                                                              dist_span_cm[i,2,0])
                ang_mom_implied_y = angular_momentum_implied(gradient_y,
                                                              dist_span_cm[i,0,0],
                                                              dist_span_cm[i,2,0])
                ang_mom_implied_z = angular_momentum_implied(gradient_z,
                                                              dist_span_cm[i,0,0],
                                                              dist_span_cm[i,1,0])
                
                ang_mom_actual_total, ang_mom_actual_xy, ang_mom_actual_xz, ang_mom_actual_yz = angular_momentum_actual(clump,mass_clump_g[i])
                
                # =============================================================
                #     STORAGE OF VALUES FOR LATER USE
                grad_x.append(gradient_x)
                grad_y.append(gradient_y)
                grad_z.append(gradient_z)
                am_implied_x.append(ang_mom_implied_x)
                am_implied_y.append(ang_mom_implied_y)
                am_implied_z.append(ang_mom_implied_z)
                am_actual_total.append(ang_mom_actual_total)
                am_actual_partial_xy.append(ang_mom_actual_xy)
                am_actual_partial_xz.append(ang_mom_actual_xz)
                am_actual_partial_yz.append(ang_mom_actual_yz)
                # =============================================================
                

            except YTErrorReshape:
                grad_x.append(np.nan)
                grad_y.append(np.nan)
                grad_z.append(np.nan)
                am_implied_x.append(np.nan)
                am_implied_y.append(np.nan)
                am_implied_z.append(np.nan)
                am_actual_total.append(np.nan)
                am_actual_partial_xy.append(np.nan)
                am_actual_partial_xz.append(np.nan)
                am_actual_partial_yz.append(np.nan)
                err_string.append('Clump Number: '+
                                  str(i)+
                                  ' , could not reshape coordinates to 256x256 array')
                break
            
            except YTErrorValue:
                grad_x.append(np.nan)
                grad_y.append(np.nan)
                grad_z.append(np.nan)
                am_implied_x.append(np.nan)
                am_implied_y.append(np.nan)
                am_implied_z.append(np.nan)
                am_actual_total.append(np.nan)
                am_actual_partial_xy.append(np.nan)
                am_actual_partial_xz.append(np.nan)
                am_actual_partial_yz.append(np.nan)
                err_string.append('Clump Number: '+
                                  str(i)+
                                  ' , has v_positions empty')
                
                break
            
            #Break out of loop here if no errors are triggered
            # makes sure that the while loop doesn't get stuck at the end
            break

    
    print('Analysis Section Completed')
    # Turning all the data into numpy arrays to convert to Column objects
    clump_number = np.arange(1,len(bregion)+1)
    am_actual_total = np.array(am_actual_total)
    am_actual_partial_xy = np.array(am_actual_partial_xy)
    am_actual_partial_xz = np.array(am_actual_partial_xz)
    am_actual_partial_yz = np.array(am_actual_partial_yz)
    am_implied_x = np.array(am_implied_x)
    am_implied_y = np.array(am_implied_y)
    am_implied_z = np.array(am_implied_z)
    grad_x = np.array(grad_x)
    grad_y = np.array(grad_y)
    grad_z = np.array(grad_z)
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    #Creation of Columns for Astropy FITS for all quantities computed above
    
    q_clump_number = fits.Column(name='Clump Number',
                                 format='D',
                                 unit='dimensionless_unscaled',
                                 array=clump_number)
    
    #Actual Data for Angular Momentum
    q_am_actual_total = fits.Column(name='Actual Total Angular Momentum',
                                    format='D',
                                    unit=str(ang_mom_unit),
                                    array=am_actual_total)
    q_am_actual_partial_xy = fits.Column(name='Actual Partial Angular Momentum (x+y components)',
                                         format='D',
                                         unit=str(ang_mom_unit),
                                         array=am_actual_partial_xy)
    q_am_actual_partial_xz = fits.Column(name='Actual Partial Angular Momentum (x+z components)',
                                         format='D',
                                         unit=str(ang_mom_unit),
                                         array=am_actual_partial_xz)
    q_am_actual_partial_yz = fits.Column(name='Actual Partial Angular Momentum (y+z components)',
                                         format='D',
                                         unit=str(ang_mom_unit),
                                         array=am_actual_partial_yz)
    
    # Implied Data for Angular Momentum
    q_am_implied_x = fits.Column(name='Implied Total Angular momentum (x LOS)',
                                 format='D',
                                 unit=str(ang_mom_unit),
                                 array=am_implied_x)
    q_am_implied_y = fits.Column(name='Implied Total Angular momentum (y LOS)',
                                 format='D',
                                 unit=str(ang_mom_unit),
                                 array=am_implied_y)
    q_am_implied_z = fits.Column(name='Implied Total Angular momentum (z LOS)',
                                 format='D',
                                 unit=str(ang_mom_unit),
                                 array=am_implied_z)
    
    # Gradients for Implied Angular Momentum Plane Fitting
    q_grad_x = fits.Column(name='Plane Fitting Gradient (x LOS)',
                           format='D',
                           unit=str(per_sec_unit),
                           array=grad_x)
    q_grad_y = fits.Column(name='Plane Fitting Gradient (y LOS)',
                           format='D',
                           unit=str(per_sec_unit),
                           array=grad_y)
    q_grad_z = fits.Column(name='Plane Fitting Gradient (z LOS)',
                           format='D',
                           unit=str(per_sec_unit),
                           array=grad_z)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    
    #Creating ColDefs Object
    coldefs = fits.ColDefs([q_clump_number,
                            q_com_x,
                            q_com_y,
                            q_com_z,
                            q_mass_clump,
                            q_dist_span_x,
                            q_dist_span_y,
                            q_dist_span_z,
                            q_am_actual_total,
                            q_am_actual_partial_xy,
                            q_am_actual_partial_xz,
                            q_am_actual_partial_yz,
                            q_am_implied_x,
                            q_am_implied_y,
                            q_am_implied_z,
                            q_grad_x,
                            q_grad_y,
                            q_grad_z])
    
    # Creating HDU Object from ColDefs
    hdu = fits.BinTableHDU.from_columns(coldefs)
    print('ColDefs Object Created')
    hdu.header['MINCLMP'] = clump_sizing
    hdu.header['STEP'] = step
    hdu.header['BETA'] = beta
    hdu.header['LENGTH'] = l
    hdu.header['CMIN'] = cmin
    
    err_string_array = np.array(err_string)
    
    #For Loop for Adding in all the Error Statements for clumps (if any)
    for i in range(0,len(err_string_array)):
        hdu.header.add_comment(err_string_array[i])    
    
    
    #INSERT STRING CONNECTED TO DATAFILE INPUT FOR SCRIPT
    hdu.writeto(save_dir_specific+"data_"+fid_str+"_"+time_stamp+".fits",
                overwrite=True)
    print('FITS FILE SAVED')
    return()
