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
from definitions import OctantSplit
from definitions import MasterClumpMaker
from definitions import ClumpFinder
from definitions import CenterOfMass
from definitions import BoundingBox
from definitions import ClumpBox
from definitions import VelocityArray
from definitions import VelocityArrayReducer
from definitions import ArrayFlattener
#Add in if needed
#from definitions import MyPlane
#from definitions import PlaneFitVisualization
from definitions import PlaneFit
from definitions import Gradient
from definitions import AngularMomentumImplied
from definitions import AngularMomentumActual
from definitions import ProjCreator

#Section for Importing Exceptions classes
from exceptions import YTErrorValue, YTErrorReshape, YTRuntimeError
from exceptions import YTPassThrough
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
import logging

module_logger = logging.getLogger("initialize.clump_code_1")

def Analyzer(filename,l,cmin,step,beta,clump_sizing,save_dir_fits):
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
    save_dir_fits: string
        Directory imput to save all FITS files
        
    Returns:
    -------- 
    FITS File contatining all the data computed.
        - Has all clump data
        - COMMENTS in the header for any clumps that returned errors
    
    
    '''
    logger = logging.getLogger("initialize.clump_code_1.Analyzer")
    
    logger.info('Simluation File Currently Working On: ',filename)
    
    #Loads data into File
    ds = yt.load(filename)
    logger.debug("File Loaded.")
    master_dist_data = int(ds.domain_dimensions[0])
    logger.debug("Master Distance Data set as: ", master_dist_data)
    err_string = [] #Used for tracking errors in functions
    
    # =========================================================================
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
    # =========================================================================
    
    
    # =========================================================================
    #Creating Sub Directory for individual Data Simulation File Timestamp
    
    #First Layer Here
    save_dir = save_dir_fits + fid_str
    logger.debug("General Save Directory string set as: ", save_dir)
    if os.path.isdir(save_dir) == True:
        logger.debug("Warning!!! Directory: " +
              save_dir +
              "is detected as a valid directory." +
              "Files will be overwritten.")
    else:
        logger.debug("Directory not detected. Creating Directory.")
        os.mkdir(save_dir)
    
    #Second Layer Here
    save_dir_specific = save_dir + '/' + time_stamp + "/"
    logger.debug("Specific Save Directory string set as: ", save_dir_specific)
    if os.path.isdir(save_dir_specific) == True:
        logger.debug("Warning!!! Directory: " +
              save_dir_specific +
              "is detected as a valid directory." +
              "FITS Files will be overwritten.")
    else:
        logger.debug("Specific Save Directory not detector. Creating Directory.")
        os.mkdir(save_dir_specific)
        
    logger.info("Directories have been set.")
    # =========================================================================
    #Creates a Data Object containing all the Simulation Data
    logger.debug("Creating 'ad' data object.")
    ad = ds.all_data()
    
    #Splits the Data into Octants
    logger.debug("Invoking OctantSplit function.")
    octant = OctantSplit(ad,ds,l)

    #Grabs each of the octants and runs them through the Clump Finding algorithm
    clumps = [] #Defining Empty List for loop
    logger.info("Clump Finding Section Started.")
    for i in range(0,len(octant)):
        logger.debug("Working on Octant: ", i+1)
        logger.debug("Invoking MasterClumpMaker function.")
        master_clump_main = MasterClumpMaker(octant[i])
        cmax = octant[i]["gas", "density"].max()
        logger.debug("cmax has been set to: ", cmax)
        logger.debug("Invoking ClumpFinder function.")
        lc = ClumpFinder(master_clump_main,clump_sizing,cmin,cmax,step)
        for j in range(0,len(lc)):
            logger.debug("Appending Leaf Clumps to clumps list.")
            clumps.append(lc[j])
    logger.info("Clump Finding Section Completed.")
    
    # =========================================================================
    # Making empty arrays to store the data, variable to allow for different
    # size of clumps list. Actual array (1,3) and (3,2) will not change as
    # that is how they are written with no difference between clump objects.
    # =========================================================================
    
    logger.debug("Data Filling Section Started.")
    
    #Creating Arrays for later export (see bottom of script for this step)
    logger.debug("Initialization com and bregion arrays.")
    com = np.zeros((len(clumps),1,3))
    bregion = np.zeros((len(clumps),3,2))
    
    #Defining Unit for Angular Momentum
    logger.debug("Setting Unit for angular momentum.")
    ang_mom_unit = u.kg * (u.m**2) / u.s
    per_sec_unit = u.s**(-1)
    
    # Creates lists for center_of_mass of clump list -> stores in com
    # Creates list for bounding boxes of clump list -> stores in bregion
    logger.debug("Filling com and bregion with data.")
    for i in range(0,len(clumps)):
    
        com[i] = CenterOfMass(clumps[i])
        bregion[i] = BoundingBox(clumps[i])

    # =========================================================================
    # CHECK IN POINT: What we have so far here in the script
    #
    #     com: arrays that contain x,y,z center of mass values for clumps
    #     bregion: x,y,z min/max value to build the boxes for clumps
    # ========================================================================        
    
    logger.debug("Initializing the following lists: data_object,",
                  "mass_clump, dist_span_cm, com_x, com_y, com_z,",
                  "dist_span_cm_x, dist_span_cm_y, dist_span_cm_z.")
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
        
    logger.debug("Creating Data Objects around Clump Objects.")
    for i in range(0,len(bregion)):
        logger.debug('Creating Clump Number:',i)
        data_object_clump.append(ClumpBox(ds,bregion[i]))
        logger.debug("Data Object Created.")
        #This loop runs to detect whether RuntimeError is triggered
        logger.debug("Appending mass to initialized list.")
        while True:
            try:
                mass_clump_g.append(np.array(data_object_clump[i].quantities.total_mass())[0])
                logger.debug("Mass Appended to List.")
            #If triggered, nan is written to mass array and 
            #error string is written with description of what happened
            except RuntimeError:
                logger.exception("RunTimeError has been caught.")
                logger.debug("Appending mass term as nan value.")
                mass_clump_g.append(np.array(np.nan))
                logger.debug("Appending err_string with relevant message.")
                err_string.append('Clump Number: '+
                                  str(i+1)+
                                  ' , has no width in one axial direction.')
                break
            break
            
        logger.debug("Creating dist_x,dist_y,dist_z from data object.")  
        dist_x = bregion[i,0,1] - bregion[i,0,0]
        dist_y = bregion[i,1,1] - bregion[i,1,0]
        dist_z = bregion[i,2,1] - bregion[i,2,0]
        
        logger.debug("Creating dist_tot from individual components.")
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
    logger.debug("Creating Tuples for Center of Mass along different LOS(s).")
    com_plotting = np.concatenate((com_x,com_y,com_z), axis=1)
    
    x_los_com = np.concatenate((com_y,com_z), axis=1) #(x,y) axis match (y,z)
    y_los_com = np.concatenate((com_z,com_x), axis=1) #(x,y) axis match (z,x)
    z_los_com = np.concatenate((com_x,com_y), axis=1) #(x,y) axis match (x,y)
    # Creatin of projection plots with markers for clump center of masses
    logger.debug("Invoking ProjCreator function.")
    ProjCreator(ds,
                ad,
                com_plotting,
                x_los_com,
                y_los_com,
                z_los_com,
                save_dir_specific,
                fid_str,
                time_stamp)
    
    #Turn into Array for later export (see bottom of script for this step)
    mass_clump_g = np.array(mass_clump_g)
    
    
    
    #Creation of Columns for Astropy FITS for all quantities computed above
    logger.debug("Creating Column Object for FITS file.")
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
    logger.debug("Initializing lists for exporting data to FITS file.")
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


    logger.debug("Data Filling Section Ended.")
    logger.debug("Plane Fitting Section Started.")
    
    for i in range(0,len(bregion)):
        clump = data_object_clump[i]
        
        logger.debug('Working on Analysis for Clump Number:',i)

        # =====================================================================
        # Computing Integrated Velocity Arrays
        logger.debug("Integrated Velocity Array Section Started.")
        while True:
            try:
                logger.debug("Invoking VelocityArray function.")
                arr_z, vz, broken = VelocityArray(clump,
                                                   'velocity_z',
                                                   'z',
                                                   master_dist_data,
                                                   l)  
                if broken==1:
                    logger.warning("Broken Value from function VelocityArray",
                                    " has been detected to be: ", broken)
                    raise YTRuntimeError
                
                logger.debug("Invoking VelocityArray function.")
                arr_x, vx, broken = VelocityArray(clump,
                                                   'velocity_x',
                                                   'x',
                                                   master_dist_data,
                                                   l)
                if broken==1:
                    logger.warning("Broken Value from function VelocityArray",
                                    " has been detected to be: ", broken)
                    raise YTRuntimeError
                
                logger.debug("Invoking VelocityArray function.")
                arr_y, vy, broken = VelocityArray(clump,
                                                   'velocity_y',
                                                   'y',
                                                   master_dist_data,
                                                   l)
                if broken==1:
                    logger.warning("Broken Value from function VelocityArray",
                                    " has been detected to be: ", broken)
                    raise YTRuntimeError
            
            #Setting Default value for detector variable
                runtime_error_detector = False
                
            except YTRuntimeError:
                logger.exception("YTRunTimeError Detected.")
                #Used below for pass through of values
                logger.debug("Setting runtime_error_detector to be TRUE.")
                runtime_error_detector = True
                break #Break out of except statement
            break       # Break out of try statement for while loop
        logger.debug("Integrated Velocity Array Section Completed.")          
        # =====================================================================
        
        # =====================================================================
        # Computing Reduced Velocity Arrays
        logger.debug("Reduced Velocity Array Section Started.")
        
        # This while loop raises errors if broken values from output of 
        # definitions file is not zero.
        
        while True:
            try:
                #This checks error value from above, if triggered, passes through
                #the code and writes nan values to the arrays
                if runtime_error_detector == True:
                    logger.debug("run_time_error_detector value is TRUE. Pass through.")
                    raise YTPassThrough
                
                arr_x_red, vx_py, vx_pz, broken = VelocityArrayReducer(arr_x,
                                                                       vx,
                                                                       'x',
                                                                       master_dist_data)
                #Check if broken statement is triggered
                if broken == 1:
                    logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                    raise YTErrorValue
                if broken == 2:
                    logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                    raise YTErrorReshape
                logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                logger.debug("X LOS Array Reducer Passed Without Issue.")
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                arr_y_red, vy_px, vy_pz, broken = VelocityArrayReducer(arr_y,
                                                                       vy,
                                                                       'y',
                                                                       master_dist_data)
                #Check if broken statement is triggered
                if broken == 1:
                    logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                    raise YTErrorValue
                if broken == 2:
                    logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                    raise YTErrorReshape
                logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                logger.debug("Y LOS Array Reducer Passed Without Issue.")
                # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
                arr_z_red, vz_px, vz_py, broken = VelocityArrayReducer(arr_z,
                                                                       vz,
                                                                       'z',
                                                                       master_dist_data)
                #Check if broken statement is triggered
                if broken == 1:
                    logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                    raise YTErrorValue
                if broken == 2:
                    logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                    raise YTErrorReshape
                logger.debug("Broken Value from VelocityArrayReducer detected to be: ", broken)
                logger.debug("Z LOS Array Reducer Passed Without Issue.")
                
                # =============================================================
                # Flattening Arrays for Plane Fitting Process
                logger.debug("Flattening Arrays for Plane Fitting.")
                logger.debug("Invoking ArrayFlattener function.")
                arr_x_red_flat = ArrayFlattener(arr_x_red)
                arr_y_red_flat = ArrayFlattener(arr_y_red)
                arr_z_red_flat = ArrayFlattener(arr_z_red)
                
                vx_py_flat = ArrayFlattener(vx_py)
                vx_pz_flat = ArrayFlattener(vx_pz)
                vy_px_flat = ArrayFlattener(vy_px)
                vy_pz_flat = ArrayFlattener(vy_pz)
                vz_px_flat = ArrayFlattener(vz_px)
                vz_py_flat = ArrayFlattener(vz_py)
                # =============================================================
                
                #Plane Fitting Process
                logger.debug("Invoking PlaneFit function.")
                result_x = PlaneFit(vx_py_flat,vx_pz_flat,arr_x_red_flat)
                result_y = PlaneFit(vy_px_flat,vy_pz_flat,arr_y_red_flat)
                result_z = PlaneFit(vz_px_flat,vz_py_flat,arr_z_red_flat)
                
                #Computing Gradients from Results
                logger.debug("Computing Gradients from Results.")
                logger.debug("Invoking Gradient function.")
                gradient_x = Gradient(result_x)
                gradient_y = Gradient(result_y)
                gradient_z = Gradient(result_z)
            
                #Computing Specific Angular Momentum
                # BETA VALUE INPUT TO DEFINITION WILL DEFAULT TO ONE
                logger.debug("Computing Specific Angular Momentum")
                ang_mom_implied_x = AngularMomentumImplied(gradient_x,
                                                           dist_span_cm[i,1,0],
                                                           dist_span_cm[i,2,0])
                ang_mom_implied_y = AngularMomentumImplied(gradient_y,
                                                           dist_span_cm[i,0,0],
                                                           dist_span_cm[i,2,0])
                ang_mom_implied_z = AngularMomentumImplied(gradient_z,
                                                           dist_span_cm[i,0,0],
                                                           dist_span_cm[i,1,0])
                
                ang_mom_actual_total, ang_mom_actual_xy, ang_mom_actual_xz, ang_mom_actual_yz = AngularMomentumActual(clump,
                                                                                                                      mass_clump_g[i])
                
                # =============================================================
                #     STORAGE OF VALUES FOR LATER USE
                logger.debug("Appending found values for later use in storage.")
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
                

            except YTErrorReshape: #Reshape Error - writes to error string
                logger.exception("YTErrorReshape caught.")
                logger.debug("Setting all relevent values to nan.")
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
                logger.debug("err_string appended with relevant message.")
                err_string.append('Clump Number: '+
                                  str(i+1)+
                                  ' , could not reshape coordinates to 256x256 array.')
                break
            
            except YTErrorValue:
                logger.exception("YTErrorValue caught.")
                logger.debug("Setting all relevent values to nan.")
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
                logger.debug("err_string appended with relevant message.")
                err_string.append('Clump Number: '+
                                  str(i+1)+
                                  ' , has v_positions empty.')
                break
            
            except YTPassThrough:
                logger.exception("YTPassThrough caught.")
                logger.debug("Setting all relevent values to nan.")
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
                logger.debug("err_string already appended with relevant message.")
                break
            break #Break out of while loop

    
    logger.debug("Reduced Velocity Array Section Completed.")
    # Turning all the data into numpy arrays to convert to Column objects
    logger.debug("Turning all data into numpy arrays.")
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
    logger.debug("Creating Column obects for new quantities.")
    
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
    q_am_implied_x = fits.Column(name='Implied Total Angular Momentum (x LOS)',
                                 format='D',
                                 unit=str(ang_mom_unit),
                                 array=am_implied_x)
    q_am_implied_y = fits.Column(name='Implied Total Angular Momentum (y LOS)',
                                 format='D',
                                 unit=str(ang_mom_unit),
                                 array=am_implied_y)
    q_am_implied_z = fits.Column(name='Implied Total Angular Momentum (z LOS)',
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
    logger.debug("coldefs object being created.")
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
    logger.debug("Creating HDU Object with Coldefs.")
    hdu = fits.BinTableHDU.from_columns(coldefs)
    logger.debug('ColDefs Object Created')
    hdu.header['MINCLMP'] = clump_sizing
    hdu.header['STEP'] = step
    hdu.header['BETA'] = beta
    hdu.header['LENGTH'] = l
    hdu.header['CMIN'] = cmin
    
    logger.debug("Setting Array from err_string.")
    err_string_array = np.array(err_string)
    
    #For Loop for Adding in all the Error Statements for clumps (if any)
    logger.debug("Appending err_string statements to header as COMMENT(s).")
    for i in range(0,len(err_string_array)):
        hdu.header.add_comment(err_string_array[i])    
    
    
    #INSERT STRING CONNECTED TO DATAFILE INPUT FOR SCRIPT
    hdu.writeto(save_dir_specific+fid_str+"_"+time_stamp+".fits",
                overwrite=True)
    logger.debug('FITS FILE SAVED')
    logger.info("Analyzer has been run successfully.")
    return()
