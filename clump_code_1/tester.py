#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 16 10:12:57 2018

@author: sfielder
"""
import os

import yt
import numpy as np
import logging
import astropy.units as u
from astropy.io import fits

# Import Functions from definitions.py file
from definitions import OctantSplit
from definitions import MasterClumpMaker
from definitions import ClumpFinder

from definitions import PlaneFit
from definitions import ArrayFlattener
from definitions import Gradient
from definitions import AngularMomentumImplied
from definitions import KineticEnergy
from definitions import GravitationalEnergy

#Local Variables
l = 10
cmin = 5.0e-21
beta = 1
clump_sizing = 30
step = 100


filename = '/Users/sfielder/Documents/Astro_Data/Fiducial00/data.0060.3d.hdf5'
save_dir_fits = '/Users/sfielder/Documents/Astro_Data/Output/'

logger = logging.getLogger("initialize.analyzer.Analyzer")

logger.info('Simulation File Current Working On: %s', filename)

#Loading File into Dataset
ds = yt.load(filename)
logger.debug('File Loaded')
master_dist_data = int(ds.domain_dimensions[0])
logger.debug('Master Distance Data set as: %s', master_dist_data)

#Might be Deprecated
err_string = []

main_string = filename.split("/")

#Make this modular in the future.
sim_string_true = [(('Fiducial' in k) or ('Design' in k)) for k in main_string]

#Grabbing Indice in main_string that corresponds to where Input Simulation string is.
sim_string_id = [i for i, x in enumerate(sim_string_true) if x]

sim_str = main_string[sim_string_id[0]] #sim_string_id should only have one entry!

out_string = main_string[-1].split(".") # Splits the data file name by periods
time_stamp = out_string[1]

#First Layer Here
save_dir = save_dir_fits + sim_str
logger.debug("General Save Directory string set as: %s", save_dir)
if os.path.isdir(save_dir) == True:
    logger.debug("Warning!!! Directory: %s , is detected as a valid directory. Files will be overwritten.",
                 save_dir)
else:
    logger.debug("Directory not detected. Creating Directory.")
    os.mkdir(save_dir)

    #Second Layer Here
save_dir_specific = save_dir + '/' + time_stamp + "/"
logger.debug("Specific Save Directory string set as: %s", save_dir_specific)
if os.path.isdir(save_dir_specific) == True:
    logger.debug("Warning!!! Directory: %s , is detected as a valid directory. FITS Files will be overwritten.",
                 save_dir_specific)
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
    logger.info("Working on Octant: %s", i+1)
    logger.debug("Invoking MasterClumpMaker function.")
    master_clump_main = MasterClumpMaker(octant[i])
    cmax = octant[i]["gas", "density"].max()
    logger.debug("cmax has been set to: %s", cmax)
    logger.debug("Invoking ClumpFinder function.")
    lc = ClumpFinder(master_clump_main,clump_sizing,cmin,cmax,step)

    #If lc returns empty the next for-loop will not be triggered. Write out message to err_string
    if len(lc) == 0:
        err_string.append('ClumpFinder has failed for Octant: ' + str(i+1))
        logger.info("Error String for FITS file has been written to. ClumpFinder failed.")
    
    for j in range(0,len(lc)):
        logger.debug("Appending Leaf Clumps to clumps list.")
        clumps.append(lc[j])

logger.info("Clump Finding Section Completed.")

logger.info("Data Filling Section Started.")

#Creating Arrays for later export (see bottom of script for this step)
logger.debug("Initialization data arrays.")
clump_number = np.arange(1,len(clumps)+1)
com = np.zeros((len(clumps),1,3))
com_x = np.zeros((len(clumps),1))
com_y = np.zeros((len(clumps),1))
com_z = np.zeros((len(clumps),1))
mass = np.zeros((len(clumps),1))
actual_angular_momentum = np.zeros((len(clumps),3))
actual_angular_momentum_x = np.zeros((len(clumps),1))
actual_angular_momentum_y = np.zeros((len(clumps),1))
actual_angular_momentum_z = np.zeros((len(clumps),1))
actual_angular_momentum_total = np.zeros((len(clumps),1))
actual_angular_momentum_par_xy = np.zeros((len(clumps),1))
actual_angular_momentum_par_xz = np.zeros((len(clumps),1))
actual_angular_momentum_par_yz = np.zeros((len(clumps),1))
implied_angular_momentum_x_los = np.zeros((len(clumps),1))
implied_angular_momentum_y_los = np.zeros((len(clumps),1))
implied_angular_momentum_z_los = np.zeros((len(clumps),1))
x_length = np.zeros((len(clumps),1))
y_length = np.zeros((len(clumps),1))
z_length = np.zeros((len(clumps),1))
volume_pix = np.zeros((len(clumps),1))
gradient_x_los = np.zeros((len(clumps),1))
gradient_y_los = np.zeros((len(clumps),1))
gradient_z_los = np.zeros((len(clumps),1))
bulk_velocity = np.zeros((len(clumps),3))
kinetic_energy = np.zeros((len(clumps),1))
gravitational_energy = np.zeros((len(clumps),1))
boundedness = np.zeros((len(clumps),1))
logger.info("Initialization of Data Arrays completed.")

#This value is used to correct for pixel resolution on length computation.
pixel = ((10 * u.pc).to(u.cm) / 256).value

#%%

for i in range(0,1):
        while True: #This is for catching IndexError if YT finds clump with no data
            try:
                logger.info("Currently Working on Clump Number %s out of %s",
                            str(i+1),
                            str(len(clumps)))
                logger.info("Computing Center of Mass Values.")
                com = clumps[i].quantities.center_of_mass()
                if str(com[0].value) == 'nan':
                    """
                    This will check to see if there is any data stored in the clump
                    this will raise an error and divert the code to the appropriate section
                    for dealing with this clump object.
                    Will return to the top of the for loop for the next clump.
                    """
                    logger.info("Center of Mass x detected to be nan. Raising error.")
                    raise ValueError
                
                logger.info("Setting COM and Mass Values.")
                com_x[i] = com[0].value
                com_y[i] = com[1].value
                com_z[i] = com[2].value
                logger.debug("Center_Of_Mass x Quantity found to be %s",
                             str(com_x[i]))
                logger.debug("Center_Of_Mass y Quantity found to be %s",
                             str(com_y[i]))
                logger.debug("Center_Of_Mass z Quantity found to be %s",
                             str(com_z[i]))
                '''
                Using the following technique to compute total mass.
                Originally used clumps[i].quantities.total_mass()
                This way takes more memory and time, the following is equivalent in terms
                of output, for actual data runs.
                '''
                mass[i] = clumps[i]['cell_mass'].sum()
                logger.debug("Mass found to be: %s", str(mass[i]))
                
                logger.info("Setting Coords and Length values.")
                x_coords = clumps[i]['x']
                logger.debug("x_coords found to be: %s", str(x_coords))
                y_coords = clumps[i]['y']
                logger.debug("y_coords found to be: %s", str(y_coords))
                z_coords = clumps[i]['z']
                logger.debug("z_coords found to be: %s", str(z_coords))
                
                x_length[i] = (x_coords.max()-x_coords.min()).value + pixel
                y_length[i] = (y_coords.max()-y_coords.min()).value + pixel
                z_length[i] = (z_coords.max()-z_coords.min()).value + pixel
                logger.debug("x_length found to be: %s", str(x_length[i]))
                logger.debug("y_length found to be: %s", str(y_length[i]))
                logger.debug("z_length found to be: %s", str(z_length[i]))
                
                #Using x_coords by convention to determine number of pixels in clump.
                volume_pix[i] = len(x_coords)
                logger.debug("volume_pix found to be %s", volume_pix[i])
        
                x_velocity = clumps[i]['velocity_x']
                y_velocity = clumps[i]['velocity_y']
                z_velocity = clumps[i]['velocity_z']
                logger.debug("x_velocity found to be: %s", str(x_velocity))
                logger.debug("y_velocity found to be: %s", str(y_velocity))
                logger.debug("z_velocity found to be: %s", str(z_velocity))
        
                #Add in for basic visualization of data.
            #    fig = plt.figure()
            #    ax = fig.add_subplot(111, projection='3d')
            #    ax.scatter(x_coords, y_coords, z_coords)
            
            
                x_coords_flat = ArrayFlattener(x_coords)
                y_coords_flat = ArrayFlattener(y_coords)
                z_coords_flat = ArrayFlattener(z_coords)
                x_velocity_flat = ArrayFlattener(x_velocity)
                y_velocity_flat = ArrayFlattener(y_velocity)
                z_velocity_flat = ArrayFlattener(z_velocity)
                logger.info("Coordinate and Velocity Arrays flattened.")
        
                actual_angular_momentum[i] = clumps[i].quantities.angular_momentum_vector()
                
                actual_angular_momentum_x[i] = actual_angular_momentum[i,0]
                actual_angular_momentum_y[i] = actual_angular_momentum[i,1]
                actual_angular_momentum_z[i] = actual_angular_momentum[i,2]
                
                #Using vectorial addition because of pos/neg quantities found.
                actual_angular_momentum_par_xy[i] = ((actual_angular_momentum_x[i]**2) +
                                              (actual_angular_momentum_y[i]**2))**(1/2)
                actual_angular_momentum_par_xz[i] = ((actual_angular_momentum_x[i]**2) +
                                              (actual_angular_momentum_z[i]**2))**(1/2)
                actual_angular_momentum_par_yz[i] = ((actual_angular_momentum_y[i]**2) +
                                              (actual_angular_momentum_z[i]**2))**(1/2)
                actual_angular_momentum_total[i] = ((actual_angular_momentum_x[i]**2) +
                                              (actual_angular_momentum_y[i]**2) + 
                                              (actual_angular_momentum_z[i]**2))**(1/2)
                
                logger.info("Actual Angular Momentum has been found.")
                
                # Implied Angular Momentum Calulations Here
                '''
                Plane Fitting is done using actual cm coordinate space, thus output
                coefficients for both gradient terms, are already in units of 1/s.
                Thus gradient_x/y/z_los are in units of 1/s.
                '''
                logger.info("Plane Fitting Commencing.")
                #For X LOS
                coeffs_x_los = PlaneFit(y_coords_flat,z_coords_flat,x_velocity_flat)
                gradient_x_los[i] = Gradient(coeffs_x_los)
                
                #For Y LOS
                coeffs_y_los = PlaneFit(z_coords_flat,x_coords_flat,y_velocity_flat)
                gradient_y_los[i] = Gradient(coeffs_y_los)
                
                #For Z LOS
                coeffs_z_los = PlaneFit(x_coords_flat, y_coords_flat, z_velocity_flat)
                gradient_z_los[i] = Gradient(coeffs_z_los)
                
                logger.info("Plane Fitting Completed.")
                
                implied_angular_momentum_x_los[i] = AngularMomentumImplied(gradient_x_los[i],
                                                                            y_length[i],
                                                                            z_length[i])
                implied_angular_momentum_y_los[i] = AngularMomentumImplied(gradient_y_los[i],
                                                                            x_length[i],
                                                                            z_length[i])
                implied_angular_momentum_z_los[i] = AngularMomentumImplied(gradient_z_los[i],
                                                                            x_length[i],
                                                                            y_length[i])
                
                ## Computation of Kinetic and Gravitational Energy
                bulk_velocity[i] = clumps[i].quantities.bulk_velocity(use_gas=True)
                kinetic_energy[i] = KineticEnergy(clumps[i],bulk_velocity[i])
                if volume_pix[i] < 10000:
                    gravitational_energy[i] = GravitationalEnergy(clumps[i],
                                        kinetic_energy[i])
                    if gravitational_energy[i] > kinetic_energy[i]:
                        logger.info("Clump is Gravitationall Bound.")
                        boundedness[i] = True
                    else:
                        boundedness[i] = False
                        logger.info("Clump is Not Gravitationally Bound.")
                else:
                    gravitational_energy[i] = 0
                    boundedness[i] = False
                    err_string.append("Pixel Volume Exceeds Threshold for Clump Number: " +
                                      str(i+1) + ". Gravitational Energy not being computed.")
                    
            except (ValueError, IndexError):
                err_string.append("Clump Number: " +
                                  str(i+1) +
                                  " has no data. Setting all values to nan.")
                logger.info("Setting all values to nan.")
                com_x[i] = np.nan
                com_y[i] = np.nan
                com_z[i] = np.nan
                mass[i] = np.nan
                actual_angular_momentum_x[i] = np.nan
                actual_angular_momentum_y[i] = np.nan
                actual_angular_momentum_z[i] = np.nan
                actual_angular_momentum_total[i] = np.nan
                actual_angular_momentum_par_xy[i] = np.nan
                actual_angular_momentum_par_xz[i] = np.nan
                actual_angular_momentum_par_yz[i] = np.nan
                implied_angular_momentum_x_los[i] = np.nan
                implied_angular_momentum_y_los[i] = np.nan
                implied_angular_momentum_z_los[i] = np.nan
                x_length[i] = np.nan
                y_length[i] = np.nan
                z_length[i] = np.nan
                volume_pix[i] = np.nan
                gradient_x_los[i] = np.nan
                gradient_y_los[i] = np.nan
                gradient_z_los[i] = np.nan
                bulk_velocity[i] = np.nan
                kinetic_energy[i] = np.nan
                gravitational_energy[i] = np.nan
                boundedness[i] = False # Needs to be a boolean for FITS File
                break # Out of Except and Restart the Loop at next iteration.
            break #Out of While Loop

