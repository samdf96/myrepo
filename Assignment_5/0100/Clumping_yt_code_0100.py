#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 19:52:29 2018

@author: sfielder
"""

import numpy as np

import yt
from yt.analysis_modules.level_sets.api import *
from yt.utilities.physical_constants import kboltz, mh
from yt.units import Kelvin

#Loading Data Set into Script
ds = yt.load("~/bigdata/Fiducial00/data.0100.3d.hdf5")

#All Data Object for Clumping Variable
ad = ds.all_data()
'''
for i in sorted(ds.derived_field_list):
    print(i)
'''
#Comment In to add in a Thermal Energy Density field for clumping
'''
def _therm(field, data):
   return(1.5 * (kboltz) * 10 * Kelvin / (2.32 * mh))
yt.add_field(("gas","thermal_energy"), function=_therm, units="erg/g", force_override=True)
'''


#Used for naming scheme of dbox_x areas
namespace = globals()
#Splitting Arrays for Octant Beginning and Ending Values for L length space
l = 10 #Data Length of One edge of original ds data set

x1 = np.array((-l/2,0,-l/2,0,-l/2,0,-l/2,0))
x2 = np.array((0,l/2,0,l/2,0,l/2,0,l/2))
y1 = np.array((-l/2,-l/2,0,0,-l/2,-l/2,0,0))
y2 = np.array((0,0,l/2,l/2,0,0,l/2,l/2))
z1 = np.array((-l/2,-l/2,-l/2,-l/2,0,0,0,0))
z2 = np.array((0,0,0,0,l/2,l/2,l/2,l/2))


#Creating 4 Loop for 8 Octants of data space to reduce computation overload
dbox_array = [] #Creating Empty List

for i in range(0, 1):
    dbox_array.append(ds.r[(x1[i],'pc'):(x2[i],'pc'), (y1[i],'pc'):(y2[i],'pc'), (z1[i],'pc'):(z2[i],'pc')])

    # the field to be used for contouring
    field = ("gas", "density")

    master_clump = Clump(dbox_array[i], ("gas", "density")) #Makes the first big clump
    clump_sizing = 50  #Any Clump size smaller than this value get eliminated
    master_clump.add_validator("min_cells", clump_sizing)
    '''
    master_clump.add_validator("gravitationally_bound",
                           use_particles=False,
                           use_thermal_energy=False)
    '''
    
    master_clump.add_info_item("center_of_mass") #Adds Center of Mass info for Clumps
    #Setting Limits for Clump Finding 
    c_min = 5e-21
    c_max = dbox_array[i]["gas", "density"].max()
    step = 1000
    find_clumps(master_clump, 5e-21, c_max, step) #Finds Clumps

    prj = yt.ProjectionPlot(ds,"z",field, data_source=dbox_array[i])    #Projcetion Plot
    prj.set_zlim('density', 1e-3, 3e-2) #Setting limits for all pictures to be the same
    prj.annotate_clumps(get_lowest_clumps(master_clump))    #Draws Contours on top of projection
    prj.save('Clump_0100_ClumpSizing_' + str(clump_sizing) + '_Octant_' + str(i) + '.png')  #Saves Fig

    #Grabs Lowest Tree Values (Leaves)
    lc = get_lowest_clumps(master_clump)
 
    center_of_mass = np.zeros((len(lc),3))
    box_region = np.zeros((len(lc),3,2))
    max_x_distance = np.zeros((len(lc),1))
    for j in range(0,len(lc)):
        com = np.array(lc[j].quantities.center_of_mass())   #Writes current com into array
        center_of_mass[j] = com     #Writes to center_of_mass array
        #print(lc[j].quantities.angular_momentum_vector())
        bregions = np.array(lc[j].quantities.extrema([('gas','x'),('gas','y'),('gas','z')]))
        box_region[j] = bregions # Creates Boundaries for clump box
        
        
        
        
        #print(lc[j].quantities.velocity_z())
    '''
    #Saves the Clumping Data Set
    #fn = master_clump.save_as_dataset('~/bigdata/Fiducial00/Clump_Data_Fielder/'+'%s_clump_%d' % (str(ds), i), fields=["density"])
    #cds = yt.load(fn)
    #leaf_clumps_reloaded = cds.leaves
    #print (cds.leaves)
    #print (cds.tree["clump", "center_of_mass"])

# We can traverse the clump hierarchy to get a list of all of the 'leaf' clumps
leaf_clumps = get_lowest_clumps(master_clump)



# Reload the clump dataset.
cds = yt.load(fn)

# Clump annotation can also be done with the reloaded clump dataset.

# Remove the original clump annotation
prj.annotate_clear()

# Get the leaves and add the callback.
leaf_clumps_reloaded = cds.leaves
prj.annotate_clumps(leaf_clumps_reloaded)
prj.save('clumps_reloaded_yt_code_0060.png')

# Query fields for clumps in the tree.
print (cds.tree["clump", "center_of_mass"])
print (cds.tree.children[0]["grid", "density"])
print (cds.tree.children[1]["all", "particle_mass"])

# Get all of the leaf clumps.
print (cds.leaves)
print (cds.leaves[0]["clump", "cell_mass"])
'''
