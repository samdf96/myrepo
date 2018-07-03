#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 25 10:57:36 2018

@author: sfielder
"""


from astropy.io import fits
import glob
import matplotlib.pyplot as plt

## Import Definitions from definitions.py here:
from definitions import data_grabber
from definitions import j_comp_plotter
from definitions import j_comp_plotter_colormap


#Main Tree Start - LOCAL TESTING PARAMETERS HERE WILL USE FOR INPUT

## THIS WILL APPEAR IN THE INITIALIZE.PY FILE WHERE TESTING IS COMPLETE
data_dir = '/Users/sfielder/Documents/Astro_Data/'
flist = glob.glob(data_dir + '**/*.fits', recursive=True)
flist.sort()
current_file = flist[1]

# Here will start the actual definition process:
    
#def Plotter(current_file):


hdu = fits.open(current_file)
    
#Grabbing BinHDUTable object here, should always be second object
hdu_table = hdu[1]
data = hdu_table.data #Grabs data stored in the table -> FITS REC   
data_dict = data_grabber(data)

## Sort the data here:

header = hdu_table.header
comment_string = str(header['COMMENT'])
comments = comment_string.split(".")
hdu.close()

#%%

#Making the Output Directory for current file:
save_dir_list = current_file.split("/")[:-1]
save_dir = '/'.join(save_dir_list)

# Seperating out Data to be used for plotting
imp_x_los = data_dict['imp_x_los']
imp_y_los = data_dict['imp_y_los']
imp_z_los = data_dict['imp_z_los']
act_tot = data_dict['actual_tot']
act_par_xy = data_dict['actual_par_xy']
act_par_xz = data_dict['actual_par_xz']
act_par_yz = data_dict['actual_par_yz']
mass = data_dict['mass']

## Want Equal Axes?
equal_axis=False
# =============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    X LOS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
axis_str = 'X'
filename = 'j_comp_X_los.pdf'
savefig_dir = save_dir + '/' + filename

#Calling Function Here for Plotting
fig = j_comp_plotter(imp_x_los,
                     act_tot,
                     act_par_yz,
                     axis_str,
                     equal_axis)
plt.tight_layout()
fig.savefig(savefig_dir, bbox_inches='tight')
plt.close(fig)

#Colormap Variant Here
filename = 'j_comp_X_los_colormapped.pdf'
savefig_dir = save_dir + '/' + filename

#Calling Function Here for Plotting
fig = j_comp_plotter_colormap(imp_x_los,
                     act_tot,
                     act_par_yz,
                     mass,
                     axis_str,
                     equal_axis)
plt.tight_layout()
fig.canvas.draw()
fig.savefig(savefig_dir, bbox_inches='tight')
plt.close(fig)
# =============================================================================

# =============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Y LOS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
axis_str = 'Y'
filename = 'j_comp_Y_los.pdf'
savefig_dir = save_dir + '/' + filename

#Calling Function Here for Plotting
fig = j_comp_plotter(imp_y_los,
                     act_tot,
                     act_par_xz,
                     axis_str,
                     equal_axis)
plt.tight_layout()
fig.savefig(savefig_dir, bbox_inches='tight')
plt.close(fig)

#Colormap Variant Here
filename = 'j_comp_Y_los_colormapped.pdf'
savefig_dir = save_dir + '/' + filename

#Calling Function Here for Plotting
fig = j_comp_plotter_colormap(imp_y_los,
                     act_tot,
                     act_par_xz,
                     mass,
                     axis_str,
                     equal_axis)
plt.tight_layout()
fig.canvas.draw()
fig.savefig(savefig_dir, bbox_inches='tight')
plt.close(fig)
# =============================================================================

# =============================================================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    Z LOS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
axis_str = 'Z'
filename = 'j_comp_Z_los.pdf'
savefig_dir = save_dir + '/' + filename

#Calling Function Here for Plotting
fig = j_comp_plotter(imp_z_los,
                     act_tot,
                     act_par_xy,
                     axis_str,
                     equal_axis)
plt.tight_layout()
fig.savefig(savefig_dir, bbox_inches='tight')
plt.close(fig)

#Colormap Variant Here
filename = 'j_comp_Z_los_colormapped.pdf'
savefig_dir = save_dir + '/' + filename

#Calling Function Here for Plotting
fig = j_comp_plotter_colormap(imp_z_los,
                     act_tot,
                     act_par_xy,
                     mass,
                     axis_str,
                     equal_axis)
plt.tight_layout()
fig.canvas.draw()
fig.savefig(savefig_dir, bbox_inches='tight')
plt.close(fig)
# =============================================================================


#    return()

