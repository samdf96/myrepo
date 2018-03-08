#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 03:45:19 2018

@author: sfielder
"""

#Advanced Plotting Script for Experimentation

import yt as yt
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset
ds = yt.load("~/bigdata/Fiducial00/data.0100.3d.hdf5")
'''
for i in sorted(ds.derived_field_list):
    print(i)
'''
# Make a density projection with annotated grids
p_annotate = yt.ProjectionPlot(ds, "z", "density")
p_annotate.annotate_grids()
p_annotate.save('z_density_annotated_grids.png')

# Density Projection with contour lines
p_contour = yt.SlicePlot(ds, "z", "density")
p_contour.annotate_contour("density")
p_contour.save('z_density_contour_density.png')

#Boxed Volume Render
sc = yt.create_scene(ds, ('gas', 'density'))
sc.annotate_domain(ds, color=[1, 1, 1, 0.1])
sc.save("%s_vr_domain.png" % ds, sigma_clip=4)
sc.annotate_grids(ds, alpha=0.01)
sc.save("%s_vr_grids.png" % ds, sigma_clip=4)
sc.annotate_axes(alpha=0.01)
sc.save("%s_vr_coords.png" % ds, sigma_clip=4)

#Transfer Function Plotting
sc = yt.create_scene(ds, lens_type='perspective')
source = sc[0]
source.tfh.set_bounds((1e-25, 1e-18))
source.tfh.set_log(True)
source.tfh.grey_opacity = True
source.tfh.plot('transfer_function.png', profile_field='density')
# save the image, flooring especially bright pixels for better contrast
sc.save('rendering_transfer.png', sigma_clip=6.0)