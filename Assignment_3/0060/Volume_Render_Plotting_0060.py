#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 19:55:20 2018

@author: sfielder
"""

import yt as yt
import matplotlib.pyplot as plt
import numpy as np

# Load the dataset
ds = yt.load("~/bigdata/Fiducial00/data.0060.3d.hdf5")

# Create a volume rendering
sc = yt.create_scene(ds, field=('gas', 'density'))

# Now increase the resolution
sc.camera.resolution = (1024, 1024)

# Set the camera focus to a position that is offset from the center of
# the domain
sc.camera.focus = ds.arr([0.3, 0.3, 0.3], 'unitary')

# Move the camera position to the other side of the dataset
sc.camera.position = ds.arr([0, 0, 0], 'unitary')

# save to disk with a custom filename and apply sigma clipping to eliminate
# very bright pixels, producing an image with better contrast.
sc.render()
sc.save('Volume_Render_Blue_0060.png', sigma_clip=4)