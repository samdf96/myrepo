#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 31 12:35:29 2018

@author: sfielder


Takes in global parameters, and writes them to a yaml file
"""

import yaml
import io

config = dict(l=10,cmin=5e-21,step=100,beta=1,clump_sizing=30)

with io.open(config_dir+'config_1.yaml', 'w') as outfile:
    yaml.dump(config, outfile, default_flow_style=True)

