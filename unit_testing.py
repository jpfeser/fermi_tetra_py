#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:47:28 2019

@author: jpfeser
"""

import numpy as np
import tetra_toolbox

# test whether get_fermi_triangles gives good interpolation of triangles
temp = np.array([(0,0,0,1),(0,0,1,2),(0,1,0,3),(1,0,0,4)])
output = tetra_toolbox.get_fermi_triangles(temp,ef=2.5)
print(output)