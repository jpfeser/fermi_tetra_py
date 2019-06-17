#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:13:00 2019

@author: jpfeser
"""
import numpy as np
import load_frmsf
import tetra_toolbox

FILENAME = "vfermi.frmsf"
# load data from frmsf file
(ng0, nb, avect, bvect, eig0, mat0, lshift) = load_frmsf.load_frmsrf(FILENAME)

#Now do some calculation for each band/tetrahedron
#for each band
for band in range(nb):
    rhomb_index = -1
    # for each starting tetrahedron index
    for i0 in range(ng0[0]): #kx direction
        for i1 in range(ng0[1]): #ky direction
            for i2 in range(ng0[2]): #kz direction
                rhomb_index += 1
                rcd = tetra_toolbox.get_corners(i0, i1, i2, bvect, ng0, eig0[band,:,:,:])
                tcd = tetra_toolbox.break2tetra(rcd)
#                for each_tetra in tcd:
#                    triangle_set = tetra_toolbox.get_fermi_triangles()
                    
                # then I calculate whatever I want (lets sum E*df/dT(vg,z dot dA)) for each tetrahedra
                # for i in range(6):  sum += eval_func(k_tetra[3,4,i])
                
