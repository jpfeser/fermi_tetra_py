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
all_triangles = []
step = 2
for band in range(nb):
    rhomb_index = -1
    # for each starting tetrahedron index
    band_triangles = []
    for i0 in range(ng0[0]): #kx direction
        for i1 in range(ng0[1]): #ky direction
            for i2 in range(ng0[2]): #kz direction
                rhomb_index += 1
                # get the corners of the next rhombehedron
                rcd = tetra_toolbox.get_corners(i0, i1, i2, bvect, ng0, eig0[band,:,:,:])
                # then break it into six tetrahedra
                tcd = tetra_toolbox.break2tetra(rcd)
                for each_tetra in tcd.values():
                    # then interpolate the Fermi surface (triangles) inside each tetrahedron
                    triangle_set = tetra_toolbox.get_fermi_triangles(each_tetra,ef = 0.0)
                    for subset in triangle_set:
                        # then add that/those Fermi surf. triangles to the set of previous triangles
                        band_triangles.append(subset)
    #now plot the Fermi surface for that band
    tetra_toolbox.plot_triangle_set(band_triangles,bvect)
    all_triangles.append(band_triangles)
                # then I calculate whatever I want (lets sum E*df/dT(vg,z dot dA)) for each tetrahedra
                
