#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:13:00 2019

@author: jpfeser
"""
import numpy as np
import load_frmsf
import tetra_toolbox
from load_frmsf import hbar,kb,au2m
from numpy import pi

FILENAME = "Al.frmsf"
mat1 = 'Al'
alat =  7.628216862*au2m; 
#Al = 7.628216862
#V = 5.671348546
plots = 'on'
nint = np.array([1.,0.,0.])
tau = 1
T = 300
nspin = 1


# -----------INPUT ABOVE HERE--------------------------

if nspin is 1:
    prefactor = 2*kb**2*T/(48*pi*hbar)
elif nspin is 2:
    prefactor = kb**2*T/(48*pi*hbar)
else:
    print('error: nspin must be either 1 (non-magnetic) or 2 (magnetic)')

# load data from frmsf files
(ng0, nb, avect, bvect, eig0, mat0, lshift, vf) = load_frmsf.load_frmsrf(FILENAME,alat)

#Now extract the Fermi surface for each band
all_triangles = []
integral = 0.0
Atot = 0.0
for band in range(nb):
    rhomb_index = -1
    # for each starting tetrahedron index
    band_triangles = []
    band_Xtri = []
    for i0 in range(ng0[0]): #kx direction
        for i1 in range(ng0[1]): #ky direction
            for i2 in range(ng0[2]): #kz direction
                rhomb_index += 1
                # get the corners of each rhombehedron
                rcd = tetra_toolbox.get_corners(i0, i1, i2, bvect, ng0, eig0[band,:,:,:],vf[band,:,:,:,0],vf[band,:,:,:,1],vf[band,:,:,:,2])
                # then break it into six tetrahedra
                tcd = tetra_toolbox.break2tetra(rcd)
                for each_tetra in tcd.values():
                    # then interpolate the Fermi surface (triangles) inside each tetrahedron
                    triangle_set,Xtri_set = tetra_toolbox.get_fermi_triangles(each_tetra,ef = 0.0)
                    for i in range(len(triangle_set)):
                        # then add that/those Fermi surf. triangles to the set of previous triangles
                        band_triangles.append(triangle_set[i])
                        band_Xtri.append(Xtri_set[i])
                        
                        # and do any calculations/interpolations/integration you need to do
                        dA,vfc,n = tetra_toolbox.get_area(triangle_set[i],Xtri_set[i])
                        integral += np.abs(np.dot(vfc,nint))*dA/np.dot(vfc,n)
                        Atot += dA
    #If desired plot the Fermi surface for each band (slow AF)
    if plots=='on':
        tetra_toolbox.plot_triangle_set(band_triangles,bvect)
    #Store the triangles for each band
    all_triangles.append(band_triangles)
                
#Thermal interface conductance
result = prefactor*integral*tau

print('The calculated radiation limit for the thermal interface conductance  of %s is' % mat1)
print("G = % .3e W/m^2-K" % result)

# Now calculate the DMM transmission coefficient
