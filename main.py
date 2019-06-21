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

# Material 1
FILENAME = "Al.frmsf"
mat1 = 'Al'
alat =  7.628*au2m; 

# Material 2
FILENAME2 = "Fe.frmsf"
mat2 = 'Fe'
alat2 =  5.4167*au2m; 
#Cu = 6.80
#Al = 7.628216862
#V = 5.671348546
#Ag = 7.7217
#Au = 7.681
#Ni = 6.612
#Fe = 5.4167

plots = 'on'
plots2 = 'on'
nint = np.array([1.,1.,1.])
tau = 1.0
T = 300
kbT = kb*T
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
(ng02, nb2, avect2, bvect2, eig02, mat02, lshift2, vf2) = load_frmsf.load_frmsrf(FILENAME2,alat2)

#Now extract the Fermi surface for each band
all_triangles = []
integral = 0.0
Atot = 0.0
vfA = 0.0
vinvA = 0.0
for band in range(nb):
    # for each starting tetrahedron index
    band_triangles = []
    band_Xtri = []
    for i0 in range(ng0[0]): #kx direction
        for i1 in range(ng0[1]): #ky direction
            for i2 in range(ng0[2]): #kz direction
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
                        vftri = Xtri_set[i][0:3,0:3]
                        dA,vfc,n = tetra_toolbox.get_area(triangle_set[i],vftri)
                        integral += np.abs(np.dot(vfc,nint))*dA/np.dot(vfc,n)
                        vfA += dA*np.linalg.norm(vfc)
                        vinvA += dA/np.linalg.norm(vfc)
                        Atot += dA
    #If desired plot the Fermi surface for each band (slow AF)
    if plots=='on':
        tetra_toolbox.plot_triangle_set(band_triangles,bvect)
    #Store the triangles for each band
    all_triangles.append(band_triangles)
                
#Thermal interface conductance
radlimit = prefactor*integral*tau

print('avg. Fermi velocity of material 1 is %.3e m/s' % (vfA/Atot))

Z1= 4*prefactor*integral
print('Z1 is %.3e' % Z1)

print('The calculated radiation limit for the thermal interface conductance  of %s is' % mat1)
print("G = % .3e W/m^2-K" % radlimit)



# Now calculate the DMM transmission coefficient
#Cu = 6.80
#Al = 7.628216862
#V = 5.671348546


#Now extract the Fermi surface for each band
all_triangles2 = []
integral2 = 0.0
Atot2 = 0.0
vfA2 = 0.0
vinvA2 = 0.0
for band in range(nb2):
    # for each starting tetrahedron index
    band_triangles2 = []
    band_Xtri2 = []
    for i0 in range(ng02[0]): #kx direction
        for i1 in range(ng02[1]): #ky direction
            for i2 in range(ng02[2]): #kz direction
                # get the corners of each rhombehedron
                rcd2 = tetra_toolbox.get_corners(i0, i1, i2, bvect2, ng02, eig02[band,:,:,:],vf2[band,:,:,:,0],vf2[band,:,:,:,1],vf2[band,:,:,:,2])
                # then break it into six tetrahedra
                tcd2 = tetra_toolbox.break2tetra(rcd2)
                for each_tetra in tcd2.values():
                    # then interpolate the Fermi surface (triangles) inside each tetrahedron
                    triangle_set2,Xtri_set2 = tetra_toolbox.get_fermi_triangles(each_tetra,ef = 0.0)
                    for i in range(len(triangle_set2)):
                        # then add that/those Fermi surf. triangles to the set of previous triangles
                        band_triangles2.append(triangle_set2[i])
                        band_Xtri2.append(Xtri_set2[i])
                        
                        # and do any calculations/interpolations/integration you need to do
                        vftri2 = Xtri_set2[i][0:3,0:3]
                        dA2,vfc2,n2 = tetra_toolbox.get_area(triangle_set2[i],vftri2)
                        integral2 += np.abs(np.dot(vfc2,nint))*dA2/np.dot(vfc2,n2)
                        vfA2 += dA2*np.linalg.norm(vfc2)
                        vinvA2 += dA2/np.linalg.norm(vfc2)
                        Atot2 += dA2
    #If desired plot the Fermi surface for each band (slow AF)
    if plots=='on':
        tetra_toolbox.plot_triangle_set(band_triangles2,bvect2)
    #Store the triangles for each band
    all_triangles2.append(band_triangles2)
    
tau = integral2/(integral + integral2)

result = prefactor*integral*tau

Z2= 4*prefactor*integral2
print('Z2 is %.3e' % Z2)

print('avg. Fermi velocity of material 2 is %.3e m/s' % (vfA2/Atot2))

print('the DMM transmission coefficient is\n tau = %f' % tau )

print('The calculated G for the %s-%s interface is' % (mat1,mat2))
print("G = % .3e W/m^2-K" % result)    
    
    
    