#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 13:13:00 2019

@author: jpfeser
"""
import numpy as np
import load_frmsf

filename="vfermi.frmsf"
# load data from frmsf file
(ng0,nb,avect,bvect,eig0,mat0,lshift)=load_frmsf.load_frmsrf(filename)

#Now do some calculation for each band/tetrahedron
#for each band
for band in range(nb):
    # for each starting tetrahedron index
    for i0 in range(ng0[0]): #kx direction 
        for i1 in range(ng0[1]): #ky direction 
            for i2 in range(ng0[2]): #kz direction 
#                for ii0 in range(2):
#                    for ii1 in range(2):
#                        for ii2 in range(2):
#                            kvect[]
                # we break it down to a cell, which itself can be broken into 6 tetrahedra
                #if points 0,1,2,3,4,5,6,7 define cell, then 
                #tetra1 = {0}
                
                # at this point I need the kpoints and eig/mat values associate with all eight corners.
                # (k_c[3,8],eigs_c[3,8],mat_c[3,8]) = get_corners())
                
                # then I reassign those eight corners to the six tetrahedra corners
                # k_tetr[3,4,6] <-> k_c[3,8]
                
                # then I calculate whatever I want (lets sum E*df/dT(vg,z dot dA)) for each tetrahedra
                # for i in range(6):  sum += eval_func(k_tetra[3,4,i])
                
                
            
    
#print("reassembling kmesh")
#for i in range(8):
#    for j in range(3):
#        kvec1[i,j]=kvec1[i,j]+shiftk[j]/(2*ng0[j])

    
    