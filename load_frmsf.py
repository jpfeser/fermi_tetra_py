#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Load important quantities from fermisurfer file
input:  filename := string with frmsf file location

returns:  (ng0,nb,avect,bvect,eig0,mat0,lshift) 
where   ng0:= k mesh size, [N1,N2,N3], size: [3,]
        nb:= number of band in frmsf file
        bvect:= 3x3 numpy array with each column a recip lattice vector
        avect:= 3x3 numpy array with each column a lattice vector
        eig0:= energy at each kpoint (size: [nb,N1,N2,N3] )
        lshift:= integer describing kmesh (0=Monkhorst Pack,1=contains Gamma, 2=?)
"""
import numpy as np

def load_frmsrf(filename):
    fo = open(filename,"r")
    
    #Now grid density
    print("loading grid density:\n")
    ng0line = fo.readline().split()
    ng0=np.asarray(ng0line, dtype=np.int8)
    print(ng0)
    
    #Now read shift type
    print("loading shift type (0=Monkhorst Pack,1=contains Gamma,2=?):\n")
    lshift = int(fo.readline().strip())
    #shiftk=[]
    print(lshift)
    #if (lshift == 0):
    #    print("    k point grid is the Monkhorst-Pack grid. \n");
    #    for i in np.arange(3):
    #        shiftk.append((ng0[i] + 1) % 2)
    #elif (lshift == 1):
    #    print("    k point grid starts from Gamma. \n")
    #    for i in np.arange(3):
    #        shiftk.append(0)
    #elif (lshift == 2):
    #    print("    k point grid starts from Gamma + a half grid. \n");
    #    for i in np.arange(3):
    #        shiftk.append(1);
            
    #Now read number of bands
    print("reading number of bands, nb")
    nb = int(fo.readline().strip())
    print(nb)
    
    print("Read reciprocal lattice vectors, b (in columns)")
    bvect = np.array(np.zeros([3,3]))
    bvect[:,0] = np.array(fo.readline().strip().split(),dtype='float64')
    bvect[:,1] = np.array(fo.readline().strip().split(),dtype='float64')
    bvect[:,2] = np.array(fo.readline().strip().split(),dtype='float64')
    print(bvect)
    
    print("Calculating direct lattice vectors (in columns)")
    avect = np.linalg.inv(bvect.transpose()) #from wikipedia "Reciprocal lattice"
    # verified correct using kittel test cases.
    # note that a factor 1/(2*pi) is missing  (may be important)
    print(avect)
    
    eig0=np.zeros([nb,ng0[0],ng0[1],ng0[2]])
    for ib in np.arange(nb):
        for i0 in np.arange(ng0[0]):
            if (lshift != 0):
                ii0 = i0;
            else:
                ii0 = (i0 + (ng0[0] + 1) / 2) % ng0[0];
            for i1 in np.arange(ng0[1]):
                if (lshift != 0):
                    ii1 = i1;
                else:
                    ii1 = (i1 + (ng0[1] + 1) / 2) % ng0[1];
                for i2 in np.arange(ng0[2]):
                    if (lshift != 0):
                        ii2 = i2;
                    else:
                        ii2 = (i2 + (ng0[2] + 1) / 2) % ng0[2];
                    eig0[ib,ii0,ii1,ii2]= float(fo.readline().strip());
    
    print("Energies have been loads into array with shape\n (nbands,nk1,nk2,nk3)=",eig0.shape)                
    
    mat0=np.zeros([nb,ng0[0],ng0[1],ng0[2]])
    for ib in np.arange(nb):
        for i0 in np.arange(ng0[0]):
            if (lshift != 0):
                ii0 = i0;
            else:
                ii0 = (i0 + (ng0[0] + 1) / 2) % ng0[0];
            for i1 in np.arange(ng0[1]):
                if (lshift != 0):
                    ii1 = i1;
                else:
                    ii1 = (i1 + (ng0[1] + 1) / 2) % ng0[1];
                for i2 in np.arange(ng0[2]):
                    if (lshift != 0):
                        ii2 = i2;
                    else:
                        ii2 = (i2 + (ng0[2] + 1) / 2) % ng0[2];
                    mat0[ib,ii0,ii1,ii2]= float(fo.readline().strip());
     
    print("mat0 has been loaded into array with shape\n (nbands,nk1,nk2,nk3)=",mat0.shape)                 
    print("closing ",filename)
    return (ng0,nb,avect,bvect,eig0,mat0,lshift)
