#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Toolset for manipulating rhombedra and tetrahedra for the fermi surfer info
"""
import numpy as np

def get_corners(i0,i1,i2,bvect,ng0,eig0):
    '''
    inputs:
        i0:  index along first bvect component
        i1:  index along second bvect component
        i2:  index along third bvect component
        bvect:  recip lattice vect, arranged columnwise (i.e. b1= bvect[:,0])
        ng0:  list with number of k-mesh points in each direction [ng00 ng01 ng02]
        eig0:  numpy array with shape [nb,ng0[0],ng0[1],ng0[2]]
    output:  
        rcd: a dictionary with corner indices of the rhombehedra as the keys
            the associated values are a tuple (float kc[3],float eigc) 
    '''
    ip0 = (i0 + 1) % ng0[0]
    ip1 = (i1 + 1) % ng0[1]
    ip2 = (i2 + 1) % ng0[2]

    k0in = (i0/(ng0[0]-1))*bvect[:,0]
    k0out = (ip0/(ng0[0]-1))*bvect[:,0]
    k1in = (i1/(ng0[1]-1))*bvect[:,1]
    k1out = (ip1/(ng0[1]-1))*bvect[:,1]
    k2in = (i2/(ng0[2]-1))*bvect[:,2]
    k2out = (ip2/(ng0[2]-1))*bvect[:,2]
    
    rcd={}
    rcd[0]=( k0in+k1in+k2in, eig0[i0,i1,i2])
    rcd[1]=( k0out+k1in+k2in, eig0[ip0,i1,i2])
    rcd[2]=( k0in+k1out+k2in, eig0[i0,ip1,i2])
    rcd[3]=( k0in+k1in+k2out, eig0[i0,i1,ip2])
    rcd[4]=( k0out+k1out+k2in, eig0[ip0,ip1,i2])
    rcd[5]=( k0in+k1out+k2out, eig0[i0,ip1,ip2])
    rcd[6]=( k0out+k1in+k2out, eig0[ip0,i1,ip2])
    rcd[7]=( k0out+k1out+k2out, eig0[ip0,ip1,ip2])
    
    return rcd

def break2tetra(rcd):
    '''
    take a rhombedra dictionary and break into into the corresponding tetrahedra.
    
    inputs:  rcd, rhombedra dictionary contains a set of 8 keys/values:  {corner number:tuple}
        where tuple is (kc0,kc1,kc2,eigvalue) where kc0,kc1,kc2 are the x,y,z 
        value of k at the corner, and eigvalue is the associated energy typically.
        
    output:  tcd a dictionary (6 keys, one for each tetrahedra) where each 
        contained dictionary has keys corresponding to the 4 corners as numbered
        in the fermisurfer ArXiv article.
        
    e.g. 
    '''
    assert isinstance(rcd,dict),'error:  input for break2tetra expected to be dictionary'
    assert len(rcd)==8,'input dictionary rcd does not have 8 corners in break2tetra'

    def sortE(tup):
        return tup[1]

    # start creating dictionary for red tetrahedra (set of points {2,6,0,3})
    redlist = [rcd[2],rcd[6],rcd[0],rcd[3]]
    redlist.sort(key=sortE)
    
    # start creating dictionary for blue tetrahedra (set of points {2,6,0,1})
    bluelist = [rcd[2],rcd[6],rcd[0],rcd[1]]
    bluelist.sort(key=sortE)
    #print(redlist)
    
    # start creating dictionary for green tetrahedra (set of points {2,6,1,4})
    greenlist = [rcd[2],rcd[6],rcd[1],rcd[4]]
    greenlist.sort(key=sortE)
    #print(redlist)
    
    # start creating dictionary for magenta tetrahedra 
    #(set of points {2,6,4,7}
    maglist = [rcd[2],rcd[6],rcd[4],rcd[7]]
    maglist.sort(key=sortE)
    #print(redlist)
    
    # start creating dictionary for yellow tetrahedra 
    #(set of points {2,6,5,7}
    yellowlist = [rcd[2],rcd[6],rcd[5],rcd[7]]
    yellowlist.sort(key=sortE)
    
    # start creating dictionary for cyan tetrahedra 
    #(set of points {2,6,3,5}
    cyanlist = [rcd[2],rcd[6],rcd[3],rcd[5]]
    cyanlist.sort(key=sortE)
    
    #now what should I return?
    
    
    
        