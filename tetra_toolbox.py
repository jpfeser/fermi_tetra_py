#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Toolset for manipulating rhombedra and tetrahedra for the fermi surfer info
"""
import numpy as np
import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import matplotlib.pyplot as plt

def get_corners(i0, i1, i2, bvect, ng0, eig0, *vargin):
    '''
    inputs:
        i0:  index along first bvect component
        i1:  index along second bvect component
        i2:  index along third bvect component
        bvect:  recip lattice vect, arranged columnwise (i.e. b1= bvect[:,0])
        ng0:  list with number of k-mesh points in each direction [ng00 ng01 ng02]
        eig0:  numpy array with shape [nb,ng0[0],ng0[1],ng0[2]]
    output:  
        *vargin: additional quantites with same shape as eig0 to be broken into corners of rhombehdron
        rcd: a dictionary with corner indices of the rhombehedra as the keys
            the associated values are a tuple (float kc[3],float eigc) 
    '''
    ip0 = (i0 + 1) % ng0[0]
    ip1 = (i1 + 1) % ng0[1]
    ip2 = (i2 + 1) % ng0[2]
    

#    k0in = (i0/(ng0[0]-1))*bvect[:,0]
#    k0out = ((i0 + 1)/(ng0[0]-1))*bvect[:,0]
#    k1in = (i1/(ng0[1]-1))*bvect[:,1]
#    k1out = ((i1 + 1)/(ng0[1]-1))*bvect[:,1]
#    k2in = (i2/(ng0[2]-1))*bvect[:,2]
#    k2out = ((i2 + 1)/(ng0[2]-1))*bvect[:,2]
    
    k0in = (i0/(ng0[0]))*bvect[:,0]
    k0out = ((i0 + 1)/(ng0[0]))*bvect[:,0]
    k1in = (i1/(ng0[1]))*bvect[:,1]
    k1out = ((i1 + 1)/(ng0[1]))*bvect[:,1]
    k2in = (i2/(ng0[2]))*bvect[:,2]
    k2out = ((i2 + 1)/(ng0[2]))*bvect[:,2]
    
    rcd = {}
    rcd[0]= (k0in+k1in+k2in, eig0[i0,i1,i2])
    rcd[1]= (k0out+k1in+k2in, eig0[ip0,i1,i2])
    rcd[2]= (k0in+k1out+k2in, eig0[i0,ip1,i2])
    rcd[3]= (k0in+k1in+k2out, eig0[i0,i1,ip2])
    rcd[4]= (k0out+k1out+k2in, eig0[ip0,ip1,i2])
    rcd[5]= (k0in+k1out+k2out, eig0[i0,ip1,ip2])
    rcd[6]= (k0out+k1in+k2out, eig0[ip0,i1,ip2])
    rcd[7]= (k0out+k1out+k2out, eig0[ip0,ip1,ip2])
    
    for j in range(len(vargin)):
        rcd[0] += (vargin[j][i0,i1,i2],)
        rcd[1] += (vargin[j][ip0,i1,i2],)
        rcd[2] += (vargin[j][i0,ip1,i2],)
        rcd[3] += (vargin[j][i0,i1,ip2],)
        rcd[4] += (vargin[j][ip0,ip1,i2],)
        rcd[5] += (vargin[j][i0,ip1,ip2],)
        rcd[6] += (vargin[j][ip0,i1,ip2],)
        rcd[7] += (vargin[j][ip0,ip1,ip2],)  
    
    return rcd

def break2tetra(rcd):
    '''
    take a rhombedra dictionary and break into into the corresponding tetrahedra.
    
    inputs:  rcd, rhombedra dictionary contains a set of 8 keys/values:  {corner number:tuple}
        where tuple is (kc0,kc1,kc2,eigvalue,...) where kc0,kc1,kc2 are the x,y,z 
        value of k at the corner, and eigvalue is the associated energy, 
        varg can be additional tuple entries and are treated similarly to eigvalue otherwise.
        
    output:  tcd is a dictionary with one key for each tetrahedron (6 total),
           the values are a list (length 4 + len(varg)) corresponding to the information at corner.
           At each corner (i.e. list value), a tuple is given (kc[0:2],energy,vargs).  
           The list is sorted such that the first tuple's energy is the 
           smallest (ascending).  For definitions of the colors see the ArXiv
           article (https://arxiv.org/pdf/1811.06177.pdf)
        
    e.g. 
    '''
    assert isinstance(rcd,dict), 'error:  input for break2tetra expected to be dictionary'
    assert len(rcd)==8, 'input dictionary rcd does not recieve exactly 8 corners in break2tetra'

    def sort_energy(tup): 
        '''pulls out energy from rcd tuple'''
        return tup[1]
    
    #note: values of rcd of form: (float kc[3],float eigc)

    # start creating dictionary for red tetrahedra (set of points {2,6,0,3})
    redlist = [rcd[2], rcd[6], rcd[0], rcd[3]]
    redlist.sort(key = sort_energy)
    
    # start creating dictionary for blue tetrahedra (set of points {2,6,0,1})
    bluelist = [rcd[2], rcd[6], rcd[0], rcd[1]]
    bluelist.sort(key = sort_energy)
    #print(redlist)
    
    # start creating dictionary for green tetrahedra (set of points {2,6,1,4})
    greenlist = [rcd[2], rcd[6], rcd[1], rcd[4]]
    greenlist.sort(key = sort_energy)
    #print(redlist)
    
    # start creating dictionary for magenta tetrahedra 
    #(set of points {2,6,4,7}
    maglist = [rcd[2], rcd[6], rcd[4], rcd[7]]
    maglist.sort(key = sort_energy)
    #print(redlist)
    
    # start creating dictionary for yellow tetrahedra 
    #(set of points {2,6,5,7}
    yellowlist = [rcd[2], rcd[6], rcd[5], rcd[7]]
    yellowlist.sort(key = sort_energy)
    
    # start creating dictionary for cyan tetrahedra 
    #(set of points {2,6,3,5}
    cyanlist = [rcd[2], rcd[6], rcd[3], rcd[5]]
    cyanlist.sort(key = sort_energy)
    
    # Q: now what should I return? 
    # A: a dictionary of property lists, one key for each tetrahedron
    #       the values are the tuples for each of four corners (kc[0:2],energy)
    #       sorted such that the first tuple energy is smallest.
    tcd = {}
    tcd['r'] = redlist
    tcd['b'] = bluelist
    tcd['g'] = greenlist
    tcd['m'] = maglist
    tcd['y'] = yellowlist
    tcd['c'] = cyanlist
    
    return tcd    

def get_fermi_triangles(tcd_list,ef=0.0):
    '''
    get_fermi_triangles(tcd_list,Ef) determine the location of the triangle 
    intersections on the tetrahedron.
    
    Inputs:
        tcd_list:  a list (length 4) corresponding to the information at corner.
       At each corner (i.e. list value), a tuple is given (kc[0:2],energy).  
       The list is sorted such that the first tuple's energy is the 
       smallest (ascending).  ex. [([0,0,0],1),([0,0,1],1),([0,1,0],1),([0,0,1],2)]
       
       ef = Fermi energy (default = 0.0)
       
    Outputs:
        triangles = list (length = # triangles) where each value in the 
                          list is a 3 x 3 np array (Xtri) with the corners of 
                          the triangle.  An empty list is returned if there 
                          are no triangles to return. 
                          rows are triangle #, 
                          columns are x,y,z values respectively
                          
                          ex. 
                          [0 0 1;
                          1 0 1;
                          0 1 0;]
                          
                          means 
                          point 0: [kx,ky,kz]=[0,0,1]
                          point 1: [kx,ky,kz]=[1,0,1]
                          point 2: [kx,ky,kz]=[0,1,0]
                          
    '''
    e = np.zeros(4)
    k = np.zeros((4,3))
    ktri = np.zeros((3,3))
    a = np.zeros((4,4))
    F = np.zeros((3,4))
    Nval = len(tcd_list[0])-2
    Xval = np.zeros((4,Nval)); #note assuming extra entries are floats for now

    
    #readout the values at the four corners (numbering should be presorted in 
    # ascending energy
    for i in range(4):
        k[i,:]=tcd_list[i][0] # first element of tuple is k vector
        e[i]=tcd_list[i][1] # second element of tuple is energy
        # if there are additional entries, read them out too.
        for j in range(Nval):
            Xval[i,j] = tcd_list[i][2+j]
        
    triangles = []
    vertices = []
    if ef>e[0] and ef<=e[1]:
        for i in range(4):
            for j in range(4):
                de = e[i]-e[j]
                if de == 0.0:
                    a[i,j]=0.5
                else:
                    a[i,j]=(ef-e[j])/de
                
        F =np.array([[a[0,1], a[1,0],     0.0, 0.0], \
                     [a[0,2],    0.0,  a[2,0], 0.0], \
                     [a[0,3],    0.0,   0.0,   a[3,0]]])
        
        ktri = np.matmul(F,k)
        triangles.append(ktri)
        
        # 3 x N (n=0 -> energy, others depend on tcd xtras)
        Xtri = np.matmul(F,Xval)
        vertices.append(Xtri)
        
        return triangles,vertices
    
    elif ef>e[1] and ef<e[2]:
        for i in range(4):
            for j in range(4):
                de = e[i]-e[j]
                if de == 0.0:
                    a[i,j]=0.5
                else:
                    #the paper has a typo: says a[i,j]=de/(ef-e[j])
                    # but that is clearly wrong!!!
                    a[i,j]=(ef-e[j])/de
                
        F =np.array([[a[0,2],    0.0, a[2,0],    0.0], \
                     [a[0,3],    0.0,    0.0, a[3,0]], \
                     [   0.0, a[1,3],    0.0, a[3,1]]])
        
        ktri = np.matmul(F,k)
        triangles.append(ktri)
        # 3 x N (n=0 -> energy, others depend on tcd xtras)
        Xtri = np.matmul(F,Xval)
        vertices.append(Xtri)
        
        F =np.array([[a[0,2],    0.0, a[2,0],    0.0], \
                     [   0.0, a[1,2], a[2,1],    0.0], \
                     [   0.0, a[1,3],    0.0, a[3,1]]])
        
        ktri = np.matmul(F,k)
        triangles.append(ktri)
                # 3 x N (n=0 -> energy, others depend on tcd xtras)
        Xtri = np.matmul(F,Xval)
        vertices.append(Xtri)
        
        return triangles,vertices
    
    elif ef>=e[2] and ef<e[3]:
        for i in range(4):
            for j in range(4):
                de = e[i]-e[j]
                if de == 0.0:
                    a[i,j]=0.5
                else:
                    a[i,j]=(ef-e[j])/de

# this is what the paper says, but I think it's a typo               
#        F =np.array([[a[0,3],    0.0,     0.0, a[3,0]], \
#                     [a[0,2], a[1,3],     0.0, a[3,1]], \
#                     [a[0,1],    0.0,  a[2,3], a[3,2]]])

# I think it should be this: proceed with caution                              
        F =np.array([[a[0,3],    0.0,     0.0, a[3,0]], \
                     [   0.0, a[1,3],     0.0, a[3,1]], \
                     [   0.0,    0.0,  a[2,3], a[3,2]]])
        
        ktri = np.matmul(F,k)
        triangles.append(ktri)
        # 3 x N (n=0 -> energy, others depend on tcd xtras)
        Xtri = np.matmul(F,Xval)
        vertices.append(Xtri)
        return triangles,vertices

    else:  #(ef < e[0]) or (ef >= e[3]):
        return triangles,vertices

# END get_fermi_triangles(tcd_list,ef) function definition
        
def plot_triangle_set(triangle_set,*argv):
    '''
    Input: triangle_set = a list of 3x3 numpy arrays.  Each numpy array describes the 3 kpoints of a triangle. 
                                      rows are triangle #, 
                                      
                          columns are x,y,z values respectively
                          
                          ex. 
                          [0 0 1;
                          1 0 1;
                          0 1 0;]
                          
                          means 
                          point 0: [kx,ky,kz]=[0,0,1]
                          point 1: [kx,ky,kz]=[1,0,1]
                          point 2: [kx,ky,kz]=[0,1,0]
                          
    Output: 3D plot of Fermi surface (no return)
    '''

    
    fig = plt.figure()
    ax = a3.Axes3D(fig)
    if len(argv) is 1:
        bvect = argv[0]
        for i in range(3):
            start_pt = bvect[:,i]
            end_pt = start_pt + bvect[:,(i + 1) % 3]
            ax.plot([start_pt[0],end_pt[0]],[start_pt[1],end_pt[1]],[start_pt[2],end_pt[2]])
            
            end_pt = start_pt + bvect[:,(i + 2) % 3]
            ax.plot([start_pt[0],end_pt[0]],[start_pt[1],end_pt[1]],[start_pt[2],end_pt[2]])
            
            end_pt = 0*start_pt
            ax.plot([start_pt[0],end_pt[0]],[start_pt[1],end_pt[1]],[start_pt[2],end_pt[2]])

            start_pt = bvect[:,0] + bvect[:,1] + bvect[:,2]
            end_pt = start_pt - bvect[:,i]
            ax.plot([start_pt[0],end_pt[0]],[start_pt[1],end_pt[1]],[start_pt[2],end_pt[2]])
        
    for vtx_array in triangle_set:
        tri = a3.art3d.Poly3DCollection([vtx_array])
        tri.set_edgecolor('none')
        tri.set_edgecolor('yellow')
        ax.add_collection3d(tri)
    plt.show()
    return

def get_area(ktri,*varg):
    '''
    input:  ktri = 3 x 3 numpy array (i = corner #, j=x,y,z) specifying a 3D triangle
            vf (optional): 3 x 3 numpy array with fermi velocity (i = corner #, j=x,y,z)
    output: dA = area of triangle
            vfc = Fermi velocity at centroid (3,1) if vf is entered.
            n = unit normal based on interpolated vf if vf is entered.
    '''
    a = ktri[1,:]-ktri[0,:];
    b = ktri[2,:]-ktri[0,:];
    dA = np.linalg.norm(np.cross(a,b))/2;
    nalt = np.cross(a,b)/(2*dA) # this does not seem to be same as n *gulp*
    
    if len(varg) is 0:
        return dA
    elif len(varg) is 1:
        vf = varg[0]
        #compute the avg. value of vf at centroid
        #position of centroid
        kc = (ktri[0,:]+ktri[1,:]+ktri[2,:])/3
        #linear intepolation
        mid12 = (ktri[1,:]+ktri[2,:])/2
        d0c = np.linalg.norm(kc-ktri[0,:])
        dcmid12 = np.linalg.norm(kc-mid12)
        dtot = d0c + dcmid12
        vfmid12 = (vf[1,:]+vf[2,:])/2
        vfc = vfmid12*dcmid12/dtot + vf[0,:]*d0c/dtot
        n = vfc/np.linalg.norm(vfc)
        return dA,vfc,n
    else:
        print('error in get_area - too many arguments')

        