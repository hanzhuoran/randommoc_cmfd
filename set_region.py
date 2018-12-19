from ray import Ray
import numpy as np
from geom import *
import math
from util import *
from moc import *
import matplotlib.pyplot as plt
import copy
from material import *
import time
import os

def set_geometry(flag):
    #######################################
    #############set up geometry###########
    #######################################



    ############Set up all planes##########
    xbd = np.empty(NX+1,dtype=Plane)
    for i in range(NX+1):
        if i == 0:
            xbd[i] = Plane(1,0,0,xbd_list[i],xmbd,i)
        elif i == NX:
            xbd[i] = Plane(1,0,0,xbd_list[i],xpbd,i)
        else: 
            xbd[i] = Plane(1,0,0,xbd_list[i],"transmission",i)

    ybd = np.empty(NY+1,dtype=Plane)
    for j in range(NY+1):
        if j == 0:
            ybd[j] = Plane(0,1,0,ybd_list[j],ymbd,NX+j+1)
        elif j == NY:
            ybd[j] = Plane(0,1,0,ybd_list[j],ypbd,NX+j+1)
        else: 
            ybd[j] = Plane(0,1,0,ybd_list[j],"transmission",NX+j+1)

    zbd = np.empty(2,dtype=Plane)
    zbd[0] = Plane(0,0,1,sidelengthZ/2,"reflection",NX+1+NY+1)
    zbd[1] = Plane(0,0,1,-sidelengthZ/2,"reflection",NX+1+NY+2)

    pln_list = np.append(xbd,ybd)
    pln_list = np.append(pln_list,zbd)
    boxlist = determine_boxlist(xbd,ybd)
    cyl_list = np.empty((NY,NX,n_rings),dtype=Surface)
    
    pos_x_array = np.linspace(-(NX-1)/2,(NX-1)/2,NX)
    pos_y_array = np.linspace((NY-1)/2,-(NY-1)/2,NY)
    index = 0

    # print("n rings")
    # print(n_rings)
    if n_rings > 0:
        for i, x in enumerate(pos_x_array):
            for j,y in enumerate(pos_y_array):
                x0 = x*pitch
                y0 = y*pitch
                for k, r in enumerate(rlist):
                    cyl_list[j,i,k] = Z_Cylinder(x0,y0,r,"transmission",index)
                    index +=1
        all_surf = np.append(pln_list,cyl_list)
    elif n_rings == 0:
        all_surf = pln_list

    # print("all_sruf_size")
    # print(all_surf.shape)
    ##########set up regions############

    Regions = np.empty((NY,NX,n_rings+1),dtype=Region)

    id = 0
    for j in range(NY):
        for i in range(NX): 
            for k in range(n_rings+1):
                if n_rings > 0:
                    if k == 0:
                        surf_list = [cyl_list[j,i,k]]
                        surf_dir = [-1]
                        Regions[j,i,k] = Region(surf_list,surf_dir,id)
                    elif k == n_rings:
                        surf_list = np.append(cyl_list[j,i,k-1],boxlist[j,i,:])
                        surf_dir = [1, 1, -1 ,-1 ,1]
                        Regions[j,i,k] = Region(surf_list,surf_dir,id)
                    else:
                        surf_list = cyl_list[j,i,(k-1):k+1]
                        surf_dir = [1, -1] 
                        Regions[j,i,k] = Region(surf_list,surf_dir,id)
                elif n_rings == 0:
                    surf_list = boxlist[j,i,:]
                    surf_dir = [1, -1 ,-1 ,1]
                    Regions[j,i,k] = Region(surf_list,surf_dir,id)
                id +=1
    if flag == 1:
        process_material(XSdirectory, N_GROUPS, Regions)
    elif flag == 2:
        process_homo_material(homo_XSdirectory, N_GROUPS, Regions)
    
    # array = np.linspace(0,np.pi/2,NX+1)
    # mid_x = 1/2*(array[0:NX]+array[1:NX+1])
    # # print(mid_x)
    # for j in range(NY):
    #     for i in range(NX): 
    #         for k in range(n_rings+1):
    #             val = np.sin(mid_x[i])*np.ones(N_GROUPS)
    #             # print("val")
    #             # print(val)
    #             Regions[j,i,k].set_phi(val)
    return Regions, all_surf