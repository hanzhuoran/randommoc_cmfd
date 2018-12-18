import numpy as np
import math
import matplotlib.pyplot as plt
import copy
import csv
import os

from ray import Ray
from geom import *
from util import *
from moc import *
from material import *
from CMFD_func import *
from set_region import *
################


###always start with this seed
seed(2)
print("We have", N_GROUPS, "energy groups")

CMFD = 0
plot = 0
if plot == 1:
    print("NRINGS:",n_rings)
    if n_rings > 0:
        ax=plt.gca()
        pos_x_array = np.linspace(-(NX-1)/2,(NX-1)/2,NX)
        pos_y_array = np.linspace((NY-1)/2,-(NY-1)/2,NY)


        for i, xx in enumerate(pos_x_array):
            for j,yy in enumerate(pos_y_array):
                x0 = xx*pitch
                y0 = yy*pitch
                for r in rlist:
                    circle = plt.Circle((x0, y0), r, color='r',fill=False)
                    ax.add_artist(circle)

        plt.axis('scaled')
        ax.set_xlim((xbd_list[0]-0.5, xbd_list[-1]+0.5))
        ax.set_ylim((ybd_list[-1]-0.5, ybd_list[0]+0.5))

    # plot x planes
    for i in range(NX+1):
         plt.plot([-NX*pitch/2+i*pitch,-NX*pitch/2+i*pitch],\
                    [ybd_list[0],ybd_list[-1]],'b-')

    # plot y planes
    for j in range(NY+1):
         plt.plot([xbd_list[0],xbd_list[-1]],
                    [NY*pitch/2-j*pitch,NY*pitch/2-j*pitch],'b-')

# plt.show()

#######################################
#############set up geometry###########
#######################################
latticebox = [sidelengthX/2, sidelengthY/2, sidelengthZ/2,
             -sidelengthX/2, -sidelengthY/2, -sidelengthZ/2]
# 1 is hetero
# 2 is homo
Regions, all_surf  = set_geometry(2)


########################################
#############ITERATION##################
########################################
k = 1
l = 0
tight = 0 
CMFD_iter = INF
converge = 0
conv2 = 1
conv1 = 1


# while not converge:
for l in range(1):
#     if l >= 100:
#         break
    surf_cnt = 0
    seed()
    l += 1
    print("#########################")
    print("#####iteration",l,"######")
    print("#########################")

    ##Calc Q and initialize for this loop
    FissionRate = 0
    for y in range(NY):
        for x in range(NX):
            for c in range(n_rings+1):
                Regions[y,x,c].calc_q(k)
                Regions[y,x,c].reset_length()
                Regions[y,x,c].phi_old = Regions[y,x,c].phi
                # newnf = np.multiply(Regions[y,x,c].Sig_nuf.transpose(), \
                                # Regions[y,x,c].phi)
                # FissionRate +=  newnf*pitch**2  
                Regions[y,x,c].reset_phi()
             
    k_pre = k


    current_filter_x = np.zeros((n_x_planes,NY,N_GROUPS))
    current_filter_y = np.zeros((n_y_planes,NX,N_GROUPS))

    #####Loop Rays#######
    Dtotal = 0
    ray_cnt = 0
    tot_hit_cnt = 0
    for i in range(np.size(all_surf)):
        all_surf[i].cnt = 0

    for ray_iter in range(N_RAYS):
        ray = Ray(latticebox)
        y, x, c = find_region(ray,Regions)
        init_psi = copy.deepcopy(Regions[y,x,c].phi_old)
        init_psi /= 4*np.pi
        init_psi = np.asarray(init_psi)
        init_psi = np.squeeze(init_psi)
        if N_GROUPS == 1:
            ray.set_psi([init_psi])
        else:
            ray.set_psi(init_psi)
        # print("init psi")
        # print(ray.psi)


    #################Start Sweeping Ray###############
        D = 0
        while (D < MAX_D):
            s = INF
            ID = -1

            ###FIND nearest distance
            for i in range(np.size(all_surf)):
                if type(all_surf[i]) is Plane:
                    d = dis2plane(all_surf[i],ray)
                elif type(all_surf[i]) is Z_Cylinder:
                    d = dis2cyl(all_surf[i],ray)
                if d < s and d > 0 :
                    ID = i
                    s = d

            moving_L = s

            y, x, c = find_region(ray,Regions)

            if (D < DZ) and (s+D) > DZ:
                s = s+D-DZ
                DeadZone(ray,Regions[y,x,c],(DZ-D))
                D = DZ  

            ###Accumulate D and check if it exceeds MAX_D
            D += s
            if D > MAX_D:
                s = s+MAX_D-D
                D = MAX_D

            #######Attenuation
            if D >= DZ:
                Regions[y,x,c].length += s
                Dtotal += s
                attenuate_segment(ray,Regions[y,x,c],s)
            else:      
                DeadZone(ray,Regions[y,x,c],s)

            if plot == 1:
	            previous_x = ray.x
	            previous_y = ray.y
         
            ray.to_next_surf(moving_L)

            #Check if relective 
            if all_surf[ID].type == "reflection":
                # Regular Reflective Boundaries
                ray.reflect(all_surf[ID])
                if  D > DZ and type(all_surf[ID]) is Plane:
                    tot_hit_cnt +=2
                    all_surf[ID].cnt +=2
            elif all_surf[ID].type == "transreflect" :
                if  D > DZ and type(all_surf[ID]) is Plane:
                    tot_hit_cnt +=2
                    all_surf[ID].cnt +=2
                    tally_current(all_surf[ID],ray,current_filter_x, \
                                        current_filter_y)
                ray.reflect(all_surf[ID]) 
                ray.set_psi(np.zeros(N_GROUPS))
            else:
                # Inner surfaces "transmission"
                if  D > DZ and type(all_surf[ID]) is Plane:
                    tot_hit_cnt +=1
                    all_surf[ID].cnt +=1
                    tally_current(all_surf[ID],ray,current_filter_x, \
                                        current_filter_y)

            ray.to_next_surf(1e-14)

            if plot == 1:
                plt.plot([previous_x, ray.x],[previous_y, ray.y], 'bo-')
            # print("===========================")

        # ray_cnt += ray.out

    for y in range(NY):
        for x in range(NX):
            for c in range(n_rings+1):
                print("pos")
                print(y,x,c)
                Regions[y,x,c].calc_volume(Dtotal)
                print(Dtotal)
                print(Regions[y,x,c].volume)
                # region.volume*Dtotal = region.length
                term1 = Regions[y,x,c].phi/(Regions[y,x,c].Sig_tot[:,0].transpose() \
                                    * Regions[y,x,c].volume*Dtotal)
                term2 = Regions[y,x,c].q*4*np.pi
                print("term1")
                print(term1)
                print("term2")
                print(term2)
                Regions[y,x,c].phi = term1+ term2

#########################################################
################# Condensation and CMFD ##################
##########################################################
    # 2 GROUP CMFD
    if CMFD == 1:
        print("USE CMFD")
        CMFD_Sig_tot = np.zeros((NY,NX,NG))
        CMFD_phi = np.zeros((NY,NX,NG))
        CMFD_Sig_abs = np.zeros((NY,NX,NG))
        CMFD_Sig_nuf = np.zeros((NY,NX,NG))
        CMFD_Chi = np.zeros((NY,NX,NG))
        CMFD_D = np.zeros((NY,NX,NG))
        CMFD_Sig_sct = np.zeros((NY,NX,NG,NG))
        condense_xs(Regions, CMFD_Sig_tot,CMFD_Sig_abs, CMFD_Sig_nuf, \
                    CMFD_Chi, CMFD_D, CMFD_Sig_sct,CMFD_phi)

#############FURTHER CONDENSATION##################
        # Condense more cells into a coarser cell
        if ratio > 1:    
            CCMFD_Sig_tot = np.zeros((NY,int(NX/ratio),NG))
            CCMFD_phi = np.zeros((NY,int(NX/ratio),NG))
            CCMFD_Sig_abs = np.zeros((NY,int(NX/ratio),NG))
            CCMFD_Sig_nuf = np.zeros((NY,int(NX/ratio),NG))
            CCMFD_Chi = CMFD_Chi[:,0::ratio,:] # Ad hoc, alwyas 1 and 0 afterall
            CCMFD_D = np.zeros((NY,int(NX/ratio),NG))
            CCMFD_Sig_sct = np.zeros((NY,int(NX/ratio),NG,NG))

            for j in range(NY):
                for i in range(int(NX/ratio)):
                    for ii in range(ratio):
                        CCMFD_phi[j,i,:] += CMFD_phi[j,i*ratio+ii,:] 
                        CCMFD_Sig_abs[j,i,:] += CMFD_Sig_abs[j,i*ratio+ii,:]*CMFD_phi[j,i*ratio+ii,:]
                        CCMFD_Sig_nuf[j,i,:] += CMFD_Sig_nuf[j,i*ratio+ii,:]*CMFD_phi[j,i*ratio+ii,:]
                        CCMFD_D[j,i,:] += CMFD_D[j,i*ratio+ii,:]*CMFD_phi[j,i*ratio+ii,:]
                        for g in range(NG):
                            for gp in range(NG):
                                CCMFD_Sig_sct[j,i,gp,g] += CMFD_Sig_sct[j,ratio*i+ii,gp,g]*CMFD_phi[j,ratio*i+ii,g]
                    CCMFD_Sig_abs[j,i,:] /=  CCMFD_phi[j,i,:]
                    CCMFD_Sig_nuf[j,i,:] /=  CCMFD_phi[j,i,:]
                    CCMFD_Sig_sct[j,i,:] /=  CCMFD_phi[j,i,:]
                    CCMFD_D[j,i,:] /=  CCMFD_phi[j,i,:]
                    CCMFD_phi[j,i,:] /= ratio


            CCMFD_Sig_R = calc_sig_r(CCMFD_Sig_abs,CCMFD_Sig_sct)

            JJx = current_filter_x[0::ratio,:,:]
            JJy = current_filter_y[:,0::ratio,:]+current_filter_y[:,1::ratio,:]

            Delta_y = pitch
            Delta_x = ratio*pitch
            CCMFD_J_x, CCMFD_J_y = condense_J(JJx/Dtotal,JJy/Dtotal,Delta_x,Delta_y)
            CCMFD_Dhat_x,CCMFD_Dhat_y = get_Dhat(CCMFD_D)
            CCMFD_Dtil_x,CCMFD_Dtil_y = get_Dtil(CCMFD_phi,CCMFD_Dhat_x, \
                                                CCMFD_Dhat_y,CCMFD_J_x,CCMFD_J_y)

            CM = Set_M(CCMFD_Dtil_x,CCMFD_Dtil_y,CCMFD_Dhat_x,CCMFD_Dhat_y, \
                        CCMFD_Sig_sct,CCMFD_Sig_R)
            CF =  Set_F(CCMFD_Sig_nuf,CCMFD_Chi)

            phi_CCM = CCMFD_phi.reshape(int(NX/ratio)*NY,NG)
            phi_CCM = phi_CCM.reshape((int(NX/ratio)*NY*NG,1),order='F')
            phi_CCM = phi_CCM/np.linalg.norm(phi_CCM)
            phi_CCM = np.mat(phi_CCM)
        else:
            # Coarse mesh is just the fuel pin
            CMFD_Sig_R = calc_sig_r(CMFD_Sig_abs,CMFD_Sig_sct)
            CMFD_J_x, CMFD_J_y = condense_J(current_filter_x/Dtotal,current_filter_y/Dtotal,pitch,pitch)
            CMFD_Dhat_x,CMFD_Dhat_y = get_Dhat(CMFD_D)
            CMFD_Dtil_x = np.zeros((n_x_planes,NY,N_GROUPS))
            CMFD_Dtil_y = np.zeros((n_y_planes,NX,N_GROUPS))    
            # CMFD_Dtil_x,CMFD_Dtil_y = get_Dtil(CMFD_phi,CMFD_Dhat_x, \
            #                                     CMFD_Dhat_y,CMFD_J_x,CMFD_J_y)

            CM = Set_M(CMFD_Dtil_x,CMFD_Dtil_y,CMFD_Dhat_x,CMFD_Dhat_y, \
                        CMFD_Sig_sct,CMFD_Sig_R)
            CF =  Set_F(CMFD_Sig_nuf,CMFD_Chi)

            phi_CCM = CMFD_phi.reshape(NX*NY,NG)
            phi_CCM = phi_CCM.reshape((NX*NY*NG,1),order='F')
            phi_CCM = phi_CCM/np.linalg.norm(phi_CCM)
            phi_CCM = np.mat(phi_CCM)

        k,phi_CCM_new = diffusion_solver(k,phi_CCM,CM,CF)
        print(phi_CCM_new)
        ##########################################
        ################ prologation################
        ##########################################
        factor = np.divide(phi_CCM_new,phi_CCM)
        if ratio > 1:
            full_factor = np.mat(np.zeros((NX*NY*NG,1)))
            for i in range(factor.shape[0]):
                full_factor[ratio*i:ratio*(i+1)] = factor[i]
            norm_pre = update_CMFD_phi(full_factor,Regions)
        else:
            norm_pre = update_CMFD_phi(factor,Regions)

        ###########################################
        ################ Update phi################
        for y in range(NY):
            for x in range(NX):
                for c in range(n_rings+1):
                    Regions[y,x,c].phi =  norm_pre[y,x,c]/ \
                                            np.linalg.norm(Regions[y,x,c].phi)
############################################################

    # Calc K-eff
  
    FR = 0
    if CMFD == 0:
        tot_FR = 0
        old_FR = 0
        tot_abs = 0
   
    for y in range(NY):
        for x in range(NX):
            for c in range(n_rings+1):
                print("y,x,c", y,x,c)
                print(Regions[y,x,c].phi)
                newnf = np.multiply(Regions[y,x,c].Sig_nuf.transpose(), \
                                Regions[y,x,c].phi)
                print(newnf)
                oldnf = np.multiply(Regions[y,x,c].Sig_nuf.transpose(), \
                                Regions[y,x,c].phi_old)
                print(oldnf)
                # newabs = np.multiply(Regions[y,x,c].Sig_abs.transpose(), \
                #                 Regions[y,x,c].phi)
                if CMFD == 0:
                    tot_FR += np.sum(newnf*pitch**2)
                    old_FR += np.sum(oldnf*pitch**2)
                    # tot_abs += np.sum(newabs**pitch**2) 
    
                print("previous fr")
                print( old_FR)
                print("new FR")  
                print(tot_FR)

    if CMFD == 0:
        k = tot_FR/old_FR*k_pre
    print("k eff:",k)

    phi_new_array = np.zeros(N_REGION)
    phi_old_array = np.zeros(N_REGION)
    index = 0
    for y in range(NY):
        for x in range(NX):
            for c in range(n_rings+1):
                phi_new_array[index] = Regions[y,x,c].phi
                phi_old_array[index] = Regions[y,x,c].phi_old
                index += 1
    print("phi_old_array")
    print(phi_old_array)

    print("phi_new_array")
    print(phi_new_array)

    # RMSD = np.sqrt(FR/(NX*NY*(n_rings+1)))
    conv1 = np.abs((k-k_pre)/k)
    conv2 = np.linalg.norm((phi_new_array-phi_old_array)/phi_new_array)

    print("K_ERR",conv1)
    print("flux",conv2)
    # print("thermal_err",conv3)
    # print("RMSD:", RMSD)
    # print("k_pre:", k_pre)
    # print("k_poweriteration:", k)


    if conv1< 1e-5 and conv2 < 1e-5:
        converge = 1
        # CMFD = 0
        # if l < CMFD_iter:
        #     CMFD_iter = l
        #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        #     print("####K converged",CMFD_iter,"#####")
        #     print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        # if RMSD < conv_crit2:
        #     converge = 1

    # Normalization
    sum = 0
    for y in range(NY):
        for x in range(NX):
            for c in range(n_rings+1):
                sum += np.sum(Regions[y,x,c].phi)
    norm = sum/(NY*NX*(n_rings+1)*N_GROUPS)
    print("norm")
    print(norm)


    for y in range(NY):
        for x in range(NX):
            for c in range(n_rings+1):
                Regions[y,x,c].phi = Regions[y,x,c].phi/norm

    #Go to Next iter
##############ITERARION OVER###########
print("Finished")
if plot == 1:
    plt.show()


# print(Regions[0,0,0].Sig_abs)
# print(Regions[0,0,0].Sig_sct)
# print(Regions[0,0,0].Sig_tot)
# print(Regions[0,0,0].phi)
# print(Regions[0,0,1].phi)
# plotphi1 = Regions[0,0,0].phi
# plotphi2 = Regions[0,0,1].phi
# group_grid = [0, 0.058,0.14,0.28,0.625,4,10,40,5530,821e3,20e6]
# plt.step(plotphi1, group_grid[1:-1], where='pre', label='fuel')
# plt.step(plotphi2, group_grid[1:-1], where='pre', label='mod')
# print(Regions[0,0,0].Sig_sct)
#############PLOT#####################
# if plot == 1:
#     if n_rings > 1:
#         ax=plt.gca()
#         pos_x_array = np.linspace(-(NX-1)/2,(NX-1)/2,NX)
#         pos_y_array = np.linspace((NY-1)/2,-(NY-1)/2,NY)


#         for i, xx in enumerate(pos_x_array):
#             for j,yy in enumerate(pos_y_array):
#                 x0 = xx*pitch
#                 y0 = yy*pitch
#                 for r in rlist:
#                     circle = plt.Circle((x0, y0), r, color='r',fill=False)
#                     ax.add_artist(circle)

#         plt.axis('scaled')
#         ax.set_xlim((xbd_list[0]-0.5, xbd_list[-1]+0.5))
#         ax.set_ylim((ybd_list[-1]-0.5, ybd_list[0]+0.5))

#     # plot x planes
#     for i in range(NX+1):
#          plt.plot([-NX*pitch/2+i*pitch,-NX*pitch/2+i*pitch],\
#                     [ybd_list[0],ybd_list[-1]],'b-')

#     # plot y planes
#     for j in range(NY+1):
#          plt.plot([xbd_list[0],xbd_list[-1]],
#                     [NY*pitch/2-j*pitch,NY*pitch/2-j*pitch],'b-')



####
# if plot == 1:
#     plt.figure(2)
#     axis = 1/2*(xbd_list[0:NX]+xbd_list[1:NX+1])
#     print("axis")
#     print(axis)
#     plotphi = np.squeeze(np.array(fast))
#     plotphi /= np.linalg.norm(plotphi)
#     print("plotphi")
#     print(plotphi)
#     plt.plot(axis,plotphi)
#     plt.show()

# if plot == 1:
#     if n_rings > 1:
#         ax=plt.gca()
#         pos_x_array = np.linspace(-(NX-1)/2,(NX-1)/2,NX)
#         pos_y_array = np.linspace((NY-1)/2,-(NY-1)/2,NY)


#         for i, xx in enumerate(pos_x_array):
#             for j,yy in enumerate(pos_y_array):
#                 x0 = xx*pitch
#                 y0 = yy*pitch
#                 for r in rlist:
#                     circle = plt.Circle((x0, y0), r, color='r',fill=False)
#                     ax.add_artist(circle)

#         plt.axis('scaled')
#         ax.set_xlim((xbd_list[0]-0.5, xbd_list[-1]+0.5))
#         ax.set_ylim((ybd_list[-1]-0.5, ybd_list[0]+0.5))

#     # plot x planes
#     for i in range(NX+1):
#          plt.plot([-NX*pitch/2+i*pitch,-NX*pitch/2+i*pitch],\
#                     [ybd_list[0],ybd_list[-1]],'b-')

#     # plot y planes
#     for j in range(NY+1):
#          plt.plot([xbd_list[0],xbd_list[-1]],
#                     [NY*pitch/2-j*pitch,NY*pitch/2-j*pitch],'b-')

##########################################
# array = np.linspace(0,np.pi,NX+1)
# mid_x = 1/2*(array[0:NX]+array[1:NX+1])
# ana = np.sin(mid_x)
# ana /= np.linalg.norm(ana)
# plt.plot(mid_x,ana,label="ana")
# plt.plot(axis,plotphi,label="phi")
# plt.legend()
# plt.show()
# for y in range(NY):
#     for x in range(NX):
#         for c in range(n_rings+1):
#             region = Regions[y,x,c]
#             print(region.phi)
# n_planes,3,N_GROUPS
# for i in range(n_planes):
#     for j in range(NX):
#         directory = './J/'+str(i)+'_'+str(j)
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         filename = 'J.txt'
#         np.savetxt(directory+'/'+filename, current_filter[i,j,:])


# for y in range(NY):
#     for x in range(NX):
#         for c in range(n_rings+1):
#             #####Calculate volume of each region
#             print("--------------")
#             print("xyc",y,x,c)
#             region = Regions[y,x,c]
#             directory = './rr_phi/'+str(x)+'_'+str(y)+'_'+str(c)
#             if not os.path.exists(directory):
#                 os.makedirs(directory)
#             filename='phi.txt'
#             np.savetxt(directory+'/'+filename, region.phi)
#             filename1='volume.txt'
#             print(region.volume)
#             np.savetxt(directory+'/'+filename1, [region.volume])


###################TEST###################
# # 2 GROUP CMFD
# CMFD_phi = np.zeros((NY,NX,NG))
# CMFD_Sig_tot = np.zeros((NY,NX,NG))
# CMFD_Sig_abs = np.zeros((NY,NX,NG)) 
# CMFD_Sig_nuf = np.zeros((NY,NX,NG))
# CMFD_Chi = np.zeros((NY,NX,NG))
# CMFD_D = np.zeros((NY,NX,NG))
# CMFD_Sig_sct = np.zeros((NY,NX,NG,NG))
# condense_xs(Regions, CMFD_Sig_tot ,CMFD_Sig_abs, CMFD_Sig_nuf, CMFD_Chi, CMFD_D, CMFD_Sig_sct,\
#              CMFD_phi)
# CMFD_Sig_R = calc_sig_r(CMFD_Sig_abs,CMFD_Sig_sct)


# print("CMFD_phi")
# print(CMFD_phi)

# print("CMFD_Sig_abs")
# print(CMFD_Sig_abs)

# print("CMFD_Sig_nuf")
# print(CMFD_Sig_nuf)

# print("CMFD_Chi")
# print(CMFD_Chi)

# print("CMFD_D")
# print(CMFD_D)

# print("CMFD_Sig_sct")
# print(CMFD_Sig_sct)


# print("CMFD_Sig_R")
# print(CMFD_Sig_R)

# print(current_filter_x)
# print("------------------")
# print(current_filter_y)
# for i in range(n_x_planes):
#     for j in range(NY):
#         directory = './J/'+str(i)+'_'+str(j)
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         filename = 'J.txt'
#         np.savetxt(directory+'/'+filename, current_filter_x[i,j,:])
# for i in range(n_y_planes):
#     for j in range(NX):
#         directory = './J/'+str(i)+'_'+str(j)
#         if not os.path.exists(directory):
#             os.makedirs(directory)
#         filename = 'J.txt'
#         np.savetxt(directory+'/'+filename, current_filter_y[i,j,:])

# print(Dtotal)
# CMFD_J_x, CMFD_J_y = condense_J(current_filter_x/Dtotal,current_filter_y/Dtotal,pitch,pitch)
# print(CMFD_J_x)
# print("------------------")
# print(CMFD_J_y)

# # 2 GROUP TEST
# k_eff = 1.2802
# scatter_term = dx*dy*np.sum(CMFD_phi[0,0,0]*CMFD_Sig_sct[0,0,0,1])
# print("chi:",np.sum(CMFD_Chi[0,0,:]))
# fission_term = dx*dy*CMFD_Chi[0,0,0]*np.sum(CMFD_Sig_nuf[0,0,:]*CMFD_phi[0,0,:])
# print("scatter_term:",scatter_term)
# print("fission_term:",fission_term)
# RHS = scatter_term+1/k_eff*fission_term
# print("RHS:",RHS)

# removal_term = dx*dy*CMFD_Sig_R[0,0,0]*CMFD_phi[0,0,0]
# print("removal_term:",removal_term)
# print("flux wannabe:",RHS-removal_term)
# Jimh = CMFD_J_x[0,0,0]
# print(Jimh)
# Jiph = CMFD_J_x[1,0,0]
# print(Jiph)
# Jjmh = CMFD_J_y[0,0,0]
# print(Jjmh)
# Jjph = CMFD_J_y[1,0,0]
# print(Jjph)
# stream_term = -dy*(Jimh-Jiph)-dx*(Jjmh-Jjph)
# print("stream_term:",stream_term)
# LHS = removal_term + stream_term
# print("LHS:",LHS)



# ###############################
# CMFD_Dhat_x,CMFD_Dhat_y = get_Dhat(CMFD_D)         
# CMFD_Dtil_x,CMFD_Dtil_y = get_Dtil(CMFD_phi,CMFD_Dhat_x, \
#                                     CMFD_Dhat_y,CMFD_J_x,CMFD_J_y)

# M = Set_M(CMFD_Dtil_x,CMFD_Dtil_y,CMFD_Dhat_x,CMFD_Dhat_y, \
#             CMFD_Sig_sct,CMFD_Sig_R)
# F =  Set_F(CMFD_Sig_nuf,CMFD_Chi)
# print("M")
# # print(M)
# filename='M_ref.txt'
# np.savetxt('./'+filename, M)
# plt.spy(M)
# plt.show()
# print("F")
# print(F)
# filename='F.txt'
# np.savetxt('./'+filename, F)

# phi_CM = CMFD_phi.reshape(NX*NY,NG)
# phi_CM = phi_CM.reshape((NX*NY*NG,1),order='F')
# phi_CM = phi_CM/np.linalg.norm(phi_CM)
# phi_CM = np.mat(phi_CM)

# kk,phi_CM_new = diffusion_solver(k,phi_CM,M,F)
# print(phi_CM_new)

#################### Plots##############
# if plot == 1:
#     ax=plt.gca()
#     pos_x_array = np.linspace(-(NX-1)/2,(NX-1)/2,NX)
#     pos_y_array = np.linspace((NY-1)/2,-(NY-1)/2,NY)


#     for i, xx in enumerate(pos_x_array):
#         for j,yy in enumerate(pos_y_array):
#             x0 = xx*pitch
#             y0 = yy*pitch
#             for r in rlist:
#                 circle = plt.Circle((x0, y0), r, color='r',fill=False)
#                 ax.add_artist(circle)

#     plt.axis('scaled')
#     ax.set_xlim((xbd_list[0]-0.5, xbd_list[-1]+0.5))
#     ax.set_ylim((ybd_list[-1]-0.5, ybd_list[0]+0.5))

#     # plot x planes
#     for i in range(NX+1):
#          plt.plot([-NX*pitch/2+i*pitch,-NX*pitch/2+i*pitch],\
#                     [ybd_list[0],ybd_list[-1]],'b-')

#     # plot y planes
#     for j in range(NY+1):
#          plt.plot([xbd_list[0],xbd_list[-1]],
#                     [NY*pitch/2-j*pitch,NY*pitch/2-j*pitch],'b-')

    # plt.show()
    # plt.savefig("plot.png",dpi=700)
