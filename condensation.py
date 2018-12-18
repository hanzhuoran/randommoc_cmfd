import numpy as np
import math
import matplotlib.pyplot as plt
import copy
import csv
import os

from ray import Ray
from geom import *
from moc import *
from material import *
from CMFD_func import *
from set_region import *
from util import *

np.random.seed(1)

latticebox = [sidelengthX/2, sidelengthY/2, sidelengthZ/2,
             -sidelengthX/2, -sidelengthY/2, -sidelengthZ/2]
# cellbox = [pitch/2,pitch/2,pitch/2,-pitch/2,-pitch/2,-pitch/2]
Regions,all_surf = set_geometry()

for y in range(NY):
    for x in range(NX):
        for c in range(n_rings+1):
            # print("--------------")
            # print("xyc",x,y,c)
            if c < n_rings:
                volume = np.pi*(vlist[c+1]**2-vlist[c]**2)
            else:
                volume = pitch**2-np.pi*vlist[c]**2

            directory = './rr_phi/'+str(y)+'_'+str(x)+'_'+str(c)
            filename='phi.txt'
            with open(directory+"/"+filename,"r") as file:
                for line in file:
                    # print("line:",line)
                    phi = line.split()
                    phi = list(map(float, phi))
            Regions[y,x,c].set_phi(phi)

            # filename='volume.txt'
            # with open(directory+"/"+filename,"r") as file:
            #     for line in file:
            #         # print("line:",line)
            #         volume = line.split()
            #         volume = float(volume[0])
            # Regions[x,y,c].set_phi(phi)
            Regions[y,x,c].set_volume(volume)


#  

# CMFD_J = np.zeros((N_surf,3,NG))
# for i in range(n_planes):
#     for j in range(NX):
#         CMFD_J[i,j,0] = np.sum(J_raw[i,j,0:6])
#         CMFD_J[i,j,1] = np.sum(J_raw[i,j,6:10])
# print("CMFD_J", CMFD_J)


# #10 GROUP CMFD
# # CMFD_phi = np.zeros((NX,NY,N_GROUPS))
# # CMFD_Sig_abs = np.zeros((NX,NY,N_GROUPS)) 
# # CMFD_Sig_nuf = np.zeros((NX,NY,N_GROUPS))
# # CMFD_Chi = np.zeros((NX,NY,N_GROUPS))
# # CMFD_D = np.zeros((NX,NY,N_GROUPS))
# # CMFD_Sig_sct = np.zeros((NX,NY,N_GROUPS,N_GROUPS))
# # spatial_condense_xs(Regions, CMFD_Sig_abs, CMFD_Sig_nuf, CMFD_Chi, CMFD_D, CMFD_Sig_sct,\
# #              CMFD_phi)

# #2 GROUP CMFD
# CMFD_phi = np.zeros((NX,NY,NG))
# CMFD_Sig_abs = np.zeros((NX,NY,NG)) 
# CMFD_Sig_nuf = np.zeros((NX,NY,NG))
# CMFD_Chi = np.zeros((NX,NY,NG))
# CMFD_D = np.zeros((NX,NY,NG))
# CMFD_Sig_sct = np.zeros((NX,NY,NG,NG))
# condense_xs(Regions, CMFD_Sig_abs, CMFD_Sig_nuf, CMFD_Chi, CMFD_D, CMFD_Sig_sct,\
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

# # k_eff = 1.2818001672645503


# # 2 GROUP TEST
# # scatter_term = dx*dy*np.sum(CMFD_phi[1,1,0]*CMFD_Sig_sct[1,1,0,1])
# # print("chi:",np.sum(CMFD_Chi[1,1,:]))
# # fission_term = dx*dy*CMFD_Chi[1,1,0]*np.sum(CMFD_Sig_nuf[1,1,:]*CMFD_phi[1,1,:])
# # print("scatter_term:",scatter_term)
# # print("fission_term:",fission_term)
# # RHS = scatter_term+1/k_eff*fission_term
# # print("RHS:",RHS)

# # removal_term = dx*dy*CMFD_Sig_R[1,1,0]*CMFD_phi[1,1,0]
# # print("removal_term:",removal_term)
# # Jimh = J_raw[1,1,0]
# # print(Jimh)
# # Jiph = J_raw[2,1,0]
# # print(Jiph)
# # Jjmh = J_raw[5,1,0]
# # print(Jjmh)
# # Jjph = J_raw[6,1,0]
# # print(Jjph)
# # stream_term = -dy*(Jimh-Jiph)-dx*(Jjmh-Jjph)
# # print("stream_term:",stream_term)
# # LHS = removal_term + stream_term
# # print("LHS:",LHS)


# # 10 GROUP TEST
# # scatter_term = dx*dy*np.sum(CMFD_phi[1,1,0:9]*CMFD_Sig_sct[1,1,9,0:9])
# # print("chi:",np.sum(CMFD_Chi[1,1,:]))
# # fission_term = dx*dy*CMFD_Chi[1,1,9]*np.sum(CMFD_Sig_nuf[1,1,:]*CMFD_phi[1,1,:])
# # print("scatter_term:",scatter_term)
# # print("fission_term:",fission_term)
# # RHS = scatter_term+1/k_eff*fission_term
# # print("RHS:",RHS)
# # removal_term = dx*dy*CMFD_Sig_R[1,1,9]*CMFD_phi[1,1,9]
# # print("removal_term:",removal_term)
# # Jimh = J_raw[1,1,9]
# # print(Jimh)
# # Jiph = J_raw[2,1,9]
# # print(Jiph)
# # Jjmh = J_raw[5,1,9]
# # print(Jjmh)
# # Jjph = J_raw[6,1,9]
# # print(Jjph)
# # stream_term = -dy*(Jimh-Jiph)-dx*(Jjmh-Jjph)
# # print("stream_term:",stream_term)
# # LHS = removal_term + stream_term
# # print("LHS:",LHS)

# ##Calc D hat on surface
# CMFD_Dhat = get_Dhat(CMFD_D)
# print("CMFD_Dhat")
# print(CMFD_Dhat)

# # CMFD_Dtil = np.zeros((N_surf,3,NG))
# CMFD_Dtil = get_Dtil(CMFD_phi,CMFD_Dhat,CMFD_J)
# print("CMFD_Dtil")
# print(CMFD_Dtil)


# #LHS matrix
# A = Set_A(CMFD_Dtil,CMFD_Dhat,CMFD_Sig_sct,CMFD_Sig_R)
# A = np.mat(A)
# print("A")
# print(A)
# filename='A_new.txt'
# np.savetxt('./'+filename,A)


# #RHS matrix
# M =  Set_M(CMFD_Sig_nuf,CMFD_Chi)
# M = np.mat(M)
# print("M")
# print(M)
# filename='M.txt'
# np.savetxt('./'+filename, M)


# # Sam's method
# # phi = CMFD_phi.reshape(9,2)
# # phi = phi.reshape((18,1),order='F')
# # phi = np.mat(phi)
# # phi = np.zeros((18,1))
# # phi[0,0] = 1
# # phi = np.mat(phi)
# # print("init phi:",phi)

# # S_old = M*phi
# # S_old /= np.sum(S_old)/(NX*NY*NG)

# # converged = 0
# # iteration = 0
# # k_eff = 1
# # while not converged:
# #     iteration +=1
# #     if iteration >= 20:
# #         break 
# #     phi = np.linalg.solve(A, S_old)
# #     print(phi)
# #     S_new = M*phi
# #     k_pre = k_eff
# #     k_eff = np.sum(S_new)/np.sum(S_old)
# #     print("k:",k_eff)
# #     S_old *= k_eff
# #     # R = (S_new-S_old)/S_new
# #     # EPS = np.sqrt(np.sum(np.square(R))/(NX*NY*NG))
# #     if  abs(k_eff-k_pre)<= 1e-4:
# #         converged =1
# #     S_new /= np.sum(S_new)/(NX*NY*NG)
# #     S_old = S_new




# phi = CMFD_phi.reshape(9,2)
# phi = phi.reshape((18,1),order='F')
# phi = phi/np.linalg.norm(phi)
# phi = np.mat(phi)
# print("init phi:",phi)

# k = 1
# converged = 0
# iteration = 0

# while not converged:
#     iteration += 1
#     if iteration >= 100:
#         break
#     b = 1/k*M*phi
#     phitemp = phi
#     phi = np.linalg.solve(A, b)
#     print("phi")
#     print(phi)
#     ktemp = k
#     k = np.linalg.norm(M*phi)/np.linalg.norm(M*phitemp)*ktemp
#     print("k")
#     print(k)
#     phi = phi/np.linalg.norm(phi)
#     converged = check_convergence(phi,k,phitemp,ktemp,M)

# print(iteration)

# # phi = CMFD_phi.reshape(9,2)
# # phi = phi.reshape((18,1),order='F')
# # phi = np.mat(phi)
# # print(phi)
# # RHS = M*phi
# # LHS = A*phi
# # print("RHS")
# # print(RHS)
# # print("LHS")
# # print(LHS)
# # print(np.divide(RHS,LHS))

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

plt.show()
plt.savefig("plot.png",dpi=500)