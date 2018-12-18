import numpy as np
from util import *


def find_D(g,i,j,D_mat_x,D_mat_y):
#        dym
# dxm          dxp
#        dyp
#
#
    # Find i,j+-1/2
    # Dy_index = i-1+NX_surf
    Dym = D_mat_y[j-1,i-1,g-1]
    Dyp = D_mat_y[j,i-1,g-1]

    # Find i+-1/2, j 
    Dxm = D_mat_x[i-1,j-1,g-1]
    Dxp = D_mat_x[i,j-1,g-1]
    return Dxm, Dxp, Dym, Dyp


def Set_F(CMFD_Sig_nuf,CMFD_Chi):
    NdimX = CMFD_Sig_nuf.shape[1]
    NdimY = CMFD_Sig_nuf.shape[0]
    n = NdimX*NdimY*NG
    m = n
    F = np.zeros((n,m))
    F = np.mat(F)
    g_list = np.linspace(1,NG,NG)
    i_list = np.linspace(1,NdimX,NdimX)
    j_list = np.linspace(1,NdimY,NdimY)

    for g in g_list.astype(int):
        for j in j_list.astype(int):
            for i in i_list.astype(int):
                INDEX = (g-1)*(NdimX*NdimY)+(j-1)*NdimX + (i-1)
                coeff = dx*dy*CMFD_Chi[j-1,i-1,g-1]
                for gp in g_list.astype(int):
                    # if gp != g:
                    F[INDEX,INDEX+(gp-g)*(NdimX*NdimY)] = coeff*\
                                                    CMFD_Sig_nuf[j-1,i-1,gp-1]
    return F


def Set_M(Dtil_x,Dtil_y,Dhat_x,Dhat_y,CMFD_Sig_S,CMFD_Sig_R):
# Create Matrix A
    NdimX = CMFD_Sig_R.shape[1]
    NdimY = CMFD_Sig_R.shape[0]
    # print(NdimX)
    # print(NdimY)
    n = NdimX*NdimY*NG
    m = n
    M = np.zeros((n,m))
    M = np.mat(M)
    g_list = np.linspace(1,NG,NG)
    i_list = np.linspace(1,NdimX,NdimX)
    j_list = np.linspace(1,NdimY,NdimY)

    for g in g_list.astype(int):
        for j in j_list.astype(int):
            for i in i_list.astype(int):
                INDEX = (g-1)*(NdimX*NdimY)+(j-1)*NdimX + (i-1)
                # print(g,j,i)
                # print("INDEX:",INDEX)
                Dhat_xm, Dhat_xp, Dhat_ym, Dhat_yp = find_D(g,i,j,Dhat_x,Dhat_y)
                Dtil_xm, Dtil_xp, Dtil_ym, Dtil_yp = find_D(g,i,j,Dtil_x,Dtil_y)
                #Process sactteing first
                for gp in g_list.astype(int):
                    if gp != g:
                        # order here might be wrong
                        M[INDEX, INDEX+(gp-g)*(NdimX*NdimY)] = -dx*dy*CMFD_Sig_S[j-1,i-1,g-1,gp-1]
                # Stencil
                #           V
                # W  X = A+B+C+D+R  Y
                #           Z
                #

                # if i == 1 and i != NdimX:
                #     W = 0
                #     Y = -dy*(Dhat_xp+Dtil_xp)
                #     B = 0
                #     C = dy*(Dhat_xp-Dtil_xp) 
                # elif i == NdimX  and i != 1:
                #     W = -dy*(Dhat_xm-Dtil_xm)
                #     Y = 0
                #     B = dy*(Dhat_xm+Dtil_xm)
                #     C = 0
                # elif i == 1 and i ==NdimX:
                #     W = 0
                #     Y = 0
                #     B = 0
                #     C = 0 
                # else:
                W = -dy*(Dhat_xm-Dtil_xm)
                Y = -dy*(Dhat_xp+Dtil_xp)
                B = dy*(Dhat_xm+Dtil_xm)
                C = dy*(Dhat_xp-Dtil_xp) 

                # if j == 1 and j != NdimY:
                #     V = 0
                #     Z = -dx*(Dhat_yp+Dtil_yp)
                #     A = 0
                #     D = dx*(Dhat_yp-Dtil_yp)
                # elif j == NdimY  and j != 1:
                #     V = -dx*(Dhat_ym-Dtil_ym)
                #     Z = 0
                #     A = dx*(Dhat_ym+Dtil_ym)
                #     D = 0
                # elif j == 1 and j == NdimY:
                #     V = 0
                #     Z = 0
                #     A = 0
                #     D = 0
                # else:
                V = -dx*(Dhat_ym-Dtil_ym)
                Z = -dx*(Dhat_yp+Dtil_yp)
                A = dx*(Dhat_ym+Dtil_ym)
                D = dx*(Dhat_yp-Dtil_yp)
                
                X = A+B+C+D+dx*dy*CMFD_Sig_R[j-1,i-1,g-1]
                M[INDEX,INDEX] = X

                if i != 1:
                    M[INDEX,INDEX-1] = W

                if i != NdimX:
                    M[INDEX,INDEX+1] = Y

                if j != 1:
                    M[INDEX,INDEX-NdimX] = V

                if j != NdimY:
                    M[INDEX,INDEX+NdimX] = Z

    return M

def get_Dhat(D_mat):
    NdimX = D_mat.shape[1]
    NdimY = D_mat.shape[0]
    # print(NdimX)
    # print(NdimY)
    Dhat_x = np.zeros((NdimX+1,NdimY,NG))
    Dhat_y = np.zeros((NdimY+1,NdimX,NG))

    g_list = (np.linspace(1,NG,NG)).astype(int)
    i_list = (np.linspace(1,NdimX,NdimX)).astype(int)
    j_list = (np.linspace(1,NdimY,NdimY)).astype(int)

    for g in g_list:
        # decide which yplanes
        if NdimY > 1:
            for j in j_list[0:-1]:
                #decide x index
                for i in i_list:
                    Dhat_y[j,i-1,g-1] = 2*D_mat[j-1,i-1,g-1]*D_mat[j,i-1,g-1] \
                                        /(dy*D_mat[j-1,i-1,g-1]+dy*D_mat[j,i-1,g-1])
        if NdimX > 1:  
        # decide which xplanes
            for i in i_list[0:-1]:
                #decide y index
                for j in j_list:
                    Dhat_x[i,j-1,g-1] = 2*D_mat[j-1,i-1,g-1]*D_mat[j-1,i,g-1] \
                                    /(dy*D_mat[j-1,i-1,g-1]+dy*D_mat[j-1,i,g-1])
    return Dhat_x, Dhat_y

def get_Dtil(phi,Dhat_x,Dhat_y,J_x,J_y):
    NdimX = phi.shape[1]
    NdimY = phi.shape[0]
    # print(NdimX)
    # print(NdimY)
    Dtil_x = np.zeros((NdimX+1,NdimY,NG))
    Dtil_y = np.zeros((NdimY+1,NdimX,NG))

    g_list = (np.linspace(1,NG,NG)).astype(int)
    i_list = (np.linspace(1,NdimX,NdimX)).astype(int)
    j_list = (np.linspace(1,NdimY,NdimY)).astype(int)

    for g in  g_list:
        # decide which y planes
        # if NdimY > 1:
        for j in range(NdimY+1):
            #decide x index
            # print("y planes")
            # print(j)
            if j != 0 and j != NdimY:
                for i in i_list:
                    Dtil_y[j,i-1,g-1] = (-Dhat_y[j,i-1,g-1]*(phi[j,i-1,g-1]- \
                                                            phi[j-1,i-1,g-1])- \
                                        J_y[j,i-1,g-1])/(phi[j,i-1,g-1]+
                                                            phi[j-1,i-1,g-1])
            elif j== 0:
                for i in i_list:
                    Dtil_y[0,i-1,g-1] = J_y[0,i-1,g-1]/phi[0,i-1,g-1]
            elif j == NdimY:
                for i in i_list:
                    Dtil_y[NdimY,i-1,g-1] = J_y[NdimY,i-1,g-1]/phi[NdimY-1,i-1,g-1]
        # decide which x planes
        # if NdimX > 1:
        for i in range(NdimX+1):
            #decide y index
            # print("x planes")
            # print(i)
            if i != 0 and i != NdimX:
                for j in j_list:
                    Dtil_x[i,j-1,g-1] = (-Dhat_x[i,j-1,g-1]*(phi[j-1,i,g-1]- \
                                                            phi[j-1,i-1,g-1])- \
                                        J_x[i,j-1,g-1])/(phi[j-1,i,g-1]+
                                                            phi[j-1,i-1,g-1])
            elif i == 0:
                for j in j_list:
                    Dtil_x[0,j-1,g-1] = J_x[0,j-1,g-1]/phi[j-1,0,g-1]

            elif i == NdimX:
                for j in j_list:
                    Dtil_x[NdimX,j-1,g-1] = J_x[NdimX,j-1,g-1]/phi[j-1,NdimX-1,g-1]

    return Dtil_x, Dtil_y


def condense_xs(Regions, CMFD_Sig_tot,CMFD_Sig_abs, CMFD_Sig_nuf, CMFD_Chi, \
                CMFD_D, CMFD_Sig_sct, CMFD_phi):

    CMFD_phiTot = np.zeros((NY,NX,NG))
    CMFD_phiAbs = np.zeros((NY,NX,NG))
    CMFD_phiNuF = np.zeros((NY,NX,NG))
    CMFD_phiV = np.zeros((NY,NX,NG))
    CMFD_phiD = np.zeros((NY,NX,NG))
    CMFD_V = np.zeros((NY,NX))
    CMFD_Chi_up = np.zeros((NY,NX,NG))
    CMFD_Chi_down = np.zeros((NY,NX))
    CMFD_phiSct = np.zeros((NY,NX,NG,NG))

    # How to calc D?
    for y in range(NY):
        for x in range(NX):            
            for c in range(n_rings+1):
                # print("coordinate:",y,x,c)
                region = Regions[y,x,c]
                #####Condensation#####
                # sum phi*A (phi*V)
                CMFD_V[y,x] += region.volume
                for g in range(N_GROUPS):
                    if g < grp_bd:
                        # fast group
                        CMFD_phiV[y,x,0] += region.phi[0,g]* region.volume
                        CMFD_phiTot[y,x,0] += region.Sig_tot[g,0]* \
                                            region.phi[0,g]*region.volume
                        CMFD_phiAbs[y,x,0] += region.Sig_abs[g,0]* \
                                            region.phi[0,g]*region.volume
                        CMFD_phiNuF[y,x,0] += region.Sig_nuf[g,0]* \
                                            region.phi[0,g]*region.volume
                        CMFD_phiD[y,x,0] += (1/(3*region.Sig_tot[g,0]))* \
                                            region.phi[0,g]*region.volume
                        for gp in range(N_GROUPS):
                            CMFD_Chi_up[y,x,0] += region.Chi[g,0]*\
                                                region.Sig_nuf[gp,0]* \
                                                region.phi[0,gp]*region.volume
                            if gp < grp_bd:
                                CMFD_phiSct[y,x,0,0] += region.Sig_sct[gp,g]* \
                                                        region.phi[0,g]* \
                                                        region.volume
                            elif gp >= grp_bd:  
                                CMFD_phiSct[y,x,1,0] += region.Sig_sct[gp,g]* \
                                                        region.phi[0,g]* \
                                                        region.volume
                    elif g >= grp_bd:
                        # thermal group
                        CMFD_phiV[y,x,1] += region.phi[0,g]* region.volume
                        CMFD_phiTot[y,x,0] += region.Sig_tot[g,0]* \
                                            region.phi[0,g]*region.volume
                        CMFD_phiAbs[y,x,1] += region.Sig_abs[g,0]* \
                                            region.phi[0,g]*region.volume
                        CMFD_phiNuF[y,x,1] += region.Sig_nuf[g,0]* \
                                            region.phi[0,g]*region.volume
                        CMFD_phiD[y,x,1] += (1/(3*region.Sig_tot[g,0]))* \
                                            region.phi[0,g]*region.volume
                        for gp in range(N_GROUPS):
                            CMFD_Chi_up[y,x,1] += region.Chi[g,0]* \
                                                region.Sig_nuf[gp,0]* \
                                                region.phi[0,gp]*region.volume
                            if gp < grp_bd:
                                CMFD_phiSct[y,x,0,1] += region.Sig_sct[gp,g]*\
                                                        region.phi[0,g]*\
                                                        region.volume
                            elif gp >= grp_bd:  
                                CMFD_phiSct[y,x,1,1] += region.Sig_sct[gp,g]*\
                                                        region.phi[0,g]* \
                                                        region.volume

                    for gpp in range(N_GROUPS):
                        CMFD_Chi_down[y,x] += region.Chi[g,0]* \
                                        region.Sig_nuf[gpp,0]*\
                                        region.phi[0,gpp]*region.volume

            CMFD_phi[y,x,:] = CMFD_phiV[y,x,:]/CMFD_V[y,x]
            CMFD_Sig_tot[y,x,:] = np.divide(CMFD_phiTot[y,x,:],CMFD_phiV[y,x,:])
            CMFD_Sig_abs[y,x,:] = np.divide(CMFD_phiAbs[y,x,:],CMFD_phiV[y,x,:])
            CMFD_Sig_nuf[y,x,:] = np.divide(CMFD_phiNuF[y,x,:],CMFD_phiV[y,x,:])
            CMFD_Chi[y,x,:] = CMFD_Chi_up[y,x,:]/CMFD_Chi_down[y,x]
            CMFD_D[y,x,:] = np.divide(CMFD_phiD[y,x,:],CMFD_phiV[y,x,:])
            CMFD_Sig_sct[y,x,:,:] = np.divide(CMFD_phiSct[y,x,:,:],\
                                    CMFD_phiV[y,x,:])

def condense_J(current_filter_x,current_filter_y,Delta_x,Delta_y):

    # current_filter_x /= Dtotal
    # current_filter_y /= Dtotal
    n_x_planes_dim = current_filter_x.shape[0]
    NdimY = current_filter_x.shape[1]
    n_y_planes_dim = current_filter_y.shape[0]
    NdimX = current_filter_y.shape[1]
    CMFD_J_x = np.zeros((n_x_planes_dim,NdimY,NG))
    CMFD_J_y = np.zeros((n_y_planes_dim,NdimX,NG))
    # print(CMFD_J_x.shape)
    # print(CMFD_J_x.shape)
    for i in range(n_x_planes_dim):
        for j in range(NdimY):
            CMFD_J_x[i,j,0] = np.sum(current_filter_x[i,j,0:grp_bd])
            CMFD_J_x[i,j,1] = np.sum(current_filter_x[i,j,grp_bd:N_GROUPS])
    
    for i in range(n_y_planes_dim):
        for j in range(NdimX):
            CMFD_J_y[i,j,0] = np.sum(current_filter_y[i,j,0:grp_bd])
            CMFD_J_y[i,j,1] = np.sum(current_filter_y[i,j,grp_bd:N_GROUPS])

    CMFD_J_x = CMFD_J_x/Delta_y
    CMFD_J_y = CMFD_J_y/Delta_x

    return CMFD_J_x,CMFD_J_y

def calc_sig_r(CMFD_Sig_abs,CMFD_Sig_sct):
    # NGG = CMFD_Sig_abs.shape[2]
    NdimX = CMFD_Sig_abs.shape[1]
    NdimY = CMFD_Sig_abs.shape[0]
    # print(NdimX)
    # print(NdimY)
    CMFD_Sig_R = np.zeros((NdimY,NdimX,NG))
    for x in range(NdimX):
        for y in range(NdimY):
            for g in range(NG):
                CMFD_Sig_R[y,x,g] = CMFD_Sig_abs[y,x,g]
                for gp in range(NG):
                    if gp != g:
                        CMFD_Sig_R[y,x,g] += CMFD_Sig_sct[y,x,gp,g]
    return CMFD_Sig_R

def check_convergence(phi,k,phitemp,ktemp,M):
    converged = 0
    tol1 = 1e-5
    tol2 = 1e-5
    tol3 = 1e-5 

    # k
    r1 = np.abs(k-ktemp)/np.abs(k)

    ###flux####

    r2 = np.linalg.norm((phi - phitemp)/phi)

    ####fission###
    Fphi = M*phi
    Fphitemp = M*phitemp

    r3 = (np.linalg.norm(Fphi-Fphitemp))/np.linalg.norm(Fphi)

    if  r1 < tol1 and r2 < tol2 and r3 < tol3:
        converged = 1
    return converged


# def spatial_condense_xs(Regions, \
#         CMFD_Sig_abs, CMFD_Sig_nuf, CMFD_Chi, CMFD_D, CMFD_Sig_sct, CMFD_phi):

#     CMFD_phiAbs = np.zeros((NY,NX,N_GROUPS))
#     CMFD_phiNuF = np.zeros((NY,NX,N_GROUPS))
#     CMFD_phiV = np.zeros((NY,NX,N_GROUPS))
#     CMFD_phiD = np.zeros((NY,NX,N_GROUPS))
#     # CMFD_phiChi = np.zeros((NY,NX,N_GROUPS))
#     CMFD_Chi_up = np.zeros((NY,NX,N_GROUPS))
#     CMFD_Chi_down = np.zeros((NY,NX))
#     CMFD_V = np.zeros((NX,N_GROUPS))
#     # CMFD_Chi_up = np.zeros((NY,NX,N_GROUPS))
#     # CMFD_Chi_down = np.zeros((NY,NX))
#     CMFD_phiSct = np.zeros((NY,NX,N_GROUPS,N_GROUPS))
#     # How to calc D?
#     for x in range(NX):
#         for y in range(NY):
#             for c in range(n_rings+1):
#                 region = Regions[y,x,c]
#                 CMFD_V[y,x] += region.volume
#                 for g in range(N_GROUPS):
#                     CMFD_phiV[y,x,g] += region.phi[0,g]* region.volume
#                     CMFD_phiAbs[y,x,g] += region.Sig_abs[g,0]* \
#                                         region.phi[0,g]*region.volume
#                     CMFD_phiNuF[y,x,g] += region.Sig_nuf[g,0]* \
#                                         region.phi[0,g]*region.volume
#                     CMFD_phiD[y,x,g] += (1/(3*region.Sig_tot[g,0]))* \
#                                         region.phi[0,g]*region.volume
#                     # CMFD_phiChi[y,x,g] += region.Chi[g,0]* \
#                     #                     region.phi[0,g]*region.volume
#                     for gp in range(N_GROUPS):
#                         CMFD_phiSct[y,x,gp,g] += region.Sig_sct[gp,g]*\
#                                                 region.phi[0,g]* \
#                                                 region.volume
#                         CMFD_Chi_up[y,x,g] += region.Chi[g,0]* \
#                                                 region.Sig_nuf[gp,0]* \
#                                                 region.phi[0,gp]*region.volume
#                     for gpp in range(N_GROUPS):
#                         CMFD_Chi_down[y,x] += region.Chi[g,0]* \
#                                         region.Sig_nuf[gpp,0]*\
#                                         region.phi[0,gpp]*region.volume


#             CMFD_phi[y,x,:] = CMFD_phiV[y,x,:]/CMFD_V[y,x]
#             CMFD_Sig_abs[y,x,:] = np.divide(CMFD_phiAbs[y,x,:],CMFD_phiV[y,x,:])
#             CMFD_Sig_nuf[y,x,:] = np.divide(CMFD_phiNuF[y,x,:],CMFD_phiV[y,x,:])
#             CMFD_Chi[y,x,:] = CMFD_Chi_up[y,x,:]/CMFD_Chi_down[y,x]
#             CMFD_D[y,x,:] = np.divide(CMFD_phiD[y,x,:],CMFD_phiV[y,x,:])
#             CMFD_Sig_sct[y,x,:,:] = np.divide(CMFD_phiSct[y,x,:,:],\
#                                     CMFD_phiV[y,x,:])

def diffusion_solver(k,phi_CM,M,F):
    phi_CM_new = phi_CM
    converged = 0
    iteration = 0
    while not converged:
        iteration += 1
        # print("CMFD iter:", iteration)
        if iteration >= 1000:
            # print("warning!")
            break
        b = 1/k*F*phi_CM_new
        phi_CM_temp = phi_CM_new
        phi_CM_new = np.linalg.solve(M, b)
        # print("phi_CM")
        # print(phi_CM_new)
        ktemp = k
        k = np.linalg.norm(F*phi_CM_new)/np.linalg.norm(F*phi_CM_temp)*ktemp
        # print("k")
        # print(k)
        phi_CM_new = phi_CM_new/np.linalg.norm(phi_CM_new)
        converged = check_convergence(phi_CM_new,k,phi_CM_temp,ktemp,F)
    print("diffusion k:", k)

    return k,phi_CM_new

def update_CMFD_phi(factor,Regions):
    factor = np.squeeze(np.asarray(factor))
    f_pro = factor.reshape(NX*NY,NG,order='F')
    f_pro = f_pro.reshape(NY,NX,NG)
    
    # print(f_pro)

    norm_pre = np.zeros((NY,NX,n_rings+1))
    for j in range (NY):
        for i in range (NX):
            for c in range(n_rings+1):
                norm_pre[j,i,c] = np.linalg.norm(Regions[j,i,c].phi)
                for g in range(N_GROUPS):
                    if g < grp_bd:
                        Regions[j,i,c].phi[0,g] *= f_pro[j,i,0]
                    elif g >= grp_bd:
                        # print(f_pro[i,j,1])
                        Regions[j,i,c].phi[0,g] *= f_pro[j,i,1]
                # if i == 0 and j ==0 and c ==0:
                #     print(f_pro[i,j,0])
                #     print(Regions[i,j,c].phi)
    return norm_pre