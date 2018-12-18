from ray import Ray
import numpy as np
from geom import *
import math
from util import *
from moc import *
import matplotlib.pyplot as plt
import copy 
def process_homo_material(directory, groups, Regions):
    full_dir = directory+str(groups)

    ###SCATTER_XS
    sct = np.empty((N_GROUPS,N_GROUPS))
    index = 0
    with open(full_dir+"/scatter","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            sct[index,:] =a
            index +=1

    ###TOTAL_XS
    tot = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/total","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            tot[index,:] =a
            index +=1

    Abs = np.empty((N_GROUPS,1))
    for i in range(N_GROUPS):
        Abs[i,0] = tot[i,0] - np.sum(sct[:,i])

    rmv = np.empty((N_GROUPS,1))
    for i in range(N_GROUPS):
        rmv[i,0] = tot[i,0] - np.sum(sct[i,i])


    ###NUF_XS
    nuf = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/nufission","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            nuf[index,:] =a
            index +=1

    ###CHI_XS
    chi = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/chi","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            chi[index,:] =a
            index +=1
            
    for j in range(NY):
        for i in range(NX): 
            for k in range(n_rings+1):
                    Regions[j,i,k].set_tot(tot)
                    Regions[j,i,k].set_sct(sct)
                    Regions[j,i,k].set_nuf(nuf)
                    Regions[j,i,k].set_chi(chi)
                    Regions[j,i,k].set_abs(Abs)
                    Regions[j,i,k].set_rmv(rmv)
                    Regions[j,i,k].NF = np.multiply(Regions[j,i,k].Sig_nuf.transpose(),
                                    Regions[j,i,k].phi)*pitch**2

def process_material(directory, groups, Regions):
    full_dir = directory+str(groups)
    ###FUEL_SCATTER_XS
    fuel_sct = np.empty((N_GROUPS,N_GROUPS))
    index = 0
    with open(full_dir+"/scatter_fuel","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            fuel_sct[index,:] =a
            index +=1

    ###FUEL_TOTAL_XS
    fuel_tot = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/total_fuel","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            fuel_tot[index,:] =a
            index +=1

    fuel_abs = np.empty((N_GROUPS,1))
    for i in range(N_GROUPS):
        fuel_abs[i,0] = fuel_tot[i,0] - np.sum(fuel_sct[:,i])

    fuel_rmv = np.empty((N_GROUPS,1))
    for i in range(N_GROUPS):
        fuel_rmv[i,0] = fuel_tot[i,0] - np.sum(fuel_sct[i,i])


    ###FUEL_NUF_XS
    fuel_nuf = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/nufission_fuel","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            fuel_nuf[index,:] =a
            index +=1

    ###FUEL_CHI_XS
    fuel_chi = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/chi_fuel","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            fuel_chi[index,:] =a
            index +=1
    
    ###mod_SCATTER_XS
    mod_sct = np.empty((N_GROUPS,N_GROUPS))
    index = 0
    with open(full_dir+"/scatter_mod","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            mod_sct[index,:] =a
            index +=1

    ###mod_TOTAL_XS
    mod_tot = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/total_mod","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            mod_tot[index,:] =a
            index +=1

    mod_abs = np.empty((N_GROUPS,1))
    for i in range(N_GROUPS):
        mod_abs[i,0] = mod_tot[i,0] - np.sum(mod_sct[:,i])

    mod_rmv = np.empty((N_GROUPS,1))
    for i in range(N_GROUPS):
        mod_rmv[i,0] = mod_tot[i,0] - np.sum(mod_sct[i,i])

    ###mod_NUF_XS
    mod_nuf = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/nufission_mod","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            mod_nuf[index,:] =a
            index +=1

    ###mod_CHI_XS
    mod_chi = np.empty((N_GROUPS,1))
    index = 0
    with open(full_dir+"/chi_mod","r") as file:
        for line in file:
            a = line.split()
            a = list(map(float, a))
            mod_chi[index,:] =a
            index +=1
            
    for j in range(NY):
        for i in range(NX): 
            for k in range(n_rings+1):
                if k < n_fuel_ring:
                    Regions[j,i,k].set_tot(fuel_tot)
                    Regions[j,i,k].set_sct(fuel_sct)
                    Regions[j,i,k].set_nuf(fuel_nuf)
                    Regions[j,i,k].set_chi(fuel_chi)
                    Regions[j,i,k].set_abs(fuel_abs)
                    Regions[j,i,k].set_rmv(fuel_rmv)
                    Regions[j,i,k].NF = np.multiply(Regions[j,i,k].Sig_nuf.transpose(),
                                    Regions[j,i,k].phi)*pitch**2
                elif k >= n_fuel_ring:
                    Regions[j,i,k].set_tot(mod_tot)
                    Regions[j,i,k].set_sct(mod_sct)
                    Regions[j,i,k].set_nuf(mod_nuf)
                    Regions[j,i,k].set_chi(mod_chi)
                    Regions[j,i,k].set_abs(mod_abs)
                    Regions[j,i,k].set_rmv(mod_rmv)
