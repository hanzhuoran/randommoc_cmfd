from ray import Ray
import numpy as np
from geom import *
import math
from util import *

def attenuate_segment(ray,region,s):
    # region.cnt = region.cnt+1
    # print("ATTEN")
    for i in range(N_GROUPS):
        d_psi = (ray.psi[i]-region.q[i])* \
                (1-np.exp(-region.Sig_tot[i,0]*s))
        region.phi[0,i] = region.phi[0,i] + (4*np.pi*d_psi)
        ray.psi[i] = ray.psi[i]-d_psi

def DeadZone(ray,region,s):
    # print("DZ")
    # region.cnt = region.cnt+1
    for i in range(N_GROUPS):
        d_psi = (ray.psi[i]-region.q[i])* \
                (1-np.exp(-region.Sig_tot[i,0]*s))
        ray.psi[i] = ray.psi[i]-d_psi
