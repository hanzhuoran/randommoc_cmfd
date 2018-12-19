import math 
import numpy as np
import random

# Parameters to change
DZ = 10
MAX_D = 200
N_RAYS = 100
n_fuel_ring = 0
n_mod_ring = 0
N_GROUPS = 1 # 10 or 2 groups or 1
NX = 10
NY = 1
xmbd = "transreflect"
# xpbd = "transreflect"
# ymbd = "transreflect"
# ypbd = "transreflect"
# xmbd = "reflection"
xpbd = "reflection"
ymbd = "reflection"
ypbd = "reflection"
# ratio = 1


# constant
n_x_planes = (NX+1)
n_y_planes = (NY+1)
n_planes = (NX+1)+(NY+1)
n_rings = n_fuel_ring+n_mod_ring
N_CELL = NX*NY
N_REGION = N_CELL*(n_rings+1)
INF = math.inf
R = 0.392
pitch = 1.26
conv_crit1 = 1e-5
conv_crit2 = 1e-5
if N_GROUPS == 10:
    grp_bd = 6
elif N_GROUPS == 2:
    grp_bd = 1


##used in 2D diffusion in 2 GROUPS###
NG = 2 
dx = pitch
dy = dx
NX_surf = NX+1
NY_surf = NY+1
N_surf = NX_surf+NY_surf
xbd_list = np.linspace(-NX*pitch/2,NX*pitch/2,NX+1)
ybd_list = np.linspace(NY*pitch/2,-NY*pitch/2,NY+1)

if n_rings > 0:
    darea = np.pi*R**2/n_fuel_ring
    area = 0
    fuel_rlist = [0]
    for i in range(n_fuel_ring):
        area += darea
        r = np.sqrt(area/np.pi)
        fuel_rlist.append(r)
    fuel_rlist = np.asarray(fuel_rlist)
    mod_rlist = np.linspace(R,pitch/2,n_mod_ring+1)[1:]
    rlist = np.append(fuel_rlist[1:], mod_rlist)# without center point
    vlist = np.append(fuel_rlist, mod_rlist) # with center point
    print("r list:",rlist)
elif n_rings == 0:
    rlist = []
    vlist = []

##########################
sidelengthX = pitch*NX
sidelengthY = pitch*NY
sidelengthZ = sidelengthY
latticebox = [sidelengthX/2, sidelengthY/2, sidelengthZ/2,
             -sidelengthX/2, -sidelengthY/2, -sidelengthZ/2]
             
vol_cell = pitch**2*sidelengthZ
VTOT= sidelengthX*sidelengthY*sidelengthZ
STOT = (sidelengthX*sidelengthZ)*2+\
        (sidelengthX*sidelengthY)*2+\
        (sidelengthY*sidelengthY)*2
S_V_ratio = STOT/VTOT
theo_left = (sidelengthY*sidelengthZ)/STOT
theo_up = (sidelengthX*sidelengthZ)/STOT
theo_top = (sidelengthX*sidelengthY)/STOT
####################
XSdirectory = "./XS/"
homo_XSdirectory = "./homo_XS/"

def seed(sd=1):
	random.seed(sd)
	return

def uniform_dis():
	return random.uniform(0,1)

def circle_dis():
	return 2*np.pi*uniform_dis()

def mu_dis():
	return 2*uniform_dis()-1


