import numpy as np
from ray import *
import math
from util import *


class Surface:
    def __init__(self,id=None):
        self.id = id

class Plane(Surface):
    def __init__(self, A, B, C, D,type="transmission",id=None):
        super().__init__(id)
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.type = type
        self.cnt = 0
        # the other type is "reflection"
    def get_hit(self):
        return self.hit
    def reset_hit(self):
        self.hit = 0
    def getA(self):
        return self.A
    def getB(self):
        return self.B
    def getC(self):
        return self.C
    def getD(self):
        return self.D

class Z_Cylinder(Surface):
    def __init__(self,x0,y0,R,type="transmission",id=None):
        super().__init__(id)
        self.x0 = x0
        self.y0 = y0
        self.R = R
        self.type = type
    def get_center(self):
        return np.asarray([self.x0, self.y0])
    def get_R(self):
        return self.R

# 0 and inf rings
c_in = Z_Cylinder(0,0,0)
c_out = Z_Cylinder(0,0,math.inf)

def determine_boxlist(xbd,ybd):
    boxlist = np.empty((NY,NX,4),dtype=Surface)
    for j in range(NY):
        for i in range(NX):
            #    3 
            #  1   2 
            #    4
            # 
            boxlist[j,i] = np.asarray([xbd[i],xbd[i+1],ybd[j],ybd[j+1]])
    return boxlist


class Region:
    # Ad hoc
    def __init__(self, surf_list = None, surf_dir=None, id = None):

        self.surf_list = np.asarray(surf_list)
        self.surf_dir = np.asarray(surf_dir)
        self.id = id
        self.volume = 0
        self.Sig_tot = None
        self.Sig_nuf = None
        self.Sig_sct = None
        self.Sig_abs = None
        self.Sig_rmv = None
        self.NF = 0
        self.Chi = None
        self.phi = np.mat(np.ones(N_GROUPS))
        self.phi_old = np.mat(np.ones(N_GROUPS))
        self.q = np.zeros(N_GROUPS)
        self.length = 0
        self.cnt = 0
        self.physic_vol = 0

    def calc_volume(self,Dtotal):
        self.volume = self.length/Dtotal

    def set_volume(self,volume):
        self.volume = volume

    def reset_length(self):
        self.length = 0

    def reset_phi(self):
        self.phi = np.mat(np.zeros(N_GROUPS))

    def set_tot(self,tot):
        self.Sig_tot = np.mat(tot)

    def set_sct(self,sct):
        self.Sig_sct = np.mat(sct)

    def set_nuf(self,nuf):
        self.Sig_nuf = np.mat(nuf)

    def set_abs(self,absxs):
        self.Sig_abs = np.mat(absxs)

    def set_rmv(self,rmvxs):
        self.Sig_rmv = np.mat(rmvxs)

    def set_chi(self,chi):
        self.Chi = np.mat(chi)
            
    def set_phi(self, phi):
        self.phi = np.mat(phi)

    def set_nf(self,nf):
        self.NF = np.mat(nf)

    def calc_q(self, k):
        q = np.empty(N_GROUPS)

        for group in range(N_GROUPS):
            # This is for a single energy group at a certain region
            # all the xs should be a 1D numpy matrix.
            # sct is g' to g
            S = np.sum(self.Sig_sct[group,:]*(self.phi.transpose()))
            # nuf is the xs at g', which is a single value
            # chi is the chi at group g, which is also a single value
            F = self.Chi[group,0]*np.sum(self.Sig_nuf.transpose()* \
                                        (self.phi.transpose()))
            q[group] = 1/(4*np.pi*self.Sig_tot[group,0])*(S+F/k)
        self.q = q


def dis2cyl(cyl,ray):
    
    x = ray.x
    y = ray.y
    z = ray.z
    u = ray.u
    v = ray.v
    w = ray.w
    x0 = cyl.x0
    y0 = cyl.y0
    R = cyl.R
    x1 = x-x0
    y1 = y-y0
    a = u**2+v**2
    k = x1*u+y1*v
    c = x1**2+y1**2-R**2

    if ( a == 0. or k**2-a*c <0):
        d = INF
    elif (c < 0):
        d = (-k+np.sqrt(k**2-a*c))/a
    elif (c > 0):
        if ((-k-np.sqrt(k**2-a*c))/a >0):
            d = (-k-np.sqrt(k**2-a*c))/a
        else: 
            d = INF
    else: 
        d = INF 
    return d 


def dis2plane(plane, ray):
    A = plane.A
    B = plane.B
    C = plane.C
    D = plane.D
    x = ray.x
    y = ray.y
    z = ray.z
    u = ray.u
    v = ray.v
    w = ray.w
    d = (D-A*x-B*y-C*z)/(A*u+B*v+C*w)
    return d

def find_region(ray,Regions):
    here = 0
    x_dim = np.shape(Regions)[1]
    y_dim = np.shape(Regions)[0]
    refine_dim = np.shape(Regions)[2]

    x = ray.x
    y = ray.y
    done = False

    pos_x_array = np.linspace(-(NX-1)/2,(NX-1)/2,NX)
    pos_y_array = np.linspace((NY-1)/2,-(NY-1)/2,NY)

    for i, xx in enumerate(pos_x_array):
        for j,yy in enumerate(pos_y_array):
            xL = (xx-1/2)*pitch
            xR = (xx+1/2)*pitch
            yT = (yy+1/2)*pitch
            yB = (yy-1/2)*pitch
            if x>xL and x<=xR and y>yB and y<=yT :
                done = True
                break
        if done:
                break

    done = False
    for k in range(refine_dim):
        surf_list = Regions[j,i,k].surf_list
        # surf_dir = Regions[j,i,k].surf_dir
        if np.size(surf_list) == 1:
            x0 = surf_list[0].x0
            y0 = surf_list[0].y0
            r = surf_list[0].R
            if ((x-x0)**2+(y-y0)**2 < r**2):
                break
        elif np.size(surf_list) == 2:
            x0 = surf_list[0].x0
            y0 = surf_list[0].y0
            r1 = surf_list[0].R
            r2 = surf_list[1].R
            if ((x-x0)**2+(y-y0)**2 > r1**2) and ((x-x0)**2+(y-y0)**2 < r2**2):
                break
    return j,i,k

def pos_plane(plane,ray):
    A = plane.A
    B = plane.B
    C = plane.C
    D = plane.D
    x = ray.x
    y = ray.y
    z = ray.z
    if A*x+b*y+c*z > D:
        return 1
    else: 
        return -1


def tally_current(pln,ray,current_filter_x,current_filter_y):

    index = 0

    if pln.B ==0 and pln.A != 0:
        #  x plane, vertical, use ybd_list
        # cos = ray.u
        # print("ray.psi")
        # print(ray.psi)
        ID = pln.id
        for i in range(NY):
            if ray.y <= ybd_list[i] and ray.y > ybd_list[i+1]:
                index = i
                break
            if ray.y == ybd_list[NY]:
                index = NY
        current_filter_x[ID,index,:] = current_filter_x[ID,index,:] + \
                                        ray.psi*ray.u*2*np.pi
        # print("current_filter_x")
        # print(current_filter_x)
        # if ray.u > 0:
        #     current_filter_x[ID,index,:] = current_filter_x[ID,index,:] + 2*np.pi*ray.psi
        # else:
        #     current_filter_x[ID,index,:] = current_filter_x[ID,index,:] - 2*np.pi*ray.psi
    elif pln.A == 0 and pln.B != 0:
        #  y plane, horizontal, use xbd_list
        #  cos = ray.v
        ID = pln.id - (NX+1)
        for i in range(NX):
            if ray.x >= xbd_list[i] and ray.x < xbd_list[i+1]:
                index = i
                break
            if ray.x == xbd_list[NX]:
                index = NX
        current_filter_y[ID,index,:] = current_filter_y[ID,index,:] + \
                                        ray.psi*ray.v*2*np.pi
        # if ray.v > 0:
        #     current_filter_y[ID,index,:] = current_filter_y[ID,index,:] + 2*np.pi*ray.psi
        # else:
        #     current_filter_y[ID,index,:] = current_filter_y[ID,index,:] - 2*np.pi*ray.psi

