import numpy as np
from geom import *
from util import *

class Ray:
    def __init__(self, box):

        xH = box[0]
        xL = box[3]
        self.x = (xH-xL)*uniform_dis()+xL
        yH = box[1]
        yL = box[4]
        self.y = (yH-yL)*uniform_dis()+yL
        zH = box[2]
        zL = box[5]
        self.z = (zH-zL)*uniform_dis()+zL

        mu = mu_dis()
        phi = circle_dis()
        self.u = np.sqrt(1-mu**2)*np.cos(phi)
        self.v = np.sqrt(1-mu**2)*np.sin(phi)
        self.w = mu

        self.out = 0
        
    def set_psi(self,init_psi):
        self.psi = np.asarray(init_psi)

    def set_pos(self,position):
        self.x = position[0]
        self.y = position[1]
        self.z = position[2]
    def set_dir(self,direction):
        self.u = direction[0]
        self.v = direction[1]
        self.w = direction[2]
    def get_pos(self):
    	return np.asarray([self.x, self.y, self.z])
    def get_dir(self):
    	return np.asarray([self.u, self.v, self.w])
    def to_next_surf(self,d):
        self.x += d*self.u
        self.y += d*self.v
        self.z += d*self.w
    def reflect(self,plane):
       	A = plane.A
        B = plane.B
        C = plane.C
        D = plane.D
        self.u -= 2*A*(A*self.u+B*self.v+C*self.w)/(A**2+B**2+C**2)
        self.v -= 2*B*(A*self.u+B*self.v+C*self.w)/(A**2+B**2+C**2)
        self.w -= 2*C*(A*self.u+B*self.v+C*self.w)/(A**2+B**2+C**2)  