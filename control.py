import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
import scipy.sparse.linalg as splg
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.sparse as sps
import scipy.interpolate as sinterp
import time as tm

import toolkit
from utils import *

from paneller import *
from body import *
from wing import *

'''ROLLPOINT() DIRECTION: IF P2 IS THE RIGHTMOST POINT, THE FLAP IN QUESTION WILL ROTATE DOWNWARDS (POSITIVE DEFLECTION)
WHEN FUNCTION IS SUMMONED ON THE FLAP'S POINTS'''

class control: #control object to be summoned directly by the user
    def __init__(self, multiplier=1.0):
        self.multiplier=1.0
        self.state=0.0

class control_axis: #definition of control axis with functions to rotate points around it
    def __init__(self, sld, p1=np.array([0.0, 0.0, 0.0]), p2=np.array([0.0, 1.0, 0.0]), control_DOF=None):
        self.sld=sld
        self.p0=p1
        self.p1=p2
        u=p2-p1
        nu=lg.norm(u)
        self.Mtosys=np.zeros((3, 3))
        self.Mtosys[2, :]=u/nu
        self.Mtosys[0, :]=np.array([1.0, 0.0, 0.0]) #LIMITATION: dont add a control axis in the x direction or too close to it
        self.Mtosys[0, :]-=self.Mtosys[2, :]*(self.Mtosys[0, :]@self.Mtosys[2, :])
        self.Mtosys[0, :]/=lg.norm(self.Mtosys[0, :])
        self.Mtosys[1, :]=np.cross(self.Mtosys[2, :], self.Mtosys[0, :])
        self.Mtouni=self.Mtosys.T
        self.Rmat=lambda th: np.array([[cos(th), sin(th), 0.0], [-sin(th), cos(th), 0.0], [0.0, 0.0, 1.0]])
        if control_DOF==None:
            print('WARNING: control object absent or set as None')
            self.delta=0.0
        else:
            self.delta=control_DOF.state*control_DOF.multiplier
        self.control_DOF=control_DOF
    def update_control(self, linelist=[], paninds=[]): #linelist contains line indexes in self.sld
        #paninds indicates the panel indexes representing the panels which will have to be updated
        d_delta=self.control_DOF.state*self.control_DOF.multiplier-self.delta
        self.delta+=d_delta
        for l in linelist:
            if self.belongs(self.sld.linelist[l, :, 0]):
                self.sld.linelist[l, :, 0]=self.rollpoint(self.sld.linelist[l, :, 0], d_delta)
            if self.belongs(self.sld.linelist[l, :, 1]):
                self.sld.linelist[l, :, 1]=self.rollpoint(self.sld.linelist[l, :, 1], d_delta)
        #update panel info for deformed panels
        self.sld.end_preprocess(paninds=paninds)
    def rollpoint(self, p, th):
        return self.Mtouni@(self.Rmat(th)@(self.Mtosys@(p-self.p0)))+self.p0
    def belongs(self, p):
        return (p[0]>(np.interp(p[1], np.array([self.p0[1], self.p1[1]]), np.array([self.p0[0], self.p1[0]]))))