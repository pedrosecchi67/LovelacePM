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

import pytoolkit
from utils import *

from paneller import *
from body import *
from wing import *

'''ROLLPOINT() DIRECTION: IF P2 IS THE RIGHTMOST POINT, THE FLAP IN QUESTION WILL ROTATE DOWNWARDS (POSITIVE DEFLECTION)
WHEN FUNCTION IS SUMMONED ON THE FLAP'S POINTS'''

def z_rotation_matrix(th):
    return np.array([[cos(th), sin(th), 0.0], [-sin(th), cos(th), 0.0], [0.0, 0.0, 1.0]])

class control_DOF: #control object to be summoned directly by the user
    def __init__(self):
        #ought to be embedded into aircraft class
        self.state=0.0

class control_axis: #definition of control axis with functions to rotate points around it
    #ought to be programmed to be into wing section class
    #to de defined and set before patchcomposing
    def __init__(self, p0=np.array([0.0, 0.0, 0.0]), p1=np.array([0.0, 1.0, 0.0])):
        Mtouni=np.zeros((3, 3))
        Mtouni[:, 2]=p1-p0
        Mtouni[:, 2]/=lg.norm(Mtouni[:, 2])
        Mtouni[:, 0]=np.array([1.0, 0.0, 0.0])
        Mtouni[:, 0]-=Mtouni[:, 2]*(Mtouni[:, 2]@Mtouni[:, 0])
        Mtouni[:, 0]/=lg.norm(Mtouni[:, 0])
        Mtouni[:, 1]=np.cross(Mtouni[:, 2], Mtouni[:, 0])
        self.p0=p0
        self.p1=p1
        self.u=p1-p0
        self.u/=lg.norm(self.u)
        self.control_rot_func=lambda pt, th: Mtouni@(z_rotation_matrix(-th)@(Mtouni.T@(pt-p0)))+p0
    def panel_hinge_moment(self, pt, F): #compute hinge moment for a force F applied at point pt
        return self.u@np.cross(pt-self.p0, F)

class control: #main control class, to be summoned attached to a control_DOF instance ran inside aircraft class instance
    def __init__(self, DOF=None, p0=np.array([0.0, 0.0, 0.0]), p1=np.array([0.0, 1.0, 0.0]), multiplier=1.0):
        self.DOF=DOF
        self.axis=self.axis=control_axis(p0=p0, p1=p1)
        self.multiplier=multiplier
        self.paninds=[]
    def panlist(self, paninds, sld): #return indexes of panels belonging to the control surface
        return [pan for pan in paninds if np.cross(sld.panels[pan].colpoint[0]-self.axis.p0, self.axis.u)@np.cross(np.array([1.0, 0.0, 0.0]), self.axis.u)>0.0]
    def hinge_moment(self, panlist, sld):
        H=0.0
        S=0.0
        for pan in self.panlist(paninds=panlist, sld=sld):
            H+=self.axis.panel_hinge_moment(sld.panels[pan].colpoint, sld.forces[pan])
            S+=sld.panels[pan].S
        return H