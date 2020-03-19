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
        self.control_rot_func=lambda pt, th: Mtouni@(z_rotation_matrix(-th)@(Mtouni.T@(pt-p0)))+p0

class control: #main control class, to be summoned attached to a control_DOF instance ran inside aircraft class instance
    def __init__(self, DOF=None, p0=np.array([0.0, 0.0, 0.0]), p1=np.array([0.0, 1.0, 0.0]), multiplier=1.0):
        self.DOF=DOF
        self.axis=self.axis=control_axis(p0=p0, p1=p1)
        self.multiplier=multiplier
        self.paninds=[]
    def addpanels(self, panlist): #function for later programming of control derivative and hinge moment computation
        self.paninds+=panlist