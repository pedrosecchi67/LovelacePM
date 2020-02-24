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

class defsect:
    def __init__(self, file='', header_lines=0, coords=np.array([])):
        if len(file)==0:
            self.base_coords=coords
        else:
            infile=open(file, 'r')
            intext=infile.read()
            linelist=intext.split('\n')
            self.base_coords=np.zeros((len(linelist)-header_lines, 2))
            for i in range(header_lines, len(linelist)):
                fllist=linelist[i].split()
                self.base_coords[i, 0]=float(fllist[0])
                self.base_coords[i, 1]=float(fllist[1])
        
        ths=np.zeros(np.size(self.base_coords, 0))
        Rs=np.zeros(np.size(self.base_coords, 0))
        for i in range(len(ths)):
            ths[i]=atan2(self.base_coords[i, 0], self.base_coords[i, 1])
            Rs[i]=lg.norm(self.base_coords[i, :])
        indord=np.argsort(ths)
        ths=ths[indord]
        Rs=Rs[indord]
        
        self.polar_rule=sinterp.CubicSpline(ths, Rs)
    def __call__(self, center=np.array([0.0, 0.0, 0.0]), R=1.0, thspacing=np.linspace(0.0, 2*pi, 20)):
        Rs=self.polar_rule(thspacing)*R
        xy=np.zeros((len(thspacing), 2))
        xy[:, 0]=np.sin(thspacing)*Rs
        xy[:, 1]=np.cos(thspacing)*Rs
        pts=np.zeros((len(thspacing), 3))
        pts[:, 1]=xy[:, 0]
        pts[:, 2]=xy[:, 1]
        for i in range(len(thspacing)):
            pts[i, :]+=center
        return pts

ths=np.linspace(0.0, 2*pi, 20)
circdefsect=defsect(coords=np.vstack((np.sin(ths), np.cos(ths))).T)

def lin_fus_surf(def1=circdefsect, def2=circdefsect, R1=1.0, R2=1.0, c1=np.array([1.0, 0.0, 0.0]), c2=np.array([0.0, 0.0, 0.0]), \
    thspacing1=np.linspace(0.0, 2*pi, 20), thspacing2=np.linspace(0.0, 2*pi, 20), xdisc=10):
    coords1=def1(center=c1, R=R1, thspacing=thspacing1) #thspacing1 and thspacing2 must have the same length
    coords2=def2(center=c2, R=R2, thspacing=thspacing2)
    ptmat=[]
    for i in range(xdisc+1):
        ptmat+=[[]]
        for j in range(len(thspacing1)):
            ptmat[-1]+=[((coords1[j, :]*i)/xdisc+(coords2[j, :]*(xdisc-i))/xdisc)]
    return ptmat

class bodyquad: #non-abuted, wingless, x-wrapped, interpolated fuselage quadrant.
    def __init__(self, sld):
        self.sld=sld
    def patchcompose(self, sects=[circdefsect], xspacing=np.linspace(0.0, 1.0, 10), strategy='spline', Rs=np.ones(10), thdisc=10, \
        centers=[]):
        if len(sects)==0:
            sects=[circdefsect]*len(Rs)
        elif len(sects)==1:
            trimlist(len(Rs), sects)