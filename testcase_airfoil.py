import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time as tm

import toolkit

from paneller import *
from utils import *

'''
Script to test the Euler solution with a NACA 0012 straight wing
'''

afl=read_afl('n0012', ext_append=True)
b=1.0
c=0.3
Uinf=0.01
afl=afl*c
ys=np.linspace(-b/2, b/2, 30)
coords=np.zeros((3, 4))
coordlist=[]
conlist=[]
n=0
last=0
for j in range(len(ys)-1):
    for i in range(np.size(afl, 0)-1):
        coords=np.zeros((3, 4))
        coords[:, 0]=np.array([afl[i, 0], ys[j], afl[i, 1]])
        coords[:, 1]=np.array([afl[i, 0], ys[j+1], afl[i, 1]])
        coords[:, 2]=np.array([afl[i+1, 0], ys[j+1], afl[i+1, 1]])
        coords[:, 3]=np.array([afl[i+1, 0], ys[j], afl[i+1, 1]])
        coordlist+=[coords]
        if i==0:
            last=n
        if i==np.size(afl, 0)-2:
            conlist+=[[n, last]]
        n+=1
sld=Solid(coordlist)
sld.genwakepanels(conlist)
sld.genvbar(Uinf)
sld.gennvv()
sld.genaicm()
sld.solve(damper=0.001)
sld.calcpress(Uinf=Uinf)
sld.plotpress(factor=0.001)