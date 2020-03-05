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
Script to test Euler solution with a NACA 4412 straight wing
'''

xdisc=20
airfoil, extra, intra=read_afl('n4412', ext_append=True, remove_TE_gap=True, disc=xdisc, extra_intra=True)
b=2.0
c=0.3
Uinf=0.2
ys=np.linspace(b/2, -b/2, 20)
#ys=np.sin(np.linspace(-pi/2, pi/2, 30))*b/2
elip_factor=np.sqrt((b**2/4-0.9*ys**2)*4/b**2)
#elip_factor=np.array([1.0]*len(ys))
totlist=[]
conlist=[]
n=0
for i in range(np.size(airfoil, 0)):
    totlist+=[[]]
    for j in range(len(ys)):
        totlist[-1]+=[np.array([airfoil[i, 0]*elip_factor[j]*c, ys[j], airfoil[i, 1]*elip_factor[j]*c])]
for i in range(len(ys)-1):
    conlist+=[[i, i+(len(ys)-1)*(np.size(airfoil, 0)-2)]]
sld=Solid(sldlist=[totlist], wraparounds=[[1]])
sld.end_preprocess()
sld.genwakepanels(wakecombs=conlist, a=radians(5.0))
sld.plotgeometry(ylim=[-b/2, b/2], xlim=[-b/2, b/2], zlim=[-b/2, b/2])
t=tm.time()
sld.genvbar(Uinf, a=radians(5.0))
sld.gennvv()
print('Pre-processing and geometry gen: '+str(tm.time()-t))
t=tm.time()
sld.genaicm()
print('AIC matrix generation: '+str(tm.time()-t))
t=tm.time()
sld.solve(damper=0.001)
sld.calcpress(Uinf=Uinf)
print('Solution and post-processing: '+str(tm.time()-t))
sld.plotgeometry(ylim=[-b/2, b/2], xlim=[-b/2, b/2], zlim=[-b/2, b/2])
#sld.plotgeometry()

xlocs1=[]
ylocs1=[]
cps1=[]
xlocs2=[]
ylocs2=[]
cps2=[]
for i in range(sld.npanels):
    if sld.panels[i].nvector[2]>0:
        xlocs1+=[sld.panels[i].colpoint[0]]
        ylocs1+=[sld.panels[i].colpoint[1]]
        cps1+=[sld.Cps[i]]
    else:
        xlocs2+=[sld.panels[i].colpoint[0]]
        ylocs2+=[sld.panels[i].colpoint[1]]
        cps2+=[sld.Cps[i]]

fig=plt.figure()
ax=plt.axes(projection='3d')

ax.scatter3D(xlocs1, ylocs1, cps1, 'red')
ax.scatter3D(xlocs2, ylocs2, cps2, 'blue')
plt.show()

xlocs1=[]
ylocs1=[]
g1=[]
xlocs2=[]
ylocs2=[]
g2=[]
for i in range(sld.npanels):
    if sld.panels[i].nvector[2]>0:
        xlocs1+=[sld.panels[i].colpoint[0]]
        ylocs1+=[sld.panels[i].colpoint[1]]
        g1+=[sld.solution[i]]
    else:
        xlocs2+=[sld.panels[i].colpoint[0]]
        ylocs2+=[sld.panels[i].colpoint[1]]
        g2+=[sld.solution[i]]

fig=plt.figure()
ax=plt.axes(projection='3d')

ax.scatter3D(xlocs1, ylocs1, g1, 'red')
ax.scatter3D(xlocs2, ylocs2, g2, 'blue')
plt.show()

sld.plotnormals(ylim=[-b/2, b/2], xlim=[-b/2, b/2], zlim=[-b/2, b/2], factor=0.01)