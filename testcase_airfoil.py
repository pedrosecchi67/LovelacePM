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

xdisc=20
airfoil, extra, intra=read_afl('n4412', ext_append=True, remove_TE_gap=True, disc=xdisc, extra_intra=True)
airfoil[:, 0]-=0.25
extra[:, 0]-=0.25
intra[:, 0]-=0.25
b=2.0
c=0.3
Uinf=0.01
ys=np.linspace(-b/2, b/2, 20)
#ys=np.sin(np.linspace(-pi/2, pi/2, 30))*b/2
#elip_factor=np.sqrt((b**2/4-ys**2)*4/b**2)
elip_factor=np.array([1.0]*len(ys))
totlist=[]
conlist=[]
n=0
for i in range(np.size(airfoil, 0)):
    totlist+=[[]]
    for j in range(len(ys)):
        totlist[-1]+=[np.array([airfoil[i, 0]*elip_factor[j]*c, ys[j], airfoil[i, 1]*elip_factor[j]*c])]
for i in range(len(ys)-1):
    conlist+=[[i, i+(len(ys)-1)*(np.size(airfoil, 0)-2)]]
sld=Solid([totlist], wraparounds=[[1]])
sld.end_preprocess()
sld.genwakepanels(conlist, a=radians(5.0))
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
sld.calcforces()
print('Forces at 10 m/s:')
print('F: '+str(66.25*sld.SCFres))
print('M: '+str(66.25*sld.SCMres))
print('C: ')
print(66.25*sld.SCFres/(b*c))
print(66.25*sld.SCMres/(b*c**2))

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

'''xlocs1=[]
ylocs1=[]
g1=[]
xlocs2=[]
ylocs2=[]
g2=[]
for i in range(sld.npanels):
    if sld.panels[i].nvector[2]>0:
        xlocs1+=[sld.panels[i].colpoint[0]]
        ylocs1+=[sld.panels[i].colpoint[1]]
        g1+=[sld.dbg[i]]
    else:
        xlocs2+=[sld.panels[i].colpoint[0]]
        ylocs2+=[sld.panels[i].colpoint[1]]
        g2+=[sld.dbg[i]]

fig=plt.figure()
ax=plt.axes(projection='3d')

ax.scatter3D(xlocs1, ylocs1, g1, 'red')
ax.scatter3D(xlocs2, ylocs2, g2, 'blue')
plt.show()'''