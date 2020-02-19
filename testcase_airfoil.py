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
airfoil, extra, intra=read_afl('n4412', ext_append=True, remove_TE_gap=False, disc=xdisc, extra_intra=True)
airfoil[:, 0]-=0.25
extra[:, 0]-=0.25
intra[:, 0]-=0.25
b=2.0
c=0.3
Uinf=0.01
#ys=np.linspace(-b/2, b/2, 20)
ys=np.sin(np.linspace(-pi/2, pi/2, 40))*b/2
#elip_factor=np.sqrt((b**2/4-ys**2)*4/b**2)
elip_factor=np.array([1.0]*len(ys))
coords=np.zeros((3, 4))
coordlist=[]
totlist=[]
conlist=[]
indlist=[]
n=0
last=0
ltipmesh=[[]]
rtipmesh=[[]]
for i in range(np.size(extra, 0)-1):
    coords=np.zeros((3, 4))
    coords[:, 0]=np.array([intra[i, 0]*elip_factor[0]*c, -b/2, (intra[i, 1]+extra[i, 1])*elip_factor[0]*c/2])
    coords[:, 1]=np.array([extra[i, 0]*elip_factor[0]*c, -b/2, extra[i, 1]*elip_factor[0]*c])
    coords[:, 2]=np.array([extra[i+1, 0]*elip_factor[0]*c, -b/2, extra[i+1, 1]*elip_factor[0]*c])
    coords[:, 3]=np.array([intra[i+1, 0]*elip_factor[0]*c, -b/2, (intra[i+1, 1]+extra[i+1, 1])*elip_factor[0]*c/2])
    ltipmesh[-1]+=[coords]
    ltipmesh+=[[]]
    coords=np.zeros((3, 4))
    coords[:, 0]=np.array([intra[i, 0]*elip_factor[0]*c, -b/2, intra[i, 1]*elip_factor[0]*c])
    coords[:, 1]=np.array([extra[i, 0]*elip_factor[0]*c, -b/2, (intra[i, 1]+extra[i, 1])*elip_factor[0]*c/2])
    coords[:, 2]=np.array([extra[i+1, 0]*elip_factor[0]*c, -b/2, (intra[i+1, 1]+extra[i+1, 1])*elip_factor[0]*c/2])
    coords[:, 3]=np.array([intra[i+1, 0]*elip_factor[0]*c, -b/2, intra[i+1, 1]*elip_factor[0]*c])
    ltipmesh[-1]+=[coords]
    coords=np.zeros((3, 4))
    coords[:, 0]=np.array([intra[i, 0]*elip_factor[0]*c, b/2, (intra[i, 1]+extra[i, 1])*elip_factor[0]*c/2])
    coords[:, 1]=np.array([extra[i, 0]*elip_factor[0]*c, b/2, extra[i, 1]*elip_factor[0]*c])
    coords[:, 2]=np.array([extra[i+1, 0]*elip_factor[0]*c, b/2, extra[i+1, 1]*elip_factor[0]*c])
    coords[:, 3]=np.array([intra[i+1, 0]*elip_factor[0]*c, b/2, (intra[i+1, 1]+extra[i+1, 1])*elip_factor[0]*c/2])
    rtipmesh[-1]+=[coords]
    rtipmesh+=[[]]
    coords=np.zeros((3, 4))
    coords[:, 0]=np.array([intra[i, 0]*elip_factor[0]*c, b/2, intra[i, 1]*elip_factor[0]*c])
    coords[:, 1]=np.array([extra[i, 0]*elip_factor[0]*c, b/2, (intra[i, 1]+extra[i, 1])*elip_factor[0]*c/2])
    coords[:, 2]=np.array([extra[i+1, 0]*elip_factor[0]*c, b/2, (intra[i+1, 1]+extra[i+1, 1])*elip_factor[0]*c/2])
    coords[:, 3]=np.array([intra[i+1, 0]*elip_factor[0]*c, b/2, intra[i+1, 1]*elip_factor[0]*c])
    rtipmesh[-1]+=[coords]
for j in range(len(ys)-1):
    afl1=airfoil*c*elip_factor[j]
    afl2=airfoil*c*elip_factor[j+1]
    coordlist=[]
    for i in range(np.size(airfoil, 0)-1):
        coords=np.zeros((3, 4))
        coords[:, 0]=np.array([afl1[i, 0], ys[j], afl1[i, 1]])
        coords[:, 1]=np.array([afl2[i, 0], ys[j+1], afl2[i, 1]])
        coords[:, 2]=np.array([afl2[i+1, 0], ys[j+1], afl2[i+1, 1]])
        coords[:, 3]=np.array([afl1[i+1, 0], ys[j], afl1[i+1, 1]])
        coordlist+=[coords]
        if i==0:
            last=n
            #conlist+=[[n, -1]]
        if i==np.size(airfoil, 0)-2:
            conlist+=[[n, last]]
            #conlist+=[[n, -1]]
        n+=1
    totlist+=[coordlist]
sld=Solid([ltipmesh, totlist, rtipmesh])
sld.genwakepanels(conlist, a=radians(5.0))
sld.plotgeometry(wake=True, ylim=[-b/2, b/2], xlim=[-b/2, b/2], zlim=[-b/2, b/2])
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
#sld.plotgeometry(ylim=[-b/2, b/2], xlim=[-b/2, b/2], zlim=[-b/2, b/2])
sld.plotgeometry()

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
        g1+=[sld.dbg[i]]
    else:
        xlocs2+=[sld.panels[i].colpoint[0]]
        ylocs2+=[sld.panels[i].colpoint[1]]
        g2+=[sld.dbg[i]]

fig=plt.figure()
ax=plt.axes(projection='3d')

ax.scatter3D(xlocs1, ylocs1, g1, 'red')
ax.scatter3D(xlocs2, ylocs2, g2, 'blue')
plt.show()