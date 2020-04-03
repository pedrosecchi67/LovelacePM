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

'''
Script to test the Euler solution modules with a generic sphere geometry
'''

Uinf=1.0
t=tm.time()
c=0.0001
thetas=np.linspace(0.0, 2*pi, 20)
phis=np.linspace(-pi/2, pi/2, 100)
R=1.0
coordlist=[]
conlist=[]
n=0
thlast=0
exc=1.0
for th in thetas:
    coordlist+=[[]]
    for phi in phis:
        coordlist[-1]+=[np.array([sin(phi)*exc, cos(phi)*cos(th), cos(phi)*sin(th)])*R]
for i in range(len(phis)-1):
    for j in range(len(thetas)-1):
        if j==0:
            thlast=n
        if j==len(thetas)-2:
            conlist+=[[thlast, n]]
        n+=1
sld=Solid(sldlist=[coordlist], wraparounds=[[1]])
sld.end_preprocess()
#sld.genwakepanels(wakecombs=conlist)
sld.genvbar(Uinf=Uinf)
sld.gennvv()
print('Solid generation and pre-processing: '+str(tm.time()-t))
t=tm.time()
sld.genaicm()
sld.gen_selfinf_mat()

#checking aicm with srivastava's paper
'''print(sld.aicm3_line)
srivastava_order=np.array([3, 2, 0, 1])
equivalences=np.array([6, 7, 4, 5])
ptinds=np.arange(0, 4, dtype='int')
aicm=(sld.aicm[0:4, equivalences]+sld.aicm[0:4, ptinds])
for i in range(sld.npanels):
    print(i)
    lines, dumpme=np.nonzero(sld.panline_matrix[:, i])
    for l in lines:
        print(sld.lines[l, :, :])
    print(sld.panels[i].colpoint)
aicm=aicm[srivastava_order, :]
aicm=aicm[:, srivastava_order]
print(aicm)'''

print('Generating AICs: '+str(tm.time()-t))
t=tm.time()

sld.solve(damper=c)
print('Solving and post-processing: '+str(tm.time()-t))

'''for l in range(sld.nlines):
    for p in range(sld.npanels):
        print('Line')
        print(sld.lines[l, :, :])
        print('to panel')
        print(sld.panels[p].colpoint)
        print('n')
        print(sld.panels[p].nvector)
        print(sld.aicm3_line[:, p, l])
        print(sld.panels[p].nvector@sld.aicm3_line[:, p, l])'''

sld.plotgeometry()
plt.scatter([p.colpoint[0] for p in sld.panels], \
    [lg.norm(sld.delphi[i, :]+sld.vbar[i, :])/Uinf \
    for i in range(len(sld.panels))])
plt.xlabel('$x$')
plt.ylabel('$v$')
plt.show()