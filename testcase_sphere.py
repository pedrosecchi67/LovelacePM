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
c=0.01
thetas=np.linspace(-pi, pi, 20)
phis=np.linspace(-pi/2, pi/2, 20)
R=1.0
coordlist=[]
conlist=[]
n=0
thlast=0
for th in thetas:
    coordlist+=[[]]
    for phi in phis:
        coordlist[-1]+=[np.array([sin(phi), cos(phi)*cos(th), cos(phi)*sin(th)])*R]
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
print('Generating AICs: '+str(tm.time()-t))
t=tm.time()
sld.solve(damper=c)
print('Solving and post-processing: '+str(tm.time()-t))

sld.plotgeometry()
plt.scatter([p.colpoint[0] for p in sld.panels], \
    [lg.norm(sld.delphi[i, :]+sld.vbar[i, :])/Uinf \
    for i in range(len(sld.panels))])
#plt.scatter([p.colpoint[0] for p in sld.panels], \
#    [sld.solution[i] \
#    for i in range(len(sld.panels))])
plt.show()

'''n=0
for p in sld.panels:
    n+=1
    print(n)
    print('lines: '+str(p.lines))
    print(sld.lines[np.array(p.lines), :, :])
    print('nvector: '+str(p.nvector))
    print('==========')
    for l in p.lines:
        print('line: '+str(sld.lines[l, :, 1]-sld.lines[l, :, 0]))
        print(p.nvector@(sld.lines[l, :, 1]-sld.lines[l, :, 0]))'''