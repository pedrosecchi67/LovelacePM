import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time as tm

from paneller import *
from utils import *
from wing import *
from body import *
from aircraft import *
from aerodynamic_output import *

Uinf=1.0
t=tm.time()
c=0.0
thetas=np.linspace(0.0, 2*pi, 100)
ys=np.linspace(-1000.0, 1000.0, 5)
R=1.0
coordlist1=[]

for th in thetas:
    coordlist1+=[[]]
    for i in range(len(ys)):
        coordlist1[-1]+=[np.array([sin(th), ys[i], cos(th)])*R]
        '''
for i in range(len(coordlist1)):
    coordlist1[i]=[np.array([0.0, ys[0], 0.0])]+coordlist1[i]+[np.array([0.0, ys[-1], 0.0])]'''
    

sld=Solid()
horzlines, vertlines, paninds, sldpts=sld.addpatch(coordlist1, wraps=[1])

sld.end_preprocess()

sld.genvbar(Uinf=Uinf)
sld.gennvv()
print('Solid generation and pre-processing: '+str(tm.time()-t))
t=tm.time()
sld.genaicm()
sld.gen_selfinf_mat()
print('Generating AICs: '+str(tm.time()-t))
t=tm.time()
sld.solve(damper=c)
print('Solving and post-processing: '+str(tm.time()-t))

sld.plotgeometry()
plt.scatter([p.colpoint[0] for p in sld.panels], \
    [lg.norm(sld.delphi[i, :]+sld.vbar[i, :])/Uinf \
    for i in range(len(sld.panels))])
plt.xlabel('$x$')
plt.ylabel('$v$')
plt.show()
