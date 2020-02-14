import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
import scipy.sparse.linalg as splg
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time as tm

import toolkit

'''Panel order when generating wake:
Upper surface: (or airfoil)
     |           |
-----p3----------p2-----
     |           |
     |           |
     |           |
-----p0----------p1-----
     w           w
     w           w
     w           w
lower surface:
     |           |
-----p0----------p1-----
     |           |
     |           |
     |           |
-----p3----------p2-----
     w           w
     w           w
     w           w
     '''

class Panel: #Panel data type
    def __init__(self, coords): #constructor to define central collocation point and normal vector
        self.coords=coords.astype(dtype='double')
        self.colpoint=np.mean(coords, axis=1)
        self.nvector=(np.cross(coords[:, 1]-coords[:, 0], coords[:, 3]-coords[:, 0])+np.cross(coords[:, 3]-coords[:, 2], coords[:, 1]-coords[:, 2]))/2
        self.nvector/=lg.norm(self.nvector)

class Solid:
    #WARNING: USE LEN(SELF.PANELS) ONLY WHEN REFERRING TO ALL (WAKE+GEOMETRY) PANELS.SELF.NPANELS AND SELF.NWAKE ARE THERE FOR OTHER PURPOSES
    def __init__(self, coordlist): #initialize solid geometry with non-wake panels specified in list. INPUT WITH 3X4 ARRAY REFERRING TO FOUR POINT COORDS
        self.panels=[]
        self.npanels=0
        self.nwake=0
        self.addto=[]
        for arr in coordlist:
            self.panels+=[Panel(arr)]
            self.npanels+=1
        self.delphi=np.array([[0.0, 0.0, 0.0]]*len(self.panels), dtype='double')
        self.vbar=np.array([[0.0, 0.0, 0.0]]*len(self.panels), dtype='double')
        self.nvv=np.array([0.0]*len(self.panels), dtype='double')
        self.solution=np.array([0.0]*len(self.panels), dtype='double')
        self.Cps=np.array([0.0]*len(self.panels), dtype='double')
    def genwakepanels(self, wakecombs, offset=1000.0, a=0.0, b=0.0): #generate wake panels based in list of TE panel combinations
        #wakecombs: list of lists, first element in sublist is upper surface panel
        coords=np.array([[0.0]*4]*3)
        for i in range(len(wakecombs)):
            self.addto+=[[wakecombs[i][0]+1, wakecombs[i][1]+1]]
            if wakecombs[i][1]!=-1:
                coords[:, 0]=(self.panels[wakecombs[i][0]].coords[:, 0]+self.panels[wakecombs[i][1]].coords[:, 3])/2
                coords[:, 1]=(self.panels[wakecombs[i][0]].coords[:, 1]+self.panels[wakecombs[i][1]].coords[:, 2])/2
            else:
                coords[:, 0]=self.panels[wakecombs[i][0]].coords[:, 0]
                coords[:, 1]=self.panels[wakecombs[i][0]].coords[:, 1]
            coords[:, 2]=coords[:, 1]+np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)])*offset
            coords[:, 3]=coords[:, 0]+np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)])*offset
            self.panels+=[Panel(coords)]
            self.nwake+=1
    def genaicm(self): #call FORTRAN backend to generate AIC matrix
        coordmat=np.array([[[0.0, 0.0, 0.0, 0.0]]*3]*len(self.panels), dtype='double')
        #nvectmat=np.array([[0.0, 0.0, 0.0]]*len(self.panels))
        colmat=np.array([[0.0, 0.0, 0.0]]*self.npanels, dtype='double')
        for i in range(len(self.panels)):
            coordmat[i, :, :]=self.panels[i].coords
            #nvectmat[i, :]=self.panels[i].nvector
            if i<self.npanels:
                colmat[i, :]=self.panels[i].colpoint
        #self.aicm3=np.array([[[0.0]*len(self.panels)]*len(self.panels)]*3)
        if self.nwake!=0:
            self.aicm3=toolkit.gen_aicm(coordmat, colmat, np.array(self.addto))
        else:
            self.aicm3=toolkit.gen_aicm_nowake(coordmat, colmat)
        self.aicm=np.array([[0.0]*self.npanels]*self.npanels, dtype='double')
        for i in range(self.npanels):
            for j in range(self.npanels):
                self.aicm[i, j]=self.aicm3[:, i, j]@self.panels[i].nvector
    def genvbar(self, Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0): #generate freestream velocity vector. Angular velocities are raw (rad/s), not normalized by Uinf or dimensions
        for i in range(self.npanels):
            self.vbar[i, :]=Uinf*np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)], dtype='double')+np.cross(np.array([-p, q, -r], dtype='double'), self.panels[i].colpoint)
    def gennvv(self): #Generate vector containing normal velocities
        for i in range(self.npanels):
            self.nvv[i]=self.panels[i].nvector@self.vbar[i, :]
    def calcpress(self, Uinf=1.0): #calculate local pressure coefficients
        for i in range(self.npanels):
            self.Cps[i]=(Uinf**2-lg.norm(self.delphi[i, :]+self.vbar[i, :])**2)/Uinf**2
    def solve(self, damper=0.0): #generate Euler solution. Inverts AIC matrix with Tikhonov regularization if "damper" is set to non-zero variable.
        #self.iaicm=toolkit.aicm_inversion(self.aicm, damper)
        if damper!=0.0:
            self.iaicm=slg.inv(self.aicm.T@self.aicm+damper*np.eye(self.npanels, self.npanels))@self.aicm.T
        else:
            self.iaicm=slg.inv(self.aicm)
        #self.solution=slg.solve(self.aicm, -self.nvv)
        self.solution=-self.iaicm@self.nvv
        for i in range(3):
            self.delphi[:, i]=self.aicm3[i, :, :]@self.solution #compute loval velocity due to panel influence
    def plotgeometry(self, wake=False): #plot geometry and local velocity vectors, either with or without wake panels
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        order=np.array([0, 1, 2, 3, 0])
        for i in range(len(self.panels)):
            if i<self.npanels or wake:
                ax.plot3D(self.panels[i].coords[0, order], self.panels[i].coords[1, order], self.panels[i].coords[2, order], 'gray')
                #ax.scatter3D(self.panels[i].colpoint[0], self.panels[i].colpoint[1], self.panels[i].colpoint[2], 'yellow')
            if i<self.npanels:
                ax.quiver(self.panels[i].colpoint[0], self.panels[i].colpoint[1], self.panels[i].colpoint[2], self.vbar[i, 0]+self.delphi[i, 0], \
                    self.vbar[i, 1]+self.delphi[i, 1], self.vbar[i, 2]+self.delphi[i, 2])
                #ax.quiver(self.panels[i].colpoint[0], self.panels[i].colpoint[1], self.panels[i].colpoint[2], self.delphi[i, 0], \
                #    self.delphi[i, 1], self.delphi[i, 2])
                #ax.quiver(self.panels[i].colpoint[0], self.panels[i].colpoint[1], self.panels[i].colpoint[2], \
                #    self.panels[i].nvector[0]*self.solution[i], \
                #    self.panels[i].nvector[1]*self.solution[i], self.panels[i].nvector[2]*self.solution[i])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    def plotpress(self, factor=1.0): #plot local Cps (scaled by a floating point factor so as to ease visualization) as normal vectors
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        order=np.array([0, 1, 2, 3, 0])
        for i in range(self.npanels):
            ax.plot3D(self.panels[i].coords[0, order], self.panels[i].coords[1, order], self.panels[i].coords[2, order], 'gray')
            ax.quiver(self.panels[i].colpoint[0], self.panels[i].colpoint[1], self.panels[i].colpoint[2], self.Cps[i]*self.panels[i].nvector[0]*factor, \
                self.Cps[i]*self.panels[i].nvector[1]*factor, self.Cps[i]*self.panels[i].nvector[2]*factor)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

'''sld=Solid([np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]).T])
#sld.plotgeometry()
sld.genaicm()
print(sld.aicm3[:, 0, 0])'''