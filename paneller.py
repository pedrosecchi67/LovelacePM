import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
import scipy.sparse.linalg as splg
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.sparse as sps
import time as tm

import toolkit
from utils import *

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
    def __init__(self, lines):
        self.lines=lines
        self.wakelines=[]

class Solid:
    def __init__(self, sldlist=[], wraparounds=[]): #solid data type
        self.panels=[]
        self.lines=[]
        self.addto=[]
        self.npanels=0
        self.nlines=0
        self.nwake=0 #number of wake panels
        self.solavailable=False
        if len(wraparounds)==0:
            wraparounds=[[]]*len(sldlist)
        horzlines=[]
        vertlines=[]
        sldcnt=0
        for sld in sldlist: #define first patches (without their interconnections)
            wraps=wraparounds[sldcnt]
            self.addpatch(sld, wraps=wraps)
            sldcnt+=1
    def addpatch(self, sld, wraps=[], prevlines={}, invlats=[]): #add panel based on point grid (list of lists)
        # returns: arrays of, respectively, horizontal and vertical line indexes; panel indexes list of lists;
        # point grid inserted as input, for external reference from high-end functions
        # lines provided with -2 index in prevlines will be accounted as new lines, and as -1, 
        # as inexistent (dismissable length) lines
        horzlines=[]
        vertlines=[]
        paninds=[]
        for i in range(len(sld)-1):
            horzlines+=[[]]
            vertlines+=[[]]
            for j in range(len(sld[0])-1):
                if i==0 and 'low' in prevlines:
                    if prevlines['low'][j]!=-2:
                        horzlines[-1]+=[prevlines['low'][j]]
                    else:
                        horzlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i][j+1])).T)]
                else:
                    horzlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i][j+1])).T)]
                if j==0 and 'right' in prevlines:
                    if prevlines['right'][i]!=-2:
                        vertlines[-1]+=[prevlines['right'][i]]
                    else:
                        vertlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i+1][j])).T)]
                else:
                    vertlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i+1][j])).T)]
        if 0 in wraps: #if wrapped in x direction (inner list layer), repeat the first lateral vertical lines
            for i in range(len(sld)-1):
                vertlines[i]+=[vertlines[i][0]]
        else:
            if 'left' in prevlines:
                for i in range(len(sld)-1):
                    if prevlines['left'][i]!=-2:
                        vertlines[i]+=[prevlines['left'][i]]
                    else:
                        vertlines[i]+=[self.addline(np.vstack((sld[i][-1], sld[i+1][-1])).T)]
            else:
                for i in range(len(sld)-1):
                    vertlines[i]+=[self.addline(np.vstack((sld[i][-1], sld[i+1][-1])).T)]
        horzlines+=[[]]
        if 1 in wraps: #same for y axis
            for i in range(len(sld[0])-1):
                horzlines[-1]+=[horzlines[0][i]]
        else:
            if 'up' in prevlines:
                for i in range(len(sld[0])-1):
                    if prevlines['up'][i]!=-2:
                        horzlines[-1]+=[prevlines['up'][i]]
                    else:
                        horzlines[-1]+=[self.addline(np.vstack((sld[-1][i], sld[-1][i+1])).T)]
            else:
                for i in range(len(sld[0])-1):
                    horzlines[-1]+=[self.addline(np.vstack((sld[-1][i], sld[-1][i+1])).T)]
        for i in range(len(sld)-1):
            paninds+=[[]]
            for j in range(len(sld[0])-1):
                l=[]
                if j==len(sld[0])-2 and 'left' in invlats:
                    l+=[3]
                if j==0 and 'right' in invlats:
                    l+=[1]
                if i==0 and 'low' in invlats:
                    l+=[0]
                if i==len(sld)-2 and 'up' in invlats:
                    l+=[2]
                paninds[-1]+=[self.addpanel([horzlines[i][j], vertlines[i][j], horzlines[i+1][j], vertlines[i][j+1]], invs=l)]
        return horzlines, vertlines, paninds, sld
    def addline(self, coords, tolerance=0.00005):
        #neglect lines of neglectible size
        if lg.norm(coords[:, 1]-coords[:, 0])>tolerance:
            if np.size(self.lines, 0)==0:
                self.lines=np.array([list(coords)])
            else:
                self.lines=np.append(self.lines, np.array([list(coords)]), axis=0)
            self.nlines+=1
            return self.nlines-1
        else:
            return -1
    def addpanel(self, lines, invs=[]): #line list: includes indexes summed by one, inverted to indicate oposite direction
        #for vortex line segment
        self.npanels+=1
        llist=[]
        n=0
        for l in lines:
            if l!=-1:
                if (n==0 or n==3):
                    llist+=[-1-l]
                else:
                    llist+=[l+1]
                if n in invs:
                    llist[-1]=-llist[-1]
            n+=1
        self.panels+=[Panel(llist)]
        return self.npanels-1
    def addwakepanel(self, refup, refdown, indup=0, indown=2, offsetleft=np.array([1000.0, 0.0, 0.0]), offsetright=np.array([]), \
        tolerance=0.00005): #add line segments correspondent to panel wake, and remove "inup-th"/"indown-th" line segment in the panel's list
        if len(offsetright)==0:
            offsetright=offsetleft
        self.nwake+=1
        if indup<len(self.panels[refup].lines):
            te_linind=self.panels[refup].lines.pop(indup)
        else:
            te_linind=self.panels[refup].wakelines.pop(indup-len(self.panels[refup].lines))
        if te_linind>0:
            #first wake vortex
            p0=self.lines[te_linind-1, :, 0]
            p1=self.lines[te_linind-1, :, 0]+offsetleft
            indsup=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)]
            #farfield wake vortex
            p0=p1
            p1=self.lines[te_linind-1, :, 1]+offsetright
            indsup=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)]
            #second wake vortex
            p0=p1
            p1=self.lines[te_linind-1, :, 1]
            indsup=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)]
        else:
            #first wake vortex
            p0=self.lines[-te_linind-1, :, 1]
            p1=self.lines[-te_linind-1, :, 1]+offsetright
            indsup=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)]
            #farfield wake vortex
            p0=p1
            p1=self.lines[-te_linind-1, :, 0]+offsetleft
            indsup=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)]
            #second wake vortex
            p0=p1
            p1=self.lines[-te_linind-1, :, 0]
            indsup=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)]
        self.panels[refup].wakelines+=indsup
        if refdown!=-1:
            self.panels[refdown].lines.pop(indown)
            self.panels[refdown].wakelines+=[-l for l in indsup]
    def genwakepanels(self, wakecombs=[], wakeinds=[], offset=1000.0, a=0.0, b=0.0): #generate wake panels based in list of TE panel combinations
        #wakecombs: list of lists, first element in sublist is upper surface panel
        #wakeinds: list of lists indicating corresponding wake vortex line segment index to apply kutta condition to
        if len(wakeinds)==0:
            wakeinds=[[0, 2]]
        trimlist(len(wakecombs), wakeinds)
        for i in range(len(wakecombs)):
            self.addto+=[[wakecombs[i][0], wakecombs[i][1]]]
            self.addwakepanel(wakecombs[i][0], wakecombs[i][1], indup=wakeinds[i][0], indown=wakeinds[i][0], \
                offsetleft=np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)])*offset)
        #straight, alpha and beta defined single wake panel as default
    def panel_getcoords(self, p): #get coordinates for panel based on line vectors and lists
        if p.lines[0]<0:
            coords=self.lines[-1-p.lines[0], :, :].T
        else:
            coords=self.lines[p.lines[0]-1, :, :].T
        for i in range(1, len(p.lines)):
            if p.lines[i]<0:
                coords=np.vstack((coords, self.lines[-1-p.lines[i], :, 1]))
            else:
                coords=np.vstack((coords, self.lines[p.lines[i]-1, :, 1]))
        return coords.T
    def line_getcoords(self, ind): #return line coordinates in order presented in panel, based on panel.lines element
        if ind>0:
            return self.lines[ind-1, :, :]
        else:
            return self.lines[-1-ind, :, np.array([1, 0])]
    def line_getvec(self, ind): #return vector linking points in line
        if ind>0:
            return self.lines[ind-1, :, 1]-self.lines[ind-1, :, 0]
        else:
            return self.lines[-ind-1, :, 0]-self.lines[-ind-1, :, 1]
    def line_midpoint(self, ind): #return midpoint of line
        return np.mean(self.lines[abs(ind)-1, :, :], axis=1)
    def nvect_diradjust(self, patchinds, vect): #adjust normal vectors in patch to follow a certain direction
        for patchlist in patchinds:
            for i in patchlist:
                if self.panels[i].nvector@vect<0.0:
                    self.panels[i].nvector*=-1
    def nvect_radadjust(self, patchinds, center, inwards=False): #adjust normal vectors in patch to follow radial convergent or 
        #divergent distribution
        if inwards:
            for patchlist in patchinds:
                for i in patchlist:
                    if self.panels[i].nvector@(self.panels[i].colpoint-center)>0.0:
                        self.panels[i].nvector*=-1
        else:
            for patchlist in patchinds:
                for i in patchlist:
                    if self.panels[i].nvector@(self.panels[i].colpoint-center)<0.0:
                        self.panels[i].nvector*=-1
    def lineadjust(self, patchinds=[]): #adjust line fortran index signals to comply with dextrogyre panel convention
        if len(patchinds)==0:
            patchinds=[list(range(self.npanels))]
        p=np.array([0.0, 0.0, 0.0])
        u=np.array([0.0, 0.0, 0.0])
        for patchlist in patchinds:
            for i in patchlist:
                for lind in range(len(self.panels[i].lines)):
                    p=self.lines[abs(self.panels[i].lines[lind])-1, :, 0]-self.panels[i].colpoint
                    u=self.lines[abs(self.panels[i].lines[lind])-1, :, 1]-self.lines[abs(self.panels[i].lines[lind])-1, :, 0]
                    if (u@np.cross(self.panels[i].nvector, p)*self.panels[i].lines[lind]<0.0):
                        self.panels[i].lines[lind]*=-1
    def end_preprocess(self): #calculate panel areas before they are altered by wake generation. Must be run before it, and terrible
        #consequences may arise if done otherwise
        u=np.array([0.0, 0.0, 0.0])
        v=np.array([0.0, 0.0, 0.0])
        for p in self.panels:
            '''coords=self.panel_getcoords(p)
            u=coords[:, 1]-coords[:, 0]
            v=coords[:, 2]-coords[:, 1]'''
            u=self.line_getvec(p.lines[0])
            v=self.line_getvec(p.lines[1])
            p.nvector=np.cross(u, v)
            p.S=0.0
            if len(p.lines)==3:
                p.S=lg.norm(p.nvector)/2
                p.nvector/=lg.norm(p.nvector)
                p.colpoint=(self.line_midpoint(p.lines[0])+self.line_midpoint(p.lines[1])+self.line_midpoint(p.lines[2]))/3
                #p.colpoint=(coords[:, 0]+coords[:, 1]+coords[:, 2])/3
            elif len(p.lines)==4:
                p.S=(lg.norm(p.nvector)+lg.norm(np.cross(self.line_getvec(p.lines[2]), self.line_getvec(p.lines[3]))))/2
                #p.S=(lg.norm(p.nvector)+lg.norm(np.cross(coords[:, 1]-coords[:, 2], coords[:, 3]-coords[:, 2])))/2
                p.nvector/=lg.norm(p.nvector)
                #p.colpoint=(coords[:, 0]+coords[:, 1]+coords[:, 2]+coords[:, 3])/4
                p.colpoint=(self.line_midpoint(p.lines[0])+self.line_midpoint(p.lines[1])+self.line_midpoint(p.lines[2])+\
                    self.line_midpoint(p.lines[3]))/4
            else:
                p.nvector=np.array([0.0, 0.0, 0.0])
                #p.colpoint=(coords[:, 0]+coords[:, 1])/2
                p.colpoint=(self.line_midpoint(p.lines[0])+self.line_midpoint(p.lines[1]))/2
        self.lineadjust()
        #initialize result vectors
        self.delphi=np.array([[0.0, 0.0, 0.0]]*len(self.panels), dtype='double')
        self.vbar=np.array([[0.0, 0.0, 0.0]]*len(self.panels), dtype='double')
        self.nvv=np.array([0.0]*len(self.panels), dtype='double')
        self.solution=np.array([0.0]*len(self.panels), dtype='double')
        self.solavailable=False
        self.Cps=np.array([0.0]*len(self.panels), dtype='double')
    def gen_panline(self): #generate panel-line correspondence matrix
        self.panline_matrix=np.zeros((self.nlines, self.npanels), dtype='double')
        for i in range(self.npanels):
            lininds=np.array(self.panels[i].lines+self.panels[i].wakelines)
            lintemp=lininds[lininds>0]
            self.panline_matrix[lintemp-1, i]=1.0
            lintemp=lininds[lininds<0]
            self.panline_matrix[-lintemp-1, i]=-1.0
        self.panline_matrix=sps.csr_matrix(self.panline_matrix)
    def genaicm(self): #call FORTRAN backend to generate AIC matrix
        colmat=np.array([[0.0, 0.0, 0.0]]*self.npanels, dtype='double')
        for i in range(len(self.panels)):
            colmat[i, :]=self.panels[i].colpoint
        self.aicm3_line=toolkit.aicm_lines_gen(self.lines, colmat)
        self.gen_panline()
        self.aicm3=np.zeros((3, self.npanels, self.npanels))
        for i in range(3):
            self.aicm3[i, :, :]=self.aicm3_line[i, :, :]@self.panline_matrix
        self.aicm=np.zeros((self.npanels, self.npanels), dtype='double')
        for i in range(self.npanels):
            for j in range(self.npanels):
                self.aicm[i, j]=self.aicm3[:, i, j]@self.panels[i].nvector
    def genvbar(self, Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0): #generate freestream velocity vector. Angular velocities are raw (rad/s), not normalized by Uinf or dimensions
        for i in range(self.npanels):
            self.vbar[i, :]=Uinf*np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)], dtype='double')+np.cross(np.array([-p, q, -r], dtype='double'), self.panels[i].colpoint)
    def gennvv(self): #Generate vector containing normal velocities
        for i in range(self.npanels):
            self.nvv[i]=self.panels[i].nvector@self.vbar[i, :]
    def gen_selfinf(self):
        for i in range(self.npanels):
            self.delphi[i, :]+=toolkit.self_influence(self.lines, self.solution_line, self.panels[i].S, \
                self.panels[i].nvector, np.array([int(abs(l)) for l in self.panels[i].lines]), len(self.panels[i].wakelines)>0)
    def solve(self, damper=0.0): #generate Euler solution. Inverts AIC matrix with Tikhonov regularization if "damper" is set to non-zero value.
        #self.iaicm=toolkit.aicm_inversion(self.aicm, damper)
        if damper!=0.0:
            self.iaicm=slg.inv(self.aicm.T@self.aicm+damper*np.eye(self.npanels, self.npanels))@self.aicm.T
        else:
            self.iaicm=slg.inv(self.aicm)
        #self.solution=slg.solve(self.aicm, -self.nvv)
        self.solution=-self.iaicm@self.nvv
        self.solution_line=self.panline_matrix@self.solution
        for i in range(3):
            self.delphi[:, i]=self.aicm3[i, :, :]@self.solution #compute local velocity due to disturbance field
        self.gen_selfinf()
        self.solavailable=True
    def calcpress(self, Uinf=1.0):
        for i in range(self.npanels):
            self.Cps[i]=(Uinf**2-(self.delphi[i, :]+self.vbar[i, :])@(self.delphi[i, :]+self.vbar[i, :]))/Uinf**2
    def plotgeometry(self, xlim=[], ylim=[], zlim=[]):
        #plot geometry and local velocity vectors, either with or without wake panels
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        for i in range(self.nlines):
            ax.plot3D(self.lines[i, 0, :], self.lines[i, 1, :], self.lines[i, 2, :], 'gray')
        if self.solavailable:
            '''ax.quiver([p.colpoint[0] for p in self.panels], [p.colpoint[1] for p in self.panels], \
                [p.colpoint[2] for p in self.panels], [p.nvector[0]*0.005 for p in self.panels], \
                    [p.nvector[1]*0.005 for p in self.panels], \
                        [p.nvector[2]*0.005 for p in self.panels])'''
            ax.quiver([p.colpoint[0] for p in self.panels], [p.colpoint[1] for p in self.panels], \
                [p.colpoint[2] for p in self.panels], [self.vbar[i, 0]+self.delphi[i, 0] for i in range(self.npanels)], \
                    [self.vbar[i, 1]+self.delphi[i, 1] for i in range(self.npanels)], \
                        [self.vbar[i, 2]+self.delphi[i, 2] for i in range(self.npanels)])
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    def plotnormals(self, xlim=[], ylim=[], zlim=[], factor=1.0):
        #plot normal vectors of panels. Function essentially produced for geometry gen. debugging
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        for i in range(self.nlines):
            ax.plot3D(self.lines[i, 0, :], self.lines[i, 1, :], self.lines[i, 2, :], 'gray')
        ax.quiver([p.colpoint[0] for p in self.panels], [p.colpoint[1] for p in self.panels], \
            [p.colpoint[2] for p in self.panels], [p.nvector[0]*factor for p in self.panels], \
                [p.nvector[1]*factor for p in self.panels], \
                    [p.nvector[2]*factor for p in self.panels])
        '''ax.quiver([p.colpoint[0] for p in self.panels], [p.colpoint[1] for p in self.panels], \
            [p.colpoint[2] for p in self.panels], [self.vbar[i, 0]+self.delphi[i, 0] for i in range(self.npanels)], \
                [self.vbar[i, 1]+self.delphi[i, 1] for i in range(self.npanels)], \
                    [self.vbar[i, 2]+self.delphi[i, 2] for i in range(self.npanels)])'''
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
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
    def eulersolve(self, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, damper=0.0, Uinf=1.0):
        self.genvbar(Uinf, a=a, b=b, p=p, q=q, r=r)
        self.gennvv()
        self.genaicm()
        self.solve(damper=damper)
        self.calcpress(Uinf=Uinf)

'''sld=Solid(sldlist=[np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]).T])
#sld.plotgeometry()
sld.genaicm()
print(sld.aicm3[:, 0, 0])'''