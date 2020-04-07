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
import multiprocessing as mp
import os

import toolkit
from utils import *

class Panel: #Panel data type
    def __init__(self, lines):
        self.lines=lines
        self.wakelines=[]
        self.no_selfinf=[]
    def nocirc_enforce(self, linind):
        if linind!=-1:
            self.no_selfinf+=[linind]

def subprocess_genaicm(queue1, queue2): #generic function to unpack AICM calc order from multiprocessing 
    #queue and deliver it to calculating function in toolkit
    order=queue1.get()
    AICM=toolkit.aicm_lines_gen(order[2], order[3])
    queue2.put((order[0], order[1], AICM))

class WakeLine: #class encompassing info about a wake line. Made so as to compute deformation with wake rollup
    def __init__(self, ind): #recieves signed index as input
        self.ind=ind
        self.v=np.zeros(3)
        self.nabut=0
        self.updateme=True
        self.inverted=False
    def addvel(self, v): #add velocity to be averaged
        self.v+=v
        self.nabut+=1

class WakePanel:
    def __init__(self, wl, wr, center): #add a wake panel based on adjacent lines
        self.wl=wl
        self.wr=wr
        self.center=center
        self.v=np.zeros(3)

class Solid:
    def __init__(self, sldlist=[], wraparounds=[]): #solid data type
        self.panels=[]
        self.lines=[]
        self.addto=[]
        self.wakestrips=[]
        self.wakeline_inds=set([])
        self.wakelines=set([])
        self.problematic=[]
        self.npanels=0
        self.nlines=0
        self.nwake=0 #number of wake panels
        self.solavailable=False
        if len(wraparounds)==0:
            wraparounds=[[]]*len(sldlist)
        sldcnt=0
        for sld in sldlist: #define first patches (without their interconnections)
            wraps=wraparounds[sldcnt]
            self.addpatch(sld, wraps=wraps)
            sldcnt+=1
    def addpatch(self, sld, wraps=[], prevlines={}, invlats=[], tolerance=0.00005): #add panel based on point grid (list of lists)
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
                        horzlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i][j+1])).T, tolerance=tolerance)]
                else:
                    horzlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i][j+1])).T, tolerance=tolerance)]
                if j==0 and 'right' in prevlines:
                    if prevlines['right'][i]!=-2:
                        vertlines[-1]+=[prevlines['right'][i]]
                    else:
                        vertlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i+1][j])).T, tolerance=tolerance)]
                else:
                    vertlines[-1]+=[self.addline(np.vstack((sld[i][j], sld[i+1][j])).T, tolerance=tolerance)]
        if 0 in wraps: #if wrapped in x direction (inner list layer), repeat the first lateral vertical lines
            for i in range(len(sld)-1):
                vertlines[i]+=[vertlines[i][0]]
        else:
            if 'left' in prevlines:
                for i in range(len(sld)-1):
                    if prevlines['left'][i]!=-2:
                        vertlines[i]+=[prevlines['left'][i]]
                    else:
                        vertlines[i]+=[self.addline(np.vstack((sld[i][-1], sld[i+1][-1])).T, tolerance=tolerance)]
            else:
                for i in range(len(sld)-1):
                    vertlines[i]+=[self.addline(np.vstack((sld[i][-1], sld[i+1][-1])).T, tolerance=tolerance)]
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
                        horzlines[-1]+=[self.addline(np.vstack((sld[-1][i], sld[-1][i+1])).T, tolerance=tolerance)]
            else:
                for i in range(len(sld[0])-1):
                    horzlines[-1]+=[self.addline(np.vstack((sld[-1][i], sld[-1][i+1])).T, tolerance=tolerance)]
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
                paninds[-1]+=[self.addpanel([horzlines[i][j], vertlines[i][j+1], horzlines[i+1][j], vertlines[i][j]], invs=l)]
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
                if (n==3 or n==2):
                    llist+=[-1-l]
                else:
                    llist+=[l+1]
                if n in invs:
                    llist[-1]=-llist[-1]
            n+=1
        self.panels+=[Panel(llist)]
        return self.npanels-1
    def addwakepanel(self, refup, refdown, indup=0, indown=2, offsetleft=np.array([1000.0, 0.0, 0.0]), offsetright=np.array([]), \
        tolerance=0.00005, disc=10, strategy=lambda x: ((np.exp(x)-1.0)/(exp(1)-1.0))**2, prevleft=[], leftinvert=False, prevright=[], rightinvert=False):
        #add line segments correspondent to panel wake, and remove "inup-th"/"indown-th" line segment in the panel's list
        if len(offsetright)==0:
            offsetright=offsetleft
        self.nwake+=1
        wakedisc=strategy(np.linspace(0.0, 1.0, disc+1))
        if indup<len(self.panels[refup].lines):
            te_linind=self.panels[refup].lines.pop(indup)
        else:
            te_linind=self.panels[refup].wakelines.pop(indup-len(self.panels[refup].lines))
        self.panels[refup].TE_line=te_linind
        coords=self.line_getcoords(te_linind)
        self.wakestrips+=[[]]
        indsup=[]
        leftlines=[]
        rightlines=[]
        if len(prevright)==0:
            for i in range(disc):
                p0=coords[:, 0]+offsetright*wakedisc[i]
                p1=coords[:, 0]+offsetright*wakedisc[i+1]
                indsup+=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)+1]
                rightlines+=[WakeLine(indsup[-1])]
        else:
            if rightinvert:
                rightlines=[prevright[i] for i in range(len(prevright)-1, -1, -1)]
                indsup+=[-prevright[i].ind for i in range(len(prevright)-1, -1, -1)]
            else:
                indsup+=[prevright[i].ind for i in range(len(prevright))]
                rightlines=prevright
        if len(prevleft)==0:
            for i in range(disc, 0, -1):
                p0=coords[:, 1]+offsetleft*wakedisc[i]
                p1=coords[:, 1]+offsetleft*wakedisc[i-1]
                indsup+=[self.addline(np.vstack((p0, p1)).T, tolerance=tolerance)+1]
                leftlines+=[WakeLine(indsup[-1])]
        else:
            if leftinvert:
                leftlines=[prevleft[i] for i in range(len(prevleft)-1, -1, -1)]
                indsup+=[-prevleft[i].ind for i in range(len(prevleft)-1, -1, -1)]
            else:
                indsup+=[prevleft[i].ind for i in range(len(prevleft))]
                leftlines=prevleft
        self.panels[refup].wakelines+=indsup
        for l in leftlines+rightlines:
            self.wakelines.add(l)
        for lind in indsup:
            self.wakeline_inds.add(abs(lind)-1)
        if refdown!=-1:
            te_linind=self.panels[refdown].lines.pop(indown)
            self.panels[refdown].wakelines+=[-l for l in indsup]
            self.panels[refdown].TE_line=te_linind
        for i in range(disc):
            center=(self.line_midpoint(rightlines[i].ind)+self.line_midpoint(leftlines[disc-i-1].ind))/2
            self.wakestrips[-1]+=[WakePanel(leftlines[disc-i-1], rightlines[i], center)]
        return leftlines, rightlines #return left, right wake lines
    def genwakepanels(self, wakecombs=[], wakeinds=[], offset=1000.0, a=0.0, b=0.0, disc=10, strategy=lambda x: ((np.exp(x)-1.0)/(exp(1)-1.0))**2, \
        prevleft=[], prevright=[], leftinvert=False, rightinvert=False):
        #generate wake panels based in list of TE panel combinations
        #wakecombs: list of lists, first element in sublist is upper surface panel
        #wakeinds: list of lists indicating corresponding wake vortex line segment index to apply kutta condition to
        if len(wakeinds)==0:
            wakeinds=[[0, 2]]
        trimlist(len(wakecombs), wakeinds)
        self.addto+=[[wakecombs[0][0], wakecombs[0][1]]]
        left, prevleft=self.addwakepanel(wakecombs[0][0], wakecombs[0][1], indup=wakeinds[0][0], indown=wakeinds[0][0], \
            offsetleft=np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)])*offset, disc=disc, strategy=strategy, prevleft=prevleft, leftinvert=leftinvert)
        for i in range(1, len(wakecombs)-1):
            self.addto+=[[wakecombs[i][0], wakecombs[i][1]]]
            _, prevleft=self.addwakepanel(wakecombs[i][0], wakecombs[i][1], indup=wakeinds[i][0], indown=wakeinds[i][0], \
                offsetleft=np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)])*offset, disc=disc, strategy=strategy, \
                    prevleft=prevleft, leftinvert=True)
        self.addto+=[[wakecombs[-1][0], wakecombs[-1][1]]]
        _, prevleft=self.addwakepanel(wakecombs[-1][0], wakecombs[-1][1], indup=wakeinds[-1][0], indown=wakeinds[-1][0], \
            offsetleft=np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)])*offset, disc=disc, strategy=strategy, \
                prevleft=prevleft, leftinvert=True, prevright=prevright, rightinvert=rightinvert)
        return left, prevleft
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
        coords=coords.T
        if np.size(coords, 1)==3:
            temp=coords
            coords=np.zeros((3, 4))
            coords[:, 0:3]=temp
            coords[:, -1]=temp[:, -1]
        return coords
    def panel_calcSn(self, p): #update panel area and normal vector
        l=self.line_getcoords(p.lines[0])
        p.nvector=np.cross(l[:, 1]-p.colpoint, l[:, 0]-p.colpoint)
        p.S=lg.norm(p.nvector)
        p.nvector/=p.S
        p.S/=2
        for i in range(1, len(p.lines)):
            l=self.line_getcoords(p.lines[i])
            p.S+=lg.norm(np.cross(l[:, 1]-p.colpoint, l[:, 0]-p.colpoint))/2
    def line_getcoords(self, ind): #return line coordinates in order presented in panel, based on panel.lines element
        if ind>0:
            return self.lines[ind-1, :, :]
        else:
            return np.flip(self.lines[-1-ind, :, :], axis=1)
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
                    if (u@np.cross(self.panels[i].nvector, p)*self.panels[i].lines[lind]>0.0):
                        self.panels[i].lines[lind]*=-1
                        print('WARNING: line '+str(abs(self.panels[i].lines[lind])-1)+' in panel '+\
                            str(i)+' had to be inverted. Please check integrity of patchcompose() functions')
                        self.problematic+=[abs(l)-1 for l in self.panels[i].lines]
    def iscontiguous(self, patchinds=[], tolerance=5e-5): #check if all panels in solid are closed quadrilaterals
        if len(patchinds)==0:
            patchinds=[list(range(self.npanels))]
        for patchlist in patchinds:
            for i in patchlist:
                u=self.line_getcoords(self.panels[i].lines[0])
                for lind in range(1, len(self.panels[i].lines)):
                    v=u
                    u=self.line_getcoords(self.panels[i].lines[lind])
                    if np.amax(np.abs(v[:, 1]-u[:, 0]))>tolerance:
                        print('WARNING: '+str(i)+' panel is not contiguous')
                        self.problematic+=[abs(l)-1 for l in self.panels[i].lines]
                        break
    def end_preprocess(self, paninds=[], tolerance=5e-5, initialize=True, check=True):
        #calculate panel areas before they are altered by wake generation. Must be run before it, and terrible
        #consequences may arise if done otherwise
        u=np.array([0.0, 0.0, 0.0])
        v=np.array([0.0, 0.0, 0.0])
        #analyse whether provided panel index list corresponds to a certain set or, as default, indicates all panels in self
        allpans=(len(paninds)==0) or (len(paninds)==self.npanels)
        if allpans:
            panlist=self.panels
        else:
            panlist=[self.panels[i] for i in paninds]
        n=0
        for p in panlist:
            '''coords=self.panel_getcoords(p)
            u=coords[:, 1]-coords[:, 0]
            v=coords[:, 2]-coords[:, 1]'''
            u=self.line_getvec(p.lines[0])
            v=self.line_getvec(p.lines[1])
            p.nvector=np.cross(v, u)
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
                p.S=0.0
            if p.S<=tolerance**2:
                print('WARNING: panel '+str(n)+' displayed incoherent area '+str(p.S)+' with '+str(len(p.lines))+' lines')
            n+=1
        #initialize result vectors
        if initialize:
            self.delphi=np.array([[0.0, 0.0, 0.0]]*len(self.panels), dtype='double')
            self.vbar=np.array([[0.0, 0.0, 0.0]]*len(self.panels), dtype='double')
            self.nvv=np.array([0.0]*len(self.panels), dtype='double')
            self.solution=np.array([0.0]*len(self.panels), dtype='double')
            self.solavailable=False
            self.Cps=np.array([0.0]*len(self.panels), dtype='double')
            self.Cfs=np.zeros(len(self.panels), dtype='double')
            self.forces=[]
            self.moments=[]
        #adjust lines in case any is set inconsistently with respect to anti-clockwise convention in panel
        if check:
            self.lineadjust()
            #self.iscontiguous(tolerance=tolerance)
    def gen_panline(self): #generate panel-line correspondence matrix
        self.panline_matrix=np.zeros((self.nlines, self.npanels), dtype='double')
        for i in range(self.npanels):
            lininds=np.array(self.panels[i].lines+self.panels[i].wakelines)
            lintemp=lininds[lininds>0]
            self.panline_matrix[lintemp-1, i]=1.0
            lintemp=lininds[lininds<0]
            self.panline_matrix[-lintemp-1, i]=-1.0
        self.panline_matrix=sps.csr_matrix(self.panline_matrix)
    def addorder(self, queue, colmat, ind1, ind2): #add order to queue for subprocess to compute AICM from rows ind1 to ind2
        queue.put((ind1, ind2, self.lines[ind1:ind2, :, :], colmat))
    def genaicm(self): #call FORTRAN backend to generate AIC matrix
        colmat=np.zeros((self.npanels, 3))
        nvectmat=np.zeros((self.npanels, 3))
        for i in range(len(self.panels)):
            colmat[i, :]=self.panels[i].colpoint
            nvectmat[i, :]=self.panels[i].nvector
        ncpus=mp.cpu_count()
        calclims=np.linspace(0, self.nlines, ncpus+1)
        calclims=[int(c) for c in calclims]
        calclims=[[calclims[i], calclims[i+1]] for i in range(ncpus)]
        queue1=mp.Queue()
        queue2=mp.Queue()
        for i in range(len(calclims)):
            self.addorder(queue1, colmat, ind1=calclims[i][0], ind2=calclims[i][1])
        processes=[]
        for i in range(ncpus):
            processes+=[mp.Process(target=subprocess_genaicm, args=(queue1, queue2))]
        for p in processes:
            #p.daemon=True
            p.start()
        ordresults=[]
        for i in range(ncpus):
            ordresults+=[queue2.get()]
        for p in processes:
            p.join()
        self.aicm3_line=np.zeros((3, self.npanels, self.nlines))
        for ordresult in ordresults:
            self.aicm3_line[:, :, ordresult[0]:ordresult[1]]=ordresult[2]
        self.gen_panline()
        self.aicm3=np.zeros((3, self.npanels, self.npanels))
        for i in range(3):
            self.aicm3[i, :, :]=self.aicm3_line[i, :, :]@self.panline_matrix
        self.aicm=toolkit.aicm_norm_conv(self.aicm3, nvectmat)
    def gen_farfield(self, Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0): #generate generic local freestream velocity dependant on parameters
        newvec=np.zeros((self.npanels, 3))
        for i in range(self.npanels):
            newvec[i, :]=Uinf*np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)], dtype='double')+np.cross(np.array([p, -q, r], dtype='double'), self.panels[i].colpoint)
        return newvec
    def gen_wake_farfield(self, Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0): #generate generic local freestream at wake panels
        for strip in self.wakestrips:
            for pan in strip:
                pan.v=Uinf*np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)], dtype='double')+np.cross(np.array([p, -q, r], dtype='double'), pan.center)
    def gen_farfield_derivative(self, Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, par='a'): #generate generic local freestream velocity dependant on parameters
        newvec=np.zeros((self.npanels, 3))
        if par=='a':
            for i in range(self.npanels):
                newvec[i, :]=Uinf*np.array([-sin(a)*cos(b), sin(a)*sin(b), cos(a)], dtype='double')
        elif par=='b':
            for i in range(self.npanels):
                newvec[i, :]=Uinf*np.array([-cos(a)*sin(b), -cos(a)*cos(b), 0.0], dtype='double')
        elif par=='p':
            for i in range(self.npanels):
                newvec[i, :]=np.cross(np.array([1.0, 0.0, 0.0], dtype='double'), self.panels[i].colpoint)
        elif par=='q':
            for i in range(self.npanels):
                newvec[i, :]=np.cross(np.array([0.0, -1.0, 0.0], dtype='double'), self.panels[i].colpoint)
        elif par=='r':
            for i in range(self.npanels):
                newvec[i, :]=np.cross(np.array([0.0, 0.0, 1.0], dtype='double'), self.panels[i].colpoint)
        elif par=='Uinf':
            for i in range(self.npanels):
                newvec[i, :]=np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)], dtype='double')
        return newvec
    def genvbar(self, Uinf, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0): #generate freestream velocity vector. Angular velocities are raw (rad/s), not normalized by Uinf or dimensions
        self.vbar=self.gen_farfield(Uinf, a=a, b=b, p=p, q=q, r=r)
    def gennvv(self, beta=1.0): #Generate vector containing normal velocities
        for i in range(self.npanels):
            self.nvv[i]=self.panels[i].nvector@self.vbar[i, :]
    def gen_selfinf_mat(self): #generate self-influence matrix according to Srivastava's equations, with each line's vorticity as an input
        #so as to ease stability derivative calculation
        self.selfinf_mat_x=np.zeros((self.npanels, self.nlines), order='F')
        self.selfinf_mat_y=np.zeros((self.npanels, self.nlines), order='F')
        self.selfinf_mat_z=np.zeros((self.npanels, self.nlines), order='F')
        for i in range(self.npanels):
            linelist=[]
            nocirc_linelist=[]
            for l in self.panels[i].lines:
                if not int(abs(l))-1 in self.panels[i].no_selfinf:
                    linelist+=[int(abs(l))-1]
                else:
                    nocirc_linelist+=[int(abs(l))-1]
            for l in linelist:
                vdv=np.cross(self.lines[l, :, 1]-self.lines[l, :, 0], self.panels[i].nvector)/(4*self.panels[i].S)
                self.selfinf_mat_x[i, l]=vdv[0]
                self.selfinf_mat_y[i, l]=vdv[1]
                self.selfinf_mat_z[i, l]=vdv[2]
            for l in nocirc_linelist:
                vdv=np.cross(self.lines[l, :, 1]-self.lines[l, :, 0], self.panels[i].nvector)/(4*self.panels[i].S)
                vdv/=2
                self.selfinf_mat_x[i, l]=vdv[0]
                self.selfinf_mat_y[i, l]=vdv[1]
                self.selfinf_mat_z[i, l]=vdv[2]
        self.selfinf_mat_x=sps.csr_matrix(self.selfinf_mat_x)
        self.selfinf_mat_y=sps.csr_matrix(self.selfinf_mat_y)
        self.selfinf_mat_z=sps.csr_matrix(self.selfinf_mat_z)
    def gen_selfinf(self): #generate self-influence velocity according to Srivastava's equations
        self.delphi[:, 0]+=self.selfinf_mat_x@self.solution_line
        self.delphi[:, 1]+=self.selfinf_mat_y@self.solution_line
        self.delphi[:, 2]+=self.selfinf_mat_z@self.solution_line
    def solve(self, damper=0.0, target=np.array([]), wakeiter=0, Uinf=1.0, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, tolerance=1e-5, echo=True):
        #generate Euler solution. Inverts AIC matrix with Tikhonov regularization if "damper" is set to non-zero value.
        #self.iaicm=toolkit.aicm_inversion(self.aicm, damper)
        if len(target)==0:
            target=np.zeros(self.npanels) # target: nvv aimed at, for stability derivative computation and/or viscous-inviscid coupling

        if damper!=0.0:
            self.iaicm=slg.inv(self.aicm.T@self.aicm+damper*np.eye(self.npanels, self.npanels))@self.aicm.T
        else:
            self.iaicm=slg.inv(self.aicm)
        tango=target-self.nvv

        self.solution=self.iaicm@tango
        self.solution_line=self.panline_matrix@self.solution
        self.solavailable=True
        if wakeiter!=0:
            wakelininds=np.array(list(self.wakeline_inds))
            colmat=np.zeros((self.npanels, 3))
            nvectmat=np.zeros((self.npanels, 3))
            for i in range(len(self.panels)):
                colmat[i, :]=self.panels[i].colpoint
                nvectmat[i, :]=self.panels[i].nvector
            if echo:
                print('%18s | %5s %5s %5s' % ('Iteration', 'vcalc', 'roll', 'solve'))
        for i in range(wakeiter): #iterate for wake generation
            t0=tm.time()
            self.gen_wake_farfield(Uinf=Uinf, a=a, b=b, p=p, q=q, r=r)
            self.add_wakevels(tolerance=tolerance)
            tvel=tm.time()
            self.wake_rollup()
            trol=tm.time()
            self.aicm3_line[:, :, wakelininds]=toolkit.aicm_lines_recalc(wakelininds+1, self.lines, colmat)
            for j in range(3):
                self.aicm3[j, :, :]=self.aicm3_line[j, :, :]@self.panline_matrix
            self.aicm=toolkit.aicm_norm_conv(self.aicm3, nvectmat)
            if damper!=0.0:
                self.iaicm=slg.inv(self.aicm.T@self.aicm+damper*np.eye(self.npanels, self.npanels))@self.aicm.T
            else:
                self.iaicm=slg.inv(self.aicm)
            self.solution=self.iaicm@tango
            self.solution_line=self.panline_matrix@self.solution
            tsol=tm.time()
            print('%14s %3d | %5f %5f %5f' % ('Wake teration', i, tvel-t0, trol-tvel, tsol-trol))
            
        for i in range(3):
            self.delphi[:, i]=self.aicm3[i, :, :]@self.solution #compute local velocity due to disturbance field
        self.gen_selfinf()
    def calcpress(self, Uinf=1.0, gamma=1.4, M=0.0):
        if M==0.0:
            for i in range(self.npanels):
                self.Cps[i]=(Uinf**2-(self.delphi[i, :]+self.vbar[i, :])@(self.delphi[i, :]+self.vbar[i, :]))/Uinf**2
        else:
            for i in range(self.npanels):
                incomp=(Uinf**2-(self.delphi[i, :]+self.vbar[i, :])@(self.delphi[i, :]+self.vbar[i, :]))/Uinf**2 #incompressible value
                self.Cps[i]=2*((1.0+(gamma-1)*M**2*incomp/2)**(gamma/(gamma-1))-1)/(M**2*gamma)
    def calcforces(self): #compute force correspondent to unitary dynamic pressure on each panel
        self.forces=[-self.panels[i].S*self.panels[i].nvector*self.Cps[i]+self.panels[i].S*self.Cfs[i]*\
            (self.vbar[i, :]+self.delphi[i, :])/lg.norm(self.vbar[i, :]+self.delphi[i, :]) for i in range(self.npanels)]
        self.moments=[np.cross(self.panels[i].colpoint, self.forces[i]) for i in range(self.npanels)]
    def calc_derivative_dv(self, Uinf, dvdksi, M=0.0, gamma=1.4): #calculate local Cp derivative by freestream derivative
        dndksi=np.array([self.panels[i].nvector@dvdksi[i, :] for i in range(self.npanels)])
        dGammadksi=-self.iaicm@dndksi
        dGamma_linedksi=self.panline_matrix@dGammadksi
        dvdksi[:, 0]+=self.aicm3[0, :, :]@dGammadksi+self.selfinf_mat_x@dGamma_linedksi
        dvdksi[:, 1]+=self.aicm3[1, :, :]@dGammadksi+self.selfinf_mat_y@dGamma_linedksi
        dvdksi[:, 2]+=self.aicm3[2, :, :]@dGammadksi+self.selfinf_mat_z@dGamma_linedksi
        dCps=np.array([-(2*(self.vbar[i, :]+self.delphi[i, :])@dvdksi[i, :])/Uinf**2 for i in range(self.npanels)])
        if M!=0.0:
            for i in range(self.npanels):
                incomp=(Uinf**2-(self.delphi[i, :]+self.vbar[i, :])@(self.delphi[i, :]+self.vbar[i, :]))/Uinf**2 #incompressible value
                dCps[i]=2*((1+(gamma-1)*M**2*incomp/2)**(1.0/(gamma-1.0)))*dCps[i]/(M**2*(gamma-1.0))
        return dCps
    def calc_derivative_dn(self, Uinf, dndksi, M=0.0, gamma=1.4): #calculate local Cp derivative by normal velocity derivative
        dGammadksi=-self.iaicm@dndksi
        dGamma_linedksi=self.panline_matrix@dGammadksi
        dvdksi=np.zeros((len(dndksi), 3))
        dvdksi[:, 0]=self.aicm3[0, :, :]@dGammadksi+self.selfinf_mat_x@dGamma_linedksi
        dvdksi[:, 1]=self.aicm3[1, :, :]@dGammadksi+self.selfinf_mat_y@dGamma_linedksi
        dvdksi[:, 2]=self.aicm3[2, :, :]@dGammadksi+self.selfinf_mat_z@dGamma_linedksi
        dCps=np.array([-(2*(self.vbar[i, :]+self.delphi[i, :])@dvdksi[i, :])/Uinf**2 for i in range(self.npanels)])
        if M!=0.0:
            for i in range(self.npanels):
                incomp=(Uinf**2-(self.delphi[i, :]+self.vbar[i, :])@(self.delphi[i, :]+self.vbar[i, :]))/Uinf**2 #incompressible value
                dCps[i]=2*((1+(gamma-1)*M**2*incomp/2)**(1.0/(gamma-1.0)))*dCps[i]/(M**2*(gamma-1.0))
        return dCps
    def add_wakevels(self, tolerance=1e-5): #calculate velocities at wake panel control points and add to local velocity variables
        for strip in self.wakestrips:
            for p in strip:
                p.v+=toolkit.get_field_influence(self.lines, self.solution_line, p.center, tolerance=tolerance)
    def wake_rollup(self): #calculate velocity at nodes located at wake line midpoints and update the lines
        for strip in self.wakestrips:
            for p in strip:
                p.wl.addvel(p.v)
                p.wr.addvel(p.v)
        self.wakelines_deform()
    def wakelines_deform(self): #deform wake lines according to computed velocity at panels between them
        for l in self.wakelines:
            l.v/=l.nabut
            l.inverted=self.lines[abs(l.ind)-1, 0, 0]>self.lines[abs(l.ind)-1, 0, 1]
            l.updateme=True
        deltaxl=0.0
        deltaxr=0.0
        dvl=np.zeros(3)
        dvr=np.zeros(3)
        origin_left=np.zeros(3)
        origin_right=np.zeros(3)
        newleft=np.zeros(3)
        newright=np.zeros(3)
        for strip in self.wakestrips:
            origin_left=self.lines[abs(strip[0].wl.ind)-1, :, 1 if strip[0].wl.inverted else 0]
            origin_right=self.lines[abs(strip[0].wr.ind)-1, :, 1 if strip[0].wr.inverted else 0]
            for i in range(len(strip)):
                if strip[i].wl.updateme:
                    deltaxl=abs(self.lines[abs(strip[i].wl.ind)-1, 0, 1]-self.lines[abs(strip[i].wl.ind)-1, 0, 0])
                    dvl=strip[i].wl.v*deltaxl/(strip[i].wl.v[0])
                    newleft=origin_left+dvl
                    if strip[i].wl.inverted:
                        self.lines[abs(strip[i].wl.ind)-1, :, 1]=origin_left
                        self.lines[abs(strip[i].wl.ind)-1, :, 0]=newleft
                    else:
                        self.lines[abs(strip[i].wl.ind)-1, :, 0]=origin_left
                        self.lines[abs(strip[i].wl.ind)-1, :, 1]=newleft
                    origin_left=newleft
                    strip[i].wl.updateme=False
                if strip[i].wr.updateme:
                    deltaxr=abs(self.lines[abs(strip[i].wr.ind)-1, 0, 1]-self.lines[abs(strip[i].wr.ind)-1, 0, 0])
                    dvr=strip[i].wr.v*deltaxr/(strip[i].wr.v[0])
                    newright=origin_right+dvr
                    if strip[i].wr.inverted:
                        self.lines[abs(strip[i].wr.ind)-1, :, 1]=origin_right
                        self.lines[abs(strip[i].wr.ind)-1, :, 0]=newright
                    else:
                        self.lines[abs(strip[i].wr.ind)-1, :, 0]=origin_right
                        self.lines[abs(strip[i].wr.ind)-1, :, 1]=newright
                    origin_right=newright
                    strip[i].wr.updateme=False
                strip[i].center=(self.line_midpoint(strip[i].wr.ind)+self.line_midpoint(strip[i].wl.ind))/2
        for l in self.wakelines:
            l.v=0.0
            l.nabut=0
    def PG_apply(self, beta, a, b): #beta indicates PG correction factor for compressibility
        PG_apmat=PG_xmult(beta, a, b)
        #convert lines
        self.lines[:, 0, :]=PG_apmat[0, 0]*self.lines[:, 0, :]+PG_apmat[0, 1]*self.lines[:, 1, :]+PG_apmat[0, 2]*self.lines[:, 2, :]
        self.lines[:, 1, :]=PG_apmat[1, 0]*self.lines[:, 0, :]+PG_apmat[1, 1]*self.lines[:, 1, :]+PG_apmat[1, 2]*self.lines[:, 2, :]
        self.lines[:, 2, :]=PG_apmat[2, 0]*self.lines[:, 0, :]+PG_apmat[2, 1]*self.lines[:, 1, :]+PG_apmat[2, 2]*self.lines[:, 2, :]
        #convert panels
        for p in self.panels:
            p.colpoint=PG_apmat@p.colpoint
            self.panel_calcSn(p)
        #convert wake panels
        for strip in self.wakestrips:
            for p in strip:
                p.center=PG_apmat@p.center
    def PG_remove(self, beta, a, b): #remove effect of previously imposed PG correction
        PG_apmat=PG_inv_xmult(beta, a, b)
        PG_apv=PG_vtouni(beta, a, b)
        #convert lines
        self.lines[:, 0, :]=PG_apmat[0, 0]*self.lines[:, 0, :]+PG_apmat[0, 1]*self.lines[:, 1, :]+PG_apmat[0, 2]*self.lines[:, 2, :]
        self.lines[:, 1, :]=PG_apmat[1, 0]*self.lines[:, 0, :]+PG_apmat[1, 1]*self.lines[:, 1, :]+PG_apmat[1, 2]*self.lines[:, 2, :]
        self.lines[:, 2, :]=PG_apmat[2, 0]*self.lines[:, 0, :]+PG_apmat[2, 1]*self.lines[:, 1, :]+PG_apmat[2, 2]*self.lines[:, 2, :]
        #convert panels
        for p in self.panels:
            p.colpoint=PG_apmat@p.colpoint
            self.panel_calcSn(p)
        #convert wake panels
        for strip in self.wakestrips:
            for p in strip:
                p.center=PG_apmat@p.center
        #convert self influence and aicm3_line matrixes
        self.selfinf_mat_x=PG_apv[0, 0]*self.selfinf_mat_x+PG_apv[0, 1]*self.selfinf_mat_y+PG_apv[0, 2]*self.selfinf_mat_z
        self.selfinf_mat_y=PG_apv[1, 0]*self.selfinf_mat_x+PG_apv[1, 1]*self.selfinf_mat_y+PG_apv[1, 2]*self.selfinf_mat_z
        self.selfinf_mat_z=PG_apv[2, 0]*self.selfinf_mat_x+PG_apv[2, 1]*self.selfinf_mat_y+PG_apv[2, 2]*self.selfinf_mat_z
        #self.aicm remains the same
        self.aicm3_line[0, :, :]=PG_apv[0, 0]*self.aicm3_line[0, :, :]+PG_apv[0, 1]*self.aicm3_line[1, :, :]+PG_apv[0, 2]*self.aicm3_line[2, :, :]
        self.aicm3_line[1, :, :]=PG_apv[1, 0]*self.aicm3_line[0, :, :]+PG_apv[1, 1]*self.aicm3_line[1, :, :]+PG_apv[1, 2]*self.aicm3_line[2, :, :]
        self.aicm3_line[2, :, :]=PG_apv[2, 0]*self.aicm3_line[0, :, :]+PG_apv[2, 1]*self.aicm3_line[1, :, :]+PG_apv[2, 2]*self.aicm3_line[2, :, :]
        self.aicm3[0, :, :]=PG_apv[0, 0]*self.aicm3[0, :, :]+PG_apv[0, 1]*self.aicm3[1, :, :]+PG_apv[0, 2]*self.aicm3[2, :, :]
        self.aicm3[1, :, :]=PG_apv[1, 0]*self.aicm3[0, :, :]+PG_apv[1, 1]*self.aicm3[1, :, :]+PG_apv[1, 2]*self.aicm3[2, :, :]
        self.aicm3[2, :, :]=PG_apv[2, 0]*self.aicm3[0, :, :]+PG_apv[2, 1]*self.aicm3[1, :, :]+PG_apv[2, 2]*self.aicm3[2, :, :]
    def plotgeometry(self, xlim=[], ylim=[], zlim=[], velfield=True, factor=1.0):
        #plot geometry and local velocity vectors, either with or without wake panels
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        for i in range(self.nlines):
            ax.plot3D(self.lines[i, 0, :], self.lines[i, 1, :], self.lines[i, 2, :], 'gray')
        for i in self.problematic:
            ax.plot3D(self.lines[i, 0, :], self.lines[i, 1, :], self.lines[i, 2, :], 'red')
        if self.solavailable and velfield:
            '''ax.quiver([p.colpoint[0] for p in self.panels], [p.colpoint[1] for p in self.panels], \
                [p.colpoint[2] for p in self.panels], [p.nvector[0]*0.005 for p in self.panels], \
                    [p.nvector[1]*0.005 for p in self.panels], \
                        [p.nvector[2]*0.005 for p in self.panels])'''
            ax.quiver([p.colpoint[0] for p in self.panels], [p.colpoint[1] for p in self.panels], \
                [p.colpoint[2] for p in self.panels], [(self.vbar[i, 0]+self.delphi[i, 0])*factor for i in range(self.npanels)], \
                    [(self.vbar[i, 1]+self.delphi[i, 1])*factor for i in range(self.npanels)], \
                        [(self.vbar[i, 2]+self.delphi[i, 2])*factor for i in range(self.npanels)])
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        ax.view_init(azim=-135, elev=30)
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
        ax.view_init(azim=-135, elev=30)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    def eulersolve(self, target=np.array([]), Uinf=1.0, M=0.0, gamma=1.4, beta=1.0, a=0.0, b=0.0, p=0.0, q=0.0, r=0.0, damper=0.0, echo=True, \
        wakeiter=0, tolerance=1e-5):
        if echo:
            print('========Euler solution=======')
            print(self.npanels, ' panels, ', self.nlines, ' lines')
        t=tm.time()
        self.genvbar(Uinf, a=a, b=b, p=p, q=q, r=r)
        self.gennvv()
        self.PG_apply(beta, a, b)
        if echo:
            print('Pre-processing: '+str(tm.time()-t)+' s')
        t=tm.time()
        self.genaicm()
        self.gen_selfinf_mat()
        if echo:
            print('AIC matrix generation: '+str(tm.time()-t)+' s')
        t=tm.time()
        self.PG_remove(beta, a, b)
        self.solve(damper=damper, wakeiter=wakeiter, Uinf=Uinf, a=a, b=b, p=p, q=q, r=r, tolerance=tolerance, echo=echo)
        self.calcpress(Uinf=Uinf, M=M, gamma=gamma)
        self.calcforces()
        if echo:
            print('Solution and post-processing: '+str(tm.time()-t)+' s')
            print('=============================')

'''sld=Solid(sldlist=[np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]).T])
#sld.plotgeometry()
sld.genaicm()
print(sld.aicm3[:, 0, 0])'''