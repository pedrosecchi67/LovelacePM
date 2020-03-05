import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
import scipy.sparse.linalg as splg
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import scipy.sparse as sps
import scipy.interpolate as sinterp
import time as tm

import toolkit
from utils import *

from paneller import *

#class to define functions and variables needed to represent body sections
class body_section:
    def __init__(self, center=np.array([0.0, 0.0, 0.0]), coords=np.vstack((np.sin(np.linspace(0.0, 2*pi, 360)), \
        np.cos(np.linspace(0.0, 2*pi, 360)))).T, R=1.0, cubic=True):
        self.R=R
        self.coords=np.vstack((np.zeros(np.size(coords, 0)), coords[:, 0], coords[:, 1])).T*R #redefine coordinate system
        self.center=center
        for i in range(np.size(coords, 0)):
            self.coords[i, :]+=center
        thetas=np.arctan2(coords[:, 0], coords[:, 1])
        Rs=lg.norm(coords, axis=1)
        order=np.argsort(thetas)
        self.thetas=thetas[order]
        self.Rs=Rs[order]
        if cubic:
            self.polar_rule=sinterp.CubicSpline(self.thetas, self.Rs, extrapolate=True)
        else:
            self.polar_rule=sinterp.interp1d(self.thetas, self.Rs)
        #define function for interpolation in polar coordinates between theta in respect to z axis as input
    def __call__(self, th):
        R=self.polar_rule(th)
        return np.array([0.0, sin(th), cos(th)])*R*self.R+self.center

class body_panel: #exclusively for intersection finding purposes
    def __init__(self, p1, p2, p3, p4, tolerance=0.00005):
        points=np.vstack((p1, p2, p3, p4)).T
        self.points=points
        self.center, self.Mtosys, self.Mtouni, self.locpoints=toolkit.body_panel_process(self.points, tolerance)
        #FORTRAN was used here for speed up purposes

def prevline_organize(queue, nlines, prevlateral=[], intra=False, right=False):
    #organize prevline dictionary element for body's patchcompose
    #prevlateral is a list including previously generated, adjacent patch edge lines
    #if no prevlateral is created, new patch edge lines will be created (-2 input argument for sld.addline)
    linelist=[]
    for i in range(len(queue)):
        if len(prevlateral)==0:
            linelist+=[-2]*nlines[i]
        else:
            linelist+=[prevlateral[j] for j in range(len(linelist), nlines[i]+len(linelist))]
        if intra:
            if right:
                linelist+=queue[i].intraright_lines
            else:
                linelist+=queue[i].intraleft_lines
        else:
            if right:
                linelist+=queue[i].extraright_lines
            else:
                linelist+=queue[i].extraleft_lines
    if len(prevlateral)==0:
        linelist+=nlines[-1]*[-2]
    else:
        linelist+=[prevlateral[i] for i in range(len(linelist), nlines[-1]+len(linelist))]
    return linelist

class body: #body definition class for center definition to obtain polar cooridinates and creation of input panels 
    #for abutments.
    #to be abutted stern to bow, clockwise when seen from the aircraft's bow
    def __init__(self, sld, sections=[], tolerance=0.00005, cubic=True):
        self.sld=sld
        self.tolerance=tolerance
        self.sections=sections
        self.sect_xpos=np.array([sect.center[0] for sect in sections])
        if cubic:
            centerfunc_y=sinterp.CubicSpline(self.sect_xpos, [sect.center[1] for sect in sections])
            centerfunc_z=sinterp.CubicSpline(self.sect_xpos, [sect.center[2] for sect in sections])
        else:
            centerfunc_y=sinterp.interp1d(self.sect_xpos, [sect.center[1] for sect in sections])
            centerfunc_z=sinterp.interp1d(self.sect_xpos, [sect.center[2] for sect in sections])
        self.center=lambda x: np.vstack((x, centerfunc_y(x), centerfunc_z(x)))
        self.body_panels=[]
        self.last_x=sinterp.interp1d(self.sect_xpos, self.sect_xpos, kind='previous')
        self.next_x=sinterp.interp1d(self.sect_xpos, self.sect_xpos, kind='next')
        for i in range(len(self.sections)-1):
            for j in range(np.size(self.sections[i].coords, 0)-1):
                self.body_panels+=[body_panel(self.sections[i].coords[j, :], \
                    self.sections[i+1].coords[j, :], self.sections[i+1].coords[j+1, :], \
                    self.sections[i].coords[j+1, :], tolerance=tolerance)]
    def find_body_intersect(self, p, u, tolerance=0.00005):
        #return abutment point for wing called function to push out wing points
        return toolkit.get_panel_contact(p, u, np.array([p.Mtosys for p in self.body_panels]), \
            np.array([list(p.Mtouni) for p in self.body_panels]), np.array([list(p.locpoints) for p in self.body_panels]), \
                np.array([list(p.center) for p in self.body_panels]), tolerance)
    def plot_input(self, fig=None, ax=None, show=False, xlim=[], \
        ylim=[], zlim=[], colour='gray'): #plot input geometry data
        if fig==None:
            fig=plt.figure()
        if ax==None:
            plt.axes(projection='3d')
        
        order=np.array([0, 1, 2, 3, 0])
        for p in self.body_panels:
            ax.plot3D(p.points[0, order], p.points[1, order], p.points[2, order], colour)
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        if show:
            plt.show()
    def surfinterp(self, th, x):
        #function to generate linear interpolations along the section surface
        lastsect=0
        nextsect=1
        while nextsect<len(self.sections) and (not (self.sections[lastsect].center[0]<=x and x<=self.sections[nextsect].center[0])):
            lastsect+=1
            nextsect+=1
        eta=(x-self.sections[lastsect].center[0])/(self.sections[nextsect].center[0]-self.sections[lastsect].center[0])
        return (1.0-eta)*self.sections[lastsect](th)+eta*self.sections[nextsect](th)
    def side_separate(self, leftqueue=[], rightqueue=[], upqueue=[], lowqueue=[], xstrategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, xdisc=100):
        #function to define x position of points along solid lateral for wing and empenage abutments

        #distinguishing interval limits in lateral line
        #left wing side
        nabut=len(leftqueue)
        if nabut==0:
            lxlims=[[self.sections[0].center[0], self.sections[-1].center[0]]]
        else:
            lxlims=[[self.sections[0].center[0], leftqueue[0].extraright_points[-1][0]]]
        for i in range(1, nabut):
            lxlims+=[[leftqueue[i-1].extraright_points[0][0], leftqueue[i].extraright_points[-1][0]]]
        if nabut!=0:
            lxlims+=[[leftqueue[nabut-1].extraright_points[0][0], self.sections[-1].center[0]]]

        #right wing side
        nabut=len(rightqueue)
        if nabut==0:
            rxlims=[[self.sections[0].center[0], self.sections[-1].center[0]]]
        else:
            rxlims=[[self.sections[0].center[0], rightqueue[0].extraleft_points[-1][0]]]
        for i in range(1, nabut):
            rxlims+=[[rightqueue[i-1].extraleft_points[0][0], rightqueue[i].extraleft_points[-1][0]]]
        if nabut!=0:
            rxlims+=[[rightqueue[nabut-1].extraleft_points[0][0], self.sections[-1].center[0]]]

        #upper fin side
        nabut=len(upqueue)
        if nabut==0:
            uxlims=[[self.sections[0].center[0], self.sections[-1].center[0]]]
        else:
            uxlims=[[self.sections[0].center[0], upqueue[0].extraright_points[-1][0]]]
        for i in range(1, nabut):
            uxlims+=[[upqueue[i-1].extraright_points[0][0], upqueue[i].extraright_points[-1][0]]]
        if nabut!=0:
            uxlims+=[[upqueue[nabut-1].extraright_points[0][0], self.sections[-1].center[0]]]

        #lower fin side
        nabut=len(lowqueue)
        if nabut==0:
            dxlims=[[self.sections[0].center[0], self.sections[-1].center[0]]]
        else:
            dxlims=[[self.sections[0].center[0], lowqueue[0].extraleft_points[-1][0]]]
        for i in range(1, nabut):
            dxlims+=[[lowqueue[i-1].extraleft_points[0][0], lowqueue[i].extraleft_points[-1][0]]]
        if nabut!=0:
            dxlims+=[[lowqueue[nabut-1].extraleft_points[0][0], self.sections[-1].center[0]]]
        
        #deduce number of points based on segment length
        ltotal=(self.sections[-1].center[0]-self.sections[0].center[0])

        #left wing side
        lxdisc=[]
        for surf in leftqueue:
            lxdisc+=[surf.extraright_points[len(surf.extraright_points)-1-i][0] for i in range(len(surf.extraright_points)-1)]
        lnlines=[]
        for i in range(len(lxlims)-1):
            lnlines+=[floor(np.interp(lxlims[i][1]-lxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        lnlines+=[xdisc-len(lxdisc)-sum(lnlines)]
        for i in range(len(lxlims)-1):
            lxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, lnlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(lxlims[i])))
        lxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, lnlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(lxlims[-1])))
        lxdisc.sort()
        lxdisc.reverse()
        lnlines.reverse()
        lxdisc=np.array(lxdisc)

        #right wing side
        rxdisc=[]
        for surf in rightqueue:
            rxdisc+=[surf.extraleft_points[len(surf.extraleft_points)-1-i][0] for i in range(len(surf.extraleft_points)-1)]
        rnlines=[]
        for i in range(len(rxlims)-1):
            rnlines+=[floor(np.interp(rxlims[i][1]-rxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        rnlines+=[xdisc-len(rxdisc)-sum(rnlines)]
        for i in range(len(rxlims)-1):
            rxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, rnlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(rxlims[i])))
        rxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, rnlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(rxlims[-1])))
        rxdisc.sort()
        rxdisc.reverse()
        rnlines.reverse()
        rxdisc=np.array(rxdisc)

        #upper fin side
        uxdisc=[]
        for surf in upqueue:
            uxdisc+=[surf.extraright_points[len(surf.extraright_points)-1-i][0] for i in range(len(surf.extraright_points)-1)]
        unlines=[]
        for i in range(len(uxlims)-1):
            unlines+=[floor(np.interp(uxlims[i][1]-uxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        unlines+=[xdisc-len(uxdisc)-sum(unlines)]
        for i in range(len(uxlims)-1):
            uxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, unlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(uxlims[i])))
        uxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, unlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(uxlims[-1])))
        uxdisc.sort()
        uxdisc.reverse()
        unlines.reverse()
        uxdisc=np.array(uxdisc)

        #lower fin side
        dxdisc=[]
        for surf in lowqueue:
            dxdisc+=[surf.extraleft_points[len(surf.extraleft_points)-1-i][0] for i in range(len(surf.extraleft_points)-1)]
        dnlines=[]
        for i in range(len(dxlims)-1):
            dnlines+=[floor(np.interp(dxlims[i][1]-dxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        dnlines+=[xdisc-len(dxdisc)-sum(dnlines)]
        for i in range(len(dxlims)-1):
            dxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, dnlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(dxlims[i])))
        dxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, dnlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(dxlims[-1])))
        dxdisc.sort()
        dxdisc.reverse()
        dnlines.reverse()
        dxdisc=np.array(dxdisc)

        return lnlines, lxdisc, rnlines, rxdisc, unlines, uxdisc, dnlines, dxdisc
    def line_surfinterp(self, x1, x2, th1, th2, xstrategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, thstrategy=lambda x: x, disc=10):
        #generate line interpolation along fuselage surface

        #for cases referent to a negative position in the z axis, in which a position in the lower, left quadrant of a body has to
        #be abutted to the low, right quadrant of it:
        #(this chunk of code pressuposes all bodies will be abutted in the conventioned direction)
        if th1<0.0 and th2>0.0:
            thinterval=np.array([th1, th2-2*pi]) #(trim_polars()/trim_polars_array function in utils.py should handle inputs out
            #of [-pi; pi] interval)
        else:
            thinterval=np.array([th1, th2])
        xspacing=np.interp(xstrategy(np.linspace(0.0, 1.0, disc+1, endpoint=True)), np.array([0.0, 1.0]), np.array([x1, x2]))
        thspacing=np.interp(thstrategy(np.linspace(0.0, 1.0, disc+1, endpoint=True)), np.array([0.0, 1.0]), thinterval)
        thspacing=trim_polars_array(thspacing)
        pts=[]
        for i in range(disc+1):
            pts+=[self.surfinterp(thspacing[i], xspacing[i])]
        return pts
    def theta_queueident(self, queue, xspacing, intra=False, right=False, queueident='l'):
        #return yz argument of points along queue at hand
        #queueident: 'l', 'r', 'u' or 'd' for left, right, up or low (down) queues
        thetas=[]
        xposit=[]

        #bow
        xposit+=[xspacing[0]]
        if queueident=='l':
            thetas+=[-pi/2]
        elif queueident=='r':
            thetas+=[pi/2]
        elif queueident=='u':
            thetas+=[0.0]
        elif queueident=='d':
            thetas+=[pi]
        
        #queued wing quadrants
        qaux=queue
        qaux.reverse()
        for surf in qaux:
            if intra:
                if right:
                    surfpts=surf.intraright_points
                else:
                    surfpts=surf.intraleft_points
            else:
                if right:
                    surfpts=surf.extraright_points
                else:
                    surfpts=surf.extraleft_points
            xposit+=[pt[0] for pt in surfpts]
            thetas+=[toolkit.pointarg(self.center(pt[0]), pt) for pt in surfpts]

        #stern
        xposit+=[xspacing[-1]]
        if queueident=='l':
            thetas+=[-pi/2]
        elif queueident=='r':
            thetas+=[pi/2]
        elif queueident=='u':
            thetas+=[0.0]
        elif queueident=='d':
            thetas+=[pi]

        thx=sinterp.interp1d(xposit, thetas) #change this later to np.interp if you can

        return thx(xspacing)
    def patchcompose(self, leftqueue=[], rightqueue=[], upqueue=[], lowqueue=[], xstrategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, xdisc=100, \
        thstrategy=lambda x: x, thdisc_upleft=20, thdisc_upright=20, thdisc_downleft=20, thdisc_downright=20):
        #function to generate patch to add panels belonging to fuselage. Queues designate abutted surfaces at a certain point
        #in the surface (e. g. upqueue to vertical fin, leftqueue to left wing (wing 1)...)
        #must be ran AFTER abutted wings's 'patchcompose's.
        lnlines, lxdisc, rnlines, rxdisc, unlines, uxdisc, dnlines, dxdisc=self.side_separate(leftqueue=leftqueue, \
            rightqueue=rightqueue, upqueue=upqueue, lowqueue=lowqueue, xstrategy=xstrategy, xdisc=xdisc)

        #left-side, 0.0 to -pi/2 patch
        thleft=self.theta_queueident(leftqueue, xspacing=lxdisc, intra=False, queueident='l', right=True)
        thright=self.theta_queueident(upqueue, xspacing=uxdisc, intra=True, right=True, queueident='u')
        ptpatch=[]
        for i in range(len(lxdisc)):
            ptpatch+=[self.line_surfinterp(uxdisc[i], lxdisc[i], thright[i], thleft[i], xstrategy=xstrategy, thstrategy=thstrategy, disc=thdisc_upleft)]
        horzlines, vertlines, paninds, sld=self.sld.addpatch(ptpatch, prevlines={'left':prevline_organize(leftqueue, lnlines, prevlateral=[], intra=False, right=True), \
            'right':prevline_organize(upqueue, unlines, prevlateral=[], intra=True, right=True)})
        uplat=[linerow[0] for linerow in vertlines]
        
        #left side, -pi/2 to -pi patch
        thleft=self.theta_queueident(lowqueue, xspacing=dxdisc, intra=True, queueident='d', right=False)
        thright=self.theta_queueident(leftqueue, xspacing=lxdisc, intra=True, right=True, queueident='l')
        ptpatch=[]
        for i in range(len(lxdisc)):
            ptpatch+=[self.line_surfinterp(lxdisc[i], dxdisc[i], thright[i], thleft[i], xstrategy=xstrategy, thstrategy=thstrategy, disc=thdisc_downleft)]
        horzlines, vertlines, paninds, sld=self.sld.addpatch(ptpatch, prevlines={'left':prevline_organize(lowqueue, dnlines, prevlateral=[], intra=True, right=False), \
            'right':prevline_organize(leftqueue, lnlines, prevlateral=[linerow[-1] for linerow in vertlines], intra=True, right=True)})
        
        #right side, pi to pi/2 patch
        thleft=self.theta_queueident(rightqueue, xspacing=rxdisc, intra=True, queueident='r', right=False)
        thright=self.theta_queueident(lowqueue, xspacing=dxdisc, intra=False, right=False, queueident='d')
        ptpatch=[]
        for i in range(len(lxdisc)):
            ptpatch+=[self.line_surfinterp(dxdisc[i], rxdisc[i], thright[i], thleft[i], xstrategy=xstrategy, thstrategy=thstrategy, disc=thdisc_downright)]
        horzlines, vertlines, paninds, sld=self.sld.addpatch(ptpatch, prevlines={'left':prevline_organize(rightqueue, rnlines, prevlateral=[], intra=True, right=False), \
            'right':prevline_organize(lowqueue, dnlines, prevlateral=[linerow[-1] for linerow in vertlines], intra=False, right=False)})
        
        #right side, pi/2 to 0.0 patch
        thleft=self.theta_queueident(upqueue, xspacing=uxdisc, intra=False, queueident='u', right=True)
        thright=self.theta_queueident(rightqueue, xspacing=rxdisc, intra=False, right=False, queueident='r')
        ptpatch=[]
        for i in range(len(lxdisc)):
            ptpatch+=[self.line_surfinterp(rxdisc[i], uxdisc[i], thright[i], thleft[i], xstrategy=xstrategy, thstrategy=thstrategy, disc=thdisc_upright)]
        horzlines, vertlines, paninds, sld=self.sld.addpatch(ptpatch, prevlines={'left':prevline_organize(upqueue, unlines, prevlateral=uplat, intra=False, right=True), \
            'right':prevline_organize(rightqueue, rnlines, prevlateral=[linerow[-1] for linerow in vertlines], intra=False, right=False)})

def circdefsect(R=1.0, center=np.array([0.0, 0.0, 0.0]), cubic=True, disc=360): #generate circular defsect based on 
    #discretization necessities
    return body_section(coords=gen_circdefsect_coords(disc), R=R, center=center, cubic=cubic)

def standard_body(sld, nose_loc=np.array([0.0, 0.0, 0.0]), head_length=0.1, head_thdisc=10, body_length=1.0, \
    body_width=0.1, tailcone_length=0.2, body_thdisc=60, cubic=False, tolerance=0.00005):
    #standard body: semi-ellipsoid head, straight cylindric body and conventional tailcone, all along x axis
    sects=[]
    head_thetas=np.linspace(0.0, pi/2, head_thdisc, endpoint=False)
    Rs=[]
    xs=[]

    for th in head_thetas:
        xs+=[(1-cos(th))*head_length]
        Rs+=[sin(th)*body_width/2]

    xs+=[head_length]
    Rs+=[body_width/2]

    xs+=[body_length-tailcone_length]
    Rs+=[body_width/2]

    xs+=[body_length]
    Rs+=[0]

    for i in range(len(xs)):
        sects+=[circdefsect(R=Rs[i], center=nose_loc+np.array([xs[i], 0.0, 0.0]), cubic=cubic, disc=body_thdisc)]
    return body(sld, sections=sects, tolerance=tolerance, cubic=cubic)