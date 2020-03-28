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

#turbulence criterion for flat plate equivalence in parasite drag estimation
def Re2e5(Re):
    return Re>2e5

#Blausius's solution laminar friction coefficient
def Blausius_Cf_l(Re):
    return 0.664/sqrt(Re)

#Prandtl's one-seventh power law for tubulent boundary layer friction coefficient
def Prandtl_1_7th(Re):
    return 0.027/(Re**(1.0/7))

#class to define functions and variables needed to represent body sections
class body_section:
    def __init__(self, center=np.array([0.0, 0.0, 0.0]), coords=np.vstack((np.sin(np.linspace(0.0, 2*pi, 360)), \
        np.cos(np.linspace(0.0, 2*pi, 360)))).T, R=1.0, y_expand=1.0, z_expand=1.0, cubic=True):
        self.R=R
        self.coords=np.vstack((np.zeros(np.size(coords, 0)), coords[:, 0], coords[:, 1])).T*R #redefine coordinate system
        self.center=center
        thetas=np.arctan2(coords[:, 0], coords[:, 1])
        if thetas[0]>0.0: #detect positive argument deduction for (0, -1)
            thetas[0]*=-1
        Rs=lg.norm(coords, axis=1)
        order=np.argsort(thetas)
        self.thetas=thetas[order]
        self.Rs=Rs[order]
        if cubic:
            self.polar_rule=sinterp.CubicSpline(self.thetas, self.Rs, extrapolate=True)
        else:
            self.polar_rule=sinterp.interp1d(self.thetas, self.Rs)
        self.y_expand=y_expand; self.z_expand=z_expand
        self.coords[:, 1]*=y_expand
        self.coords[:, 2]*=z_expand
        for i in range(np.size(coords, 0)):
            self.coords[i, :]+=center
        #define function for interpolation in polar coordinates between theta in respect to z axis as input
    def __call__(self, th):
        R=self.polar_rule(th)
        return np.array([0.0, self.y_expand*sin(th), self.z_expand*cos(th)])*R*self.R+self.center

class body_panel: #exclusively for intersection finding purposes
    def __init__(self, p1, p2, p3, p4, tolerance=0.00005):
        points=np.vstack((p1, p2, p3, p4)).T
        self.points=points
        self.center, self.Mtosys, self.Mtouni, self.locpoints, error=toolkit.body_panel_process(self.points, tolerance)
        if error:
            print('WARNING: error in body panel generation. Please check input geometry\'s integrity')
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
                linelist+=queue[-i-1].intraright_lines
            else:
                linelist+=queue[-i-1].intraleft_lines
        else:
            if right:
                linelist+=queue[-i-1].extraright_lines
            else:
                linelist+=queue[-i-1].extraleft_lines
    if len(prevlateral)==0:
        linelist+=nlines[-1]*[-2]
    else:
        linelist+=[prevlateral[i] for i in range(len(linelist), nlines[-1]+len(linelist))]
    return linelist

class body: #body definition class for center definition to obtain polar cooridinates and creation of input panels 
    #for abutments.
    #to be abutted stern to bow, clockwise when seen from the aircraft's bow
    def set_aircraft(self, acft):
        self.acft=acft #define aircraft structure related to instance
    def __init__(self, sld, sections=[], tolerance=0.00005):
        self.sld=sld
        self.tolerance=tolerance
        self.sections=sections
        self.sect_xpos=np.array([sect.center[0] for sect in sections])
        centerfunc_y=sinterp.interp1d(self.sect_xpos, [sect.center[1] for sect in sections])
        centerfunc_z=sinterp.interp1d(self.sect_xpos, [sect.center[2] for sect in sections])
        self.center=lambda x: np.vstack((x, centerfunc_y(x), centerfunc_z(x)))
        self.body_panels=[]
        self.last_x=sinterp.interp1d(self.sect_xpos, self.sect_xpos, kind='previous')
        self.next_x=sinterp.interp1d(self.sect_xpos, self.sect_xpos, kind='next')
        self.y_expand_rule=sinterp.interp1d(self.sect_xpos, np.array([sect.y_expand for sect in self.sections]))
        self.z_expand_rule=sinterp.interp1d(self.sect_xpos, np.array([sect.z_expand for sect in self.sections]))
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
        lnlines=[]
        for i in range(len(lxlims)-1):
            lnlines+=[floor(np.interp(lxlims[i][1]-lxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        lnlines+=[xdisc-sum([len(surf.extraright_lines) for surf in leftqueue])-sum(lnlines)]
        lnlines.reverse()
        lxlims.reverse()
        for sind in range(len(leftqueue)-1, -1, -1):
            i=len(leftqueue)-1-sind
            lxlims[i].reverse()
            surf=leftqueue[sind]
            lxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, lnlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(lxlims[i])))
            lxdisc+=[surf.extraright_points[i][0] for i in range(len(surf.extraright_points)-1)]
        lxlims[-1].reverse()
        lxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, lnlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(lxlims[-1])))
        lxdisc=np.array(lxdisc)

        #right wing side
        rxdisc=[]
        rnlines=[]
        for i in range(len(rxlims)-1):
            rnlines+=[floor(np.interp(rxlims[i][1]-rxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        rnlines+=[xdisc-sum([len(surf.extraleft_lines) for surf in rightqueue])-sum(rnlines)]
        rnlines.reverse()
        rxlims.reverse()
        for sind in range(len(rightqueue)-1, -1, -1):
            i=len(rightqueue)-1-sind
            rxlims[i].reverse()
            surf=rightqueue[sind]
            rxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, rnlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(rxlims[i])))
            rxdisc+=[surf.extraleft_points[i][0] for i in range(len(surf.extraleft_points)-1)]
        rxlims[-1].reverse()
        rxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, rnlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(rxlims[-1])))
        rxdisc=np.array(rxdisc)

        #upper fin side
        uxdisc=[]
        unlines=[]
        for i in range(len(uxlims)-1):
            unlines+=[floor(np.interp(uxlims[i][1]-uxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        unlines+=[xdisc-sum([len(surf.extraright_lines) for surf in upqueue])-sum(unlines)]
        unlines.reverse()
        uxlims.reverse()
        for sind in range(len(upqueue)-1, -1, -1):
            i=len(upqueue)-1-sind
            uxlims[i].reverse()
            surf=upqueue[sind]
            uxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, unlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(uxlims[i])))
            uxdisc+=[surf.extraright_points[i][0] for i in range(len(surf.extraright_points)-1)]
        uxlims[-1].reverse()
        uxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, unlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(uxlims[-1])))
        uxdisc=np.array(uxdisc)

        #lower fin side
        dxdisc=[]
        dnlines=[]
        for i in range(len(dxlims)-1):
            dnlines+=[floor(np.interp(dxlims[i][1]-dxlims[i][0], np.array([0.0, ltotal]), np.array([0, xdisc])))]
        dnlines+=[xdisc-sum([len(surf.extraleft_lines) for surf in lowqueue])-sum(dnlines)]
        dnlines.reverse()
        dxlims.reverse()
        for sind in range(len(lowqueue)-1, -1, -1):
            i=len(lowqueue)-1-sind
            dxlims[i].reverse()
            surf=lowqueue[sind]
            dxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, dnlines[i], endpoint=False)), np.array([0.0, 1.0]), np.array(dxlims[i])))
            dxdisc+=[surf.extraleft_points[i][0] for i in range(len(surf.extraleft_points)-1)]
        dxlims[-1].reverse()
        dxdisc+=list(np.interp(xstrategy(np.linspace(0.0, 1.0, dnlines[-1]+1, endpoint=True)), np.array([0.0, 1.0]), np.array(dxlims[-1])))
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
        qaux=[queue[-1-i] for i in range(len(queue))]
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
            #thetas+=[toolkit.pointarg(self.center(pt[0]), pt) for pt in surfpts]
            thetas+=[atan2((pt[1]-self.center(pt[0])[1])/self.y_expand_rule(pt[0]), (pt[2]-self.center(pt[0])[2])/self.z_expand_rule(pt[0])) \
                for pt in surfpts]

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
    def patchcompose(self, leftqueue=[], rightqueue=[], upqueue=[], lowqueue=[], xstrategy=lambda x: x, xdisc=100, \
        thstrategy=lambda x: x, thdisc_upleft=20, thdisc_upright=20, thdisc_downleft=20, thdisc_downright=20, tolerance=5e-5):
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
            'right':prevline_organize(upqueue, unlines, prevlateral=[], intra=True, right=True)}, tolerance=tolerance)
        uplat=[linerow[0] for linerow in vertlines]
        self.paninds=[]
        for panlist in paninds:
            self.paninds+=panlist
        
        #left side, -pi/2 to -pi patch
        thleft=self.theta_queueident(lowqueue, xspacing=dxdisc, intra=True, queueident='d', right=False)
        thright=self.theta_queueident(leftqueue, xspacing=lxdisc, intra=True, right=True, queueident='l')
        ptpatch=[]
        for i in range(len(lxdisc)):
            ptpatch+=[self.line_surfinterp(lxdisc[i], dxdisc[i], thright[i], thleft[i], xstrategy=xstrategy, thstrategy=thstrategy, disc=thdisc_downleft)]
        horzlines, vertlines, paninds, sld=self.sld.addpatch(ptpatch, prevlines={'left':prevline_organize(lowqueue, dnlines, prevlateral=[], intra=True, right=False), \
            'right':prevline_organize(leftqueue, lnlines, prevlateral=[linerow[-1] for linerow in vertlines], intra=True, right=True)}, tolerance=tolerance)
        for panlist in paninds:
            self.paninds+=panlist
        
        #right side, pi to pi/2 patch
        thleft=self.theta_queueident(rightqueue, xspacing=rxdisc, intra=True, queueident='r', right=False)
        thright=self.theta_queueident(lowqueue, xspacing=dxdisc, intra=False, right=False, queueident='d')
        ptpatch=[]
        for i in range(len(lxdisc)):
            ptpatch+=[self.line_surfinterp(dxdisc[i], rxdisc[i], thright[i], thleft[i], xstrategy=xstrategy, thstrategy=thstrategy, disc=thdisc_downright)]
        horzlines, vertlines, paninds, sld=self.sld.addpatch(ptpatch, prevlines={'left':prevline_organize(rightqueue, rnlines, prevlateral=[], intra=True, right=False), \
            'right':prevline_organize(lowqueue, dnlines, prevlateral=[linerow[-1] for linerow in vertlines], intra=False, right=False)}, tolerance=tolerance)
        for panlist in paninds:
            self.paninds+=panlist
        
        #right side, pi/2 to 0.0 patch
        thleft=self.theta_queueident(upqueue, xspacing=uxdisc, intra=False, queueident='u', right=True)
        thright=self.theta_queueident(rightqueue, xspacing=rxdisc, intra=False, right=False, queueident='r')
        ptpatch=[]
        for i in range(len(lxdisc)):
            ptpatch+=[self.line_surfinterp(rxdisc[i], uxdisc[i], thright[i], thleft[i], xstrategy=xstrategy, thstrategy=thstrategy, disc=thdisc_upright)]
        horzlines, vertlines, paninds, sld=self.sld.addpatch(ptpatch, prevlines={'left':prevline_organize(upqueue, unlines, prevlateral=uplat, intra=False, right=True), \
            'right':prevline_organize(rightqueue, rnlines, prevlateral=[linerow[-1] for linerow in vertlines], intra=False, right=False)}, tolerance=tolerance)
        for panlist in paninds:
            self.paninds+=panlist
    def apply_eqflatplate(self, rho=1.225, Uinf=1.0, mu=1.72*10e-5, turbulent_criterion=Re2e5, Cf_l_rule=Blausius_Cf_l, Cf_t_rule=Prandtl_1_7th):
        #all rule parameters should be provided as a function of the local Reynolds number
        for p in self.paninds:
            local_Re=((self.sld.panels[p].colpoint[0]-self.sect_xpos[0])*rho*Uinf)/mu
            if turbulent_criterion(local_Re):
                self.sld.Cfs[p]=Cf_t_rule(local_Re)
            else:
                self.sld.Cfs[p]=Cf_l_rule(local_Re)
    def bodypanel_plotnormals(self, xlim=[], ylim=[], zlim=[], factor=1.0):
        #plot normal vectors of body panels. Function essentially produced for geometry gen. debugging
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        for i in range(len(self.body_panels)):
            for j in range(4):
                ax.plot3D([self.body_panels[i].points[0, j], self.body_panels[i].points[0, (j+1)%4]], \
                    [self.body_panels[i].points[1, j], self.body_panels[i].points[1, (j+1)%4]], \
                        [self.body_panels[i].points[2, j], self.body_panels[i].points[2, (j+1)%4]], 'gray')
            ax.quiver(self.body_panels[i].center[0], self.body_panels[i].center[1], self.body_panels[i].center[2], \
                self.body_panels[i].Mtouni[0, 2]*factor, self.body_panels[i].Mtouni[1, 2]*factor, self.body_panels[i].Mtouni[2, 2]*factor)
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

def circdefsect(y_expand=1.0, z_expand=1.0, R=1.0, center=np.array([0.0, 0.0, 0.0]), cubic=True, disc=360): #generate circular defsect based on 
    #discretization necessities
    return body_section(coords=gen_circdefsect_coords(disc), R=R, y_expand=y_expand, z_expand=z_expand, center=center, cubic=cubic)

def squaredefsect(R=1.0, center=np.array([0.0, 0.0, 0.0]), cubic=True, disc=360, z_expand=1.0, y_expand=1.0): #generate circular defsect based on 
    #discretization necessities
    return body_section(coords=gen_squaredefsect_coords(disc), y_expand=y_expand, z_expand=z_expand, R=R, center=center, cubic=cubic)

def smooth_angle_defsect_function(r_1x=0.5, r_2x=0.5, r_1y=0.5, r_2y=0.5, ldisc=30, thdisc=20): #COMPLETE ME LATER
    #defsect with elliptic concordances on edges
    return lambda R=1.0, center=np.array([0.0, 0.0, 0.0]), y_expand=1.0, z_expand=1.0: \
        body_section(R=R, center=center, z_expand=z_expand, y_expand=y_expand, coords=\
            smooth_angle_defsect_coords(r_1x=r_1x, r_1y=r_1y, r_2x=r_2x, r_2y=r_2y, ldisc=ldisc, thdisc=thdisc))

def standard_body(sld, defsect=circdefsect, nose_loc=np.array([0.0, 0.0, 0.0]), nose_length=0.1, nose_thdisc=10, body_length=1.0, \
    body_width=0.1, tailcone_length=0.2, body_thdisc=60, tolerance=0.00005, nose_lift=0.0, tail_lift=0.0, z_expand=1.0, \
y_expand=1.0):
    #standard body: semi-ellipsoid head, straight cylindric body and conventional tailcone, all along x axis
    sects=[]
    head_thetas=np.linspace(0.0, pi/2, nose_thdisc, endpoint=False)
    Rs=[]
    xs=[]

    for th in head_thetas:
        xs+=[(1-cos(th))*nose_length]
        Rs+=[sin(th)*body_width/2]

    xs+=[nose_length]
    Rs+=[body_width/2]

    xs+=[body_length-tailcone_length]
    Rs+=[body_width/2]

    xs+=[body_length]
    Rs+=[0]

    zs=[nose_lift, 0.0, 0.0, tail_lift]
    zs=np.interp(np.array(xs), np.array([0.0, nose_length, body_length-tailcone_length, body_length]), np.array(zs))

    for i in range(len(xs)):
        sects+=[defsect(z_expand=z_expand, y_expand=y_expand, R=Rs[i], center=nose_loc+np.array([xs[i], 0.0, zs[i]]), disc=body_thdisc)]
    return body(sld, sections=sects, tolerance=tolerance)