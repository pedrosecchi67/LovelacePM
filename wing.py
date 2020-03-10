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
from body import *

class wing_section: #class to define wing section based on airfoil info
    def __init__(self, c=1.0, incidence=0.0, gamma=0.0, CA_position=np.array([0.0, 0.0, 0.0]), afl='n4412', \
        header_lines=1, xdisc=10, remove_TE_gap=True, inverse=False, xstrategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, closed=False):
        self.CA_position=CA_position
        self.c=c
        self.closed=closed
        self.points=wing_afl_positprocess(read_afl(afl=afl, ext_append=True, header_lines=header_lines, disc=xdisc, \
            remove_TE_gap=remove_TE_gap, incidence=incidence, inverse=inverse, strategy=xstrategy, closed=closed), gamma=gamma, c=c, xpos=CA_position[0], \
                ypos=CA_position[1], zpos=CA_position[2])
        self.inverted=inverse

class wing_quadrant: #wing region between two airfoil sections
    def __init__(self, sld, sect1=None, sect2=None):
        #by convention, set section 1 as the rightmost section in the quadrant (positive in y axis)
        #vertical fins should have sect1 as their upmost section
        self.sld=sld
        self.sect1=sect1
        self.sect2=sect2
        self.wakecomb=[]
        self.panstrips_extra=[]
        self.panstrips_intra=[]
    def trim_bybody(self, contactbody, sectside=2, tolerance=0.00005):
        #trim wing section by body contact, for abutment
        if sectside==2:
            ps=self.sect1.points
            us=self.sect2.points-self.sect1.points
            for i in range(np.size(self.sect2.points, 0)):
                newpt, errorcode=contactbody.find_body_intersect(ps[i, :], us[i, :], tolerance=tolerance)
                self.sect2.points[i, :]=newpt
                if errorcode:
                    print('An error has been detected while performing abutments. Please check geometry.')
        else:
            ps=self.sect2.points
            us=self.sect1.points-self.sect2.points
            for i in range(np.size(self.sect2.points, 0)):
                newpt, errorcode=contactbody.find_body_intersect(ps[i, :], us[i, :], tolerance=tolerance)
                self.sect1.points[i, :]=newpt
                if errorcode:
                    print('An error has been detected while performing abutments. Please check geometry.')
    def plot_input(self, fig=None, ax=None, show=False, xlim=[], \
        ylim=[], zlim=[], colour='blue'): #plot input geometry data
        if fig==None:
            fig=plt.figure()
        if ax==None:
            plt.axes(projection='3d')
        
        ax.plot3D(self.sect1.points[:, 0], self.sect1.points[:, 1], self.sect1.points[:, 2], colour)
        ax.plot3D(self.sect2.points[:, 0], self.sect2.points[:, 1], self.sect2.points[:, 2], colour)
        
        for i in range(np.size(self.sect1.points, 0)):
            ax.plot3D([self.sect1.points[i, 0], self.sect2.points[i, 0]], \
                [self.sect1.points[i, 1], self.sect2.points[i, 1]], \
                    [self.sect1.points[i, 2], self.sect2.points[i, 2]], colour)
        
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        if show:
            plt.show()
    def patchcompose(self, prevlines={}, strategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, ldisc=20):
        #code for prevlines: 'intra_left', 'intra_right', 'extra_...
        #and wing quadrant patch additions
        lspacing=strategy(np.linspace(0.0, 1.0, ldisc))
        extrasld=[]
        xdisc=int(np.size(self.sect1.points, 0)/2)+1

        #detecting closed extremities
        closedl=self.sect1.closed
        closedr=self.sect2.closed

        #handling inverted airfoil sections
        if self.sect1.inverted:
            sect1_extraorder=range(np.size(self.sect1.points, 0)-1, xdisc-2, -1)
            sect1_intraorder=range(xdisc)
        else:
            sect1_intraorder=range(np.size(self.sect1.points, 0)-1, xdisc-2, -1)
            sect1_extraorder=range(xdisc)
        if self.sect2.inverted:
            sect2_extraorder=range(np.size(self.sect2.points, 0)-1, xdisc-2, -1)
            sect2_intraorder=range(xdisc)
        else:
            sect2_intraorder=range(np.size(self.sect2.points, 0)-1, xdisc-2, -1)
            sect2_extraorder=range(xdisc)

        for i in range(xdisc):
            extrasld+=[[]]
            for eta in lspacing:
                extrasld[-1]+=[eta*self.sect1.points[sect1_extraorder[i], :]+(1.0-eta)*self.sect2.points[sect2_extraorder[i], :]]
        intrasld=[]
        for i in range(xdisc):
            intrasld+=[[]]
            for eta in lspacing:
                intrasld[-1]+=[eta*self.sect2.points[sect2_intraorder[i], :]+(1.0-eta)*self.sect1.points[sect1_intraorder[i], :]]
        
        #wake info
        self.wakecombs=[]
        
        prev={}
        if 'extra_right' in prevlines:
            prev['right']=prevlines['extra_right']
        if 'extra_left' in prevlines:
            prev['left']=prevlines['extra_left']
        horzlines, vertlines, paninds, points=self.sld.addpatch(extrasld, prevlines=prev)
        self.panstrips_extra=[[paninds[i][j] for i in range(len(paninds))] for j in range(len(paninds[0]))]
        self.paninds=[]
        for panlist in paninds:
            self.paninds+=panlist
        TE_extra=paninds[0] #for wake generation
        self.extraright_lines=[vertlines[i][0] for i in range(len(vertlines))]
        self.extraright_points=[points[i][0] for i in range(len(points))]
        self.extraleft_lines=[vertlines[i][-1] for i in range(len(vertlines))]
        self.extraleft_points=[points[i][-1] for i in range(len(points))]
        LE_lines=[horzlines[-1][i] for i in range(len(horzlines[-1]))]
        LE_lines.reverse()
        TE_lines=[horzlines[0][i] for i in range(len(horzlines[0]))]
        TE_lines.reverse()
        prev={'up':LE_lines, 'low':TE_lines}
        if 'intra_right' in prevlines:
            prev['left']=prevlines['intra_right']
        if 'intra_left' in prevlines:
            prev['right']=prevlines['intra_left']
        if closedl:
            prev['right']=self.extraleft_lines
        if closedr:
            prev['left']=self.extraright_lines
        horzlines, vertlines, paninds, points=self.sld.addpatch(intrasld, prevlines=prev, invlats=['up', 'low'])
        self.panstrips_intra=[[paninds[i][j] for i in range(len(paninds))] for j in range(len(paninds[0])-1, -1, -1)]
        for panlist in paninds:
            self.paninds+=panlist
        TE_intra=paninds[0]
        TE_intra.reverse() #for wake generation
        self.intraright_lines=[vertlines[i][-1] for i in range(len(vertlines))]
        self.intraright_points=[points[i][-1] for i in range(len(points))]
        self.intraleft_lines=[vertlines[i][0] for i in range(len(vertlines))]
        self.intraleft_points=[points[i][0] for i in range(len(points))]

        #wake info
        for i in range(len(TE_extra)):
            self.wakecombs+=[[TE_extra[i], TE_intra[i]]]
    def calc_reference(self, axis=1): #input for wing's calc reference function
        ys=np.array([self.sect1.CA_position[axis], self.sect2.CA_position[axis]])
        cs=np.array([self.sect1.c, self.sect2.c])
        yspacing=np.linspace(ys[0], ys[1], 100)
        cspacing=np.interp(yspacing, ys, cs)
        #returns S and parcel for MAC integration
        return abs(np.trapz(cs, x=ys)), abs(np.trapz(cspacing**2, x=yspacing))
    def close_tip(self, sectside=2):
        if sectside==2:
            extra_paninds=self.panstrips_extra[0]
            for i in range(len(extra_paninds)):
                self.sld.panels[extra_paninds[i]].nocirc_enforce(self.extraright_lines[i])
            intra_paninds=self.panstrips_intra[0]
            for i in range(len(intra_paninds)):
                self.sld.panels[intra_paninds[i]].nocirc_enforce(self.intraright_lines[i])
        else:
            extra_paninds=self.panstrips_extra[-1]
            for i in range(len(extra_paninds)):
                self.sld.panels[extra_paninds[i]].nocirc_enforce(self.extraleft_lines[i])
            intra_paninds=self.panstrips_intra[-1]
            for i in range(len(intra_paninds)):
                self.sld.panels[intra_paninds[i]].nocirc_enforce(self.intraleft_lines[i])
    def calc_coefs(self, alpha=0.0, axis=1): #calculate local sectional coefficients
        if axis==1: #calculate unitary, streamwise direction vectors
            u=np.array([cos(alpha), 0.0, sin(alpha)])
            v=np.array([-sin(alpha), 0.0, cos(alpha)])
            ax=np.array([0.0, 1.0, 0.0])
            perpaxis=2
        elif axis==2:
            u=np.array([cos(alpha), sin(alpha), 0.0])
            v=np.array([-sin(alpha), cos(alpha), 0.0])
            ax=np.array([0.0, 0.0, 1.0])
            perpaxis=1
        else:
            print('WARNING: error in axis provided to wingquad.calc_coefs()')
        ys=[]
        cs=[]
        CA_posits=[]
        Cls=[]
        Cms=[]
        Cds=[]
        Gammas=[]
        for i in range(len(self.panstrips_extra)):
            xposits=np.array([self.sld.panels[pind].colpoint[0] for pind in self.panstrips_extra[i]])
            yposits=np.array([self.sld.panels[pind].colpoint[axis] for pind in self.panstrips_extra[i]])
            aftextreme=min(xposits)
            rearextreme=max(xposits)
            cs+=[rearextreme-aftextreme]
            ys+=[np.mean(yposits)]
            CA_posits+=[np.array([rearextreme/4+3*aftextreme/4, ys[-1], \
                np.interp(ys[-1], np.array([self.sect1.CA_position[1], self.sect2.CA_position[1]]), \
                    np.array([self.sect1.CA_position[2], self.sect2.CA_position[2]]))])]
            extracps=np.array([self.sld.Cps[pind] for pind in self.panstrips_extra[i]])
            intracps=np.array([self.sld.Cps[pind] for pind in self.panstrips_intra[i]])
            extrax=xposits-CA_posits[-1][0]
            perposits=np.array([self.sld.panels[pind].colpoint[perpaxis] for pind in self.panstrips_intra[i]])
            perposits-=np.mean(perposits)
            intrax=np.array([self.sld.panels[pind].colpoint[0] for pind in self.panstrips_intra[i]])-CA_posits[-1][0]
            extranu=np.array([self.sld.panels[pind].nvector@u for pind in self.panstrips_extra[i]])
            extranv=np.array([self.sld.panels[pind].nvector@v for pind in self.panstrips_extra[i]])
            intranu=np.array([self.sld.panels[pind].nvector@u for pind in self.panstrips_intra[i]])
            intranv=np.array([self.sld.panels[pind].nvector@v for pind in self.panstrips_intra[i]])
            Cls+=[(np.trapz(extracps*extranv, x=extrax)+np.trapz(intracps*intranv, x=intrax))/cs[-1]]
            circulation=0.0
            gamma_loc=np.zeros(len(self.panstrips_extra[i]))
            for j in range(len(gamma_loc)):
                vloc=self.sld.delphi[self.panstrips_extra[i][j], :]+self.sld.vbar[self.panstrips_extra[i][j], :]
                vloc-=ax*(vloc@ax)
                gamma_loc[j]=lg.norm(vloc)
            circulation-=np.trapz(gamma_loc*extranv, x=extrax)
            for j in range(len(gamma_loc)):
                vloc=self.sld.delphi[self.panstrips_intra[i][j], :]+self.sld.vbar[self.panstrips_intra[i][j], :]
                vloc-=ax*(vloc@ax)
                gamma_loc[j]=lg.norm(vloc)
            circulation-=np.trapz(gamma_loc*intranv, x=intrax)
            Gammas+=[circulation]
            Cds+=[(np.trapz(extracps*extranu, x=extrax)+np.trapz(intracps*intranu, x=intrax))/cs[-1]]
            Cms+=[-(np.trapz(extracps*extranv*extrax, x=extrax)+np.trapz(intracps*intranv*intrax, x=intrax))/cs[-1]**2]
        Cls=np.array(Cls)
        Cds=np.array(Cds)
        Cms=np.array(Cms)
        ys=np.array(ys)
        cs=np.array(cs)
        Gammas=np.array(Gammas)
        '''Cls=np.gradient(Cls, ys)
        Cms=np.gradient(Cms, ys)
        Cds=np.gradient(Cds, ys)'''

        return ys, cs, Cls, Cds, Cms, Gammas


class wing:
    def __init__(self, wingquads=[]):
        self.coefavailable=False
        self.wingquads=wingquads
        self.leftclosed=False
        self.rightclosed=False
    def patchcompose(self, ystrategy=lambda x: x, ydisc=20):
        if type(ydisc)==list: #use as list if list is provided, specifying each of the quadrants
            trimlist(len(self.wingquads), ydisc)
        else: #distribute between quadrants proportionally to distance between quadrant section CAs if integer
            ltemp=[]
            for quad in self.wingquads:
                ltemp+=[lg.norm(quad.sect1.CA_position-quad.sect2.CA_position)]
            ltot=sum(ltemp)
            disctemp=[]
            for i in range(len(self.wingquads)-1):
                disctemp+=[floor(np.interp(ltemp[i]/ltot, np.array([0.0, 1.0]), np.array([0, ydisc])))]
            disctemp+=[ydisc-sum(disctemp)]
            ydisc=disctemp
        
        #generate patches one by one
        
        #first patch
        self.wingquads[0].patchcompose(strategy=ystrategy, ldisc=ydisc[0])
        #all other patches
        for i in range(1, len(self.wingquads)):
            self.wingquads[i].patchcompose(prevlines={'intra_left':self.wingquads[i-1].intraright_lines, 'extra_left':self.wingquads[i-1].extraright_lines}, \
                strategy=ystrategy, ldisc=ydisc[i])

        self.extraleft_lines=self.wingquads[0].extraleft_lines
        self.extraright_lines=self.wingquads[-1].extraright_lines
        self.intraleft_lines=self.wingquads[0].intraleft_lines
        self.intraright_lines=self.wingquads[-1].intraright_lines
        self.extraleft_points=self.wingquads[0].extraleft_points
        self.extraright_points=self.wingquads[-1].extraright_points
        self.intraleft_points=self.wingquads[0].intraleft_points
        self.intraright_points=self.wingquads[-1].intraright_points

        self.paninds=[]
        for quad in self.wingquads:
            self.paninds+=quad.paninds
    def trim_bybody(self, contactbody, sectside=2, tolerance=0.00005):
        #trim wing section by body contact, for abutment
        if sectside==2:
            self.wingquads[-1].trim_bybody(contactbody, sectside=2, tolerance=tolerance)
        else:
            self.wingquads[0].trim_bybody(contactbody, sectside=1, tolerance=tolerance)
    def calc_coefs(self, alpha=0.0, axis=1):
        self.ys=np.array([])
        self.cs=np.array([])
        self.Cls=np.array([])
        self.Cds=np.array([])
        self.Cms=np.array([])
        self.Gammas=np.array([])
        for quad in self.wingquads:
            y, c, Cl, Cd, Cm, Gamma=quad.calc_coefs(alpha=alpha, axis=axis)
            self.ys=np.hstack((self.ys, y))
            self.Cls=np.hstack((self.Cls, Cl))
            self.Cds=np.hstack((self.Cds, Cd))
            self.Cms=np.hstack((self.Cms, Cm))
            self.cs=np.hstack((self.cs, c))
            self.Gammas=np.hstack((self.Gammas, Gamma))
        order=np.argsort(self.ys)
        self.ys=self.ys[order]
        self.cs=self.cs[order]
        self.Cls=self.Cls[order]
        self.Cds=self.Cds[order]
        self.Cms=self.Cms[order]
        self.Gammas=self.Gammas[order]
        self.coefavailable=True
    def calc_reference(self): #calculate local MAC, surface and span
        #identify wing axis (1 or 2)
        disty=abs(self.wingquads[0].sect1.CA_position[1]-self.wingquads[-1].sect2.CA_position[1])
        distz=abs(self.wingquads[0].sect1.CA_position[2]-self.wingquads[-1].sect2.CA_position[2])
        if distz>disty:
            self.axis=2
            b=distz
        else:
            self.axis=1
            b=disty

        S=0.0
        mac=0.0
        for quad in self.wingquads:
            Sq, macq=quad.calc_reference(axis=self.axis)
            S+=Sq
            mac+=macq
        mac/=S
        return S, mac, b
    def close_tip(self, sectside=2):
        if sectside==2:
            self.wingquads[-1].close_tip(sectside=2)
            self.rightclosed=True
        else:
            self.wingquads[0].close_tip(sectside=1)
            self.leftclosed=True
    def plot_input(self, fig=None, ax=None, show=False, xlim=[], \
        ylim=[], zlim=[], colour='blue'): #plot input geometry data
        if fig==None:
            fig=plt.figure()
        if ax==None:
            plt.axes(projection='3d')
        
        for quad in self.wingquads:
            quad.plot_input(fig=fig, ax=ax, show=False)
        
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        if show:
            plt.show()
    def genwakepanels(self, offset=1000.0, a=0.0, b=0.0):
        for quad in self.wingquads:
            quad.sld.genwakepanels(wakecombs=quad.wakecombs, wakeinds=[[0, 0]], a=a, b=b)