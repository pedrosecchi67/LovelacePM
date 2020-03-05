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
        header_lines=1, xdisc=10, remove_TE_gap=True, inverse=False):
        self.points=wing_afl_positprocess(read_afl(afl=afl, ext_append=True, header_lines=header_lines, disc=xdisc, \
            remove_TE_gap=remove_TE_gap, incidence=incidence, inverse=inverse), gamma=gamma, c=c, xpos=CA_position[0], \
                ypos=CA_position[1], zpos=CA_position[2])

class wing_quadrant: #wing region between two airfoil sections
    def __init__(self, sld, sect1=None, sect2=None):
        #by convention, set section 1 as the rightmost section in the quadrant (positive in y axis)
        #vertical fins should have sect1 as their upmost section
        self.sld=sld
        self.sect1=sect1
        self.sect2=sect2
        self.wakecomb=[]
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
    def patchcompose(self, prevlines={}, strategy=lambda x: x, ldisc=20):
        #code for prevlines: 'intra_left', 'intra_right', 'extra_...
        #and wing quadrant patch additions
        lspacing=strategy(np.linspace(0.0, 1.0, ldisc))
        extrasld=[]
        xdisc=int(np.size(self.sect1.points, 0)/2)+1
        for i in range(xdisc):
            extrasld+=[[]]
            for eta in lspacing:
                extrasld[-1]+=[eta*self.sect1.points[i, :]+(1.0-eta)*self.sect2.points[i, :]]
        intrasld=[]
        for i in range(np.size(self.sect1.points, 0)-1, xdisc-2, -1):
            intrasld+=[[]]
            for eta in lspacing:
                intrasld[-1]+=[eta*self.sect2.points[i, :]+(1.0-eta)*self.sect1.points[i, :]]
        
        #wake info
        self.wakecombs=[]
        
        prev={}
        if 'extra_right' in prevlines:
            prev['right']=prev['extra_right']
        if 'extra_left' in prevlines:
            prev['left']=prev['extra_left']
        horzlines, vertlines, paninds, points=self.sld.addpatch(extrasld, prevlines=prev)
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
        horzlines, vertlines, paninds, points=self.sld.addpatch(intrasld, prevlines=prev, invlats=['up', 'low'])
        TE_intra=paninds[0]
        TE_intra.reverse() #for wake generation
        self.intraright_lines=[vertlines[i][-1] for i in range(len(vertlines))]
        self.intraright_points=[points[i][-1] for i in range(len(points))]
        self.intraleft_lines=[vertlines[i][0] for i in range(len(vertlines))]
        self.intraleft_points=[points[i][0] for i in range(len(points))]

        #wake info
        for i in range(len(TE_extra)):
            self.wakecombs+=[[TE_extra[i], TE_intra[i]]]