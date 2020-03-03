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
            self.polar_rule=sinterp.CubicSpline(self.thetas, self.Rs)
        else:
            self.polar_rule=sinterp.interp1d(self.thetas, self.Rs)
        #define function for interpolation in polar coordinates between theta in respect to z axis as input

class body_panel: #exclusively for intersection finding purposes
    def __init__(self, p1, p2, p3, p4, tolerance=0.00005):
        points=np.vstack((p1, p2, p3, p4)).T
        self.points=points
        self.center, self.Mtosys, self.Mtouni, self.locpoints=toolkit.body_panel_process(self.points, tolerance)
        #FORTRAN was used here for speed up purposes

class body: #body definition class for center definition to obtain polar cooridinates and creation of input panels 
    #for abutments
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
        self.center=lambda x: np.array([x, centerfunc_y(x), centerfunc_z(x)])
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