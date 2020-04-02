import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.interpolate import CubicSpline
import time as tm
import os

import toolkit

def read_afl(afl, afldir='', ext_append=False, header_lines=1, disc=0, strategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, \
    remove_TE_gap=False, extra_intra=False, incidence=0.0, inverse=False, closed=False, sweep=0.0):
    ordir=os.getcwd()
    if len(afldir)==0:
	    os.chdir(ordir)
        #os.chdir(os.path.dirname(os.path.abspath(__file__)))
    else:
	    os.chdir(afldir)
    #read arfoil data points from Selig format file
    if ext_append:
        infile=open(afl+'.dat', 'r')
    else:
        infile=open(afl, 'r')
    aflpts=[]
    alltext=infile.read()
    lines=alltext.split('\n')
    for i in range(header_lines, len(lines)):
        linelist=lines[i].split()
        if len(linelist)!=0:
            aflpts+=[[float(linelist[0]), float(linelist[1])]]
    aflpts=np.array(aflpts)
    if inverse:
        aflpts[:, 1]=-aflpts[:, 1]
    aflpts[:, 1]*=cos(sweep)
    if incidence!=0.0:
        R=np.array([[cos(incidence), sin(incidence)], [-sin(incidence), cos(incidence)]])
        if inverse:
            R=R.T
        aflpts=(R@(aflpts.T)).T/cos(incidence)
    if disc!=0:
        xpts=strategy(np.linspace(1.0, 0.0, disc))
        leading_edge_ind=np.argmin(aflpts[:, 0])
        extra=aflpts[0:leading_edge_ind, :]
        intra=aflpts[leading_edge_ind:np.size(aflpts, 0), :]
        extracs=CubicSpline(np.flip(extra[:, 0]), np.flip(extra[:, 1]))
        extra=np.vstack((xpts, extracs(xpts))).T
        xpts=np.flip(xpts)
        intracs=CubicSpline(intra[:, 0], intra[:, 1])
        intra=np.vstack((xpts, intracs(xpts))).T
        aflpts=np.vstack((extra, intra[1:np.size(intra, 0), :]))
    aflpts[:, 0]-=0.25
    intra[:, 0]-=0.25
    extra[:, 0]-=0.25
    if remove_TE_gap:
        midpoint=(aflpts[-1, :]+aflpts[0, :])/2
        aflpts[-1, :]=midpoint
        aflpts[0, :]=midpoint
    if extra_intra:
        tipind=np.argmin(aflpts[:, 0])
        extra=aflpts[0:tipind, :]
        intra=aflpts[tipind:np.size(aflpts, 0), :]
        extra=np.flip(extra, axis=0)
        return aflpts, extra, intra
    if closed:
        tipind=np.argmin(aflpts[:, 0])
        extra=aflpts[0:tipind+1, :]
        intra=aflpts[tipind:np.size(aflpts, 0), :]
        extra=np.flip(extra, axis=0)
        camberline=(intra+extra)/2
        aflpts[:, 1]=np.interp(aflpts[:, 0], camberline[:, 0], camberline[:, 1])
    os.chdir(ordir)
    return aflpts

def wing_afl_positprocess(afl, gamma=0.0, c=1.0, ypos=0.0, xpos=0.0, zpos=0.0):
    #position airfoil coordinate in 3D axis system
    R=np.array([[1.0, 0.0, 0.0], [0.0, cos(gamma), -sin(gamma)], [0.0, sin(gamma), cos(gamma)]])
    aflnew=(R@np.vstack((afl[:, 0]*c, np.zeros(np.size(afl, 0)), afl[:, 1]*c))).T
    aflnew[:, 1]+=ypos
    aflnew[:, 2]+=zpos
    aflnew[:, 0]+=xpos
    return aflnew

def trimlist(n, l): #trim list to defined length based on first element. For input handling
    if (len(l)<n) and len(l)!=0:
        l+=[l[0]]*(n-len(l))
    return l

def trim_polars(th):
    if th>pi:
        return th-2*pi
    elif th<-pi:
        return th+2*pi

def trim_polars_array(thspacing): #trim polars to eliminate congruous equivalences
    validpos=thspacing>=pi
    thspacing[validpos]-=2*pi
    validpos=thspacing<=-pi
    thspacing[validpos]+=2*pi
    return thspacing

def gen_circdefsect_coords(disc): #present input coordinates for circular defsect
    thetas=np.linspace(-pi, pi, disc)
    return np.vstack((np.sin(thetas), np.cos(thetas))).T

def linear_pts(p1, p2, n, endpoint=False):
    #function to provide linearly interpolated segments in 2D space, serves as tool for other functions in folder
    eta=np.linspace(0.0, 1.0, n, endpoint=endpoint)
    coords=np.zeros((n, 2))
    for i in range(n):
        coords[i, :]=(1.0-eta[i])*p1+eta[i]*p2
    return coords

def elliptic_pts(p1, p2, center, r_x, r_y, th1, th2, n, endpoint=False):
    #function to provide elliptically interpolated segments in 2D space, serves as tool for other functions in folder
    thspacing=np.linspace(th1, th2, n, endpoint=endpoint)
    coords=np.zeros((n, 2))
    coords[:, 0]=np.sin(thspacing)*r_x+center[0]
    coords[:, 1]=np.cos(thspacing)*r_y+center[1]
    return coords

def gen_squaredefsect_coords(disc): #present input coordinates for square defsect
    nside=int(disc/8)
    pts=np.zeros((nside*8, 2))
    #side 1
    pts[0:nside, 0]=np.linspace(0.0, -1.0, nside, endpoint=False)
    pts[0:nside, 1]=-1.0
    #side 2
    pts[nside:3*nside, 1]=np.linspace(-1.0, 1.0, 2*nside, endpoint=False)
    pts[nside:3*nside, 0]=-1.0
    #side 3
    pts[3*nside:5*nside, 0]=np.linspace(-1.0, 1.0, 2*nside, endpoint=False)
    pts[3*nside:5*nside, 1]=1.0
    #side 4
    pts[5*nside:7*nside, 1]=np.linspace(1.0, -1.0, 2*nside, endpoint=False)
    pts[5*nside:7*nside, 0]=1.0
    #side 5
    pts[7*nside:8*nside, 0]=np.linspace(1.0, 0.0, nside, endpoint=True)
    pts[7*nside:8*nside, 1]=-1.0
    
    return pts

def smooth_angle_defsect_coords(r_1x, r_2x, r_1y, r_2y, ldisc=30, thdisc=20):
    #coordinates for body.py's smooth coordinate defsect
    n_low=ldisc
    n_sides=ldisc
    n_up=ldisc
    coords=linear_pts(np.array([0.0, -1.0]), np.array([r_1x-1.0, -1.0]), n_low)
    coords=np.vstack((coords, elliptic_pts(np.array([r_1x-1.0, -1.0]), np.array([-1.0, r_1y-1.0]), np.array([r_1x-1.0, r_1y-1.0]), r_1x, r_1y, -pi, -pi/2, \
        thdisc)))
    coords=np.vstack((coords, linear_pts(np.array([-1.0, r_1y-1.0]), np.array([-1.0, 1.0-r_2y]), n_sides)))
    coords=np.vstack((coords, elliptic_pts(np.array([-1.0, 1.0-r_2y]), np.array([r_2x-1.0, 1.0]), np.array([r_1x-1.0, 1.0-r_2y]), r_2x, r_2y, -pi/2, 0.0, \
        thdisc)))
    coords=np.vstack((coords, linear_pts(np.array([r_2x-1.0, 1.0]), np.array([1.0-r_2x, 1.0]), n_up)))
    coords=np.vstack((coords, elliptic_pts(np.array([1.0-r_2x, 1.0]), np.array([1.0, 1.0-r_2y]), np.array([1.0-r_2x, 1.0-r_2y]), r_2x, r_2y, 0.0, pi/2, \
        thdisc)))
    coords=np.vstack((coords, linear_pts(np.array([1.0, 1.0-r_2y]), np.array([1.0, r_1y-1.0]), n_sides)))
    coords=np.vstack((coords, elliptic_pts(np.array([1.0, r_1y-1.0]), np.array([1.0-r_1x, -1.0]), np.array([1.0-r_1x, r_1y-1.0]), r_1x, r_1y, pi/2, pi, \
        thdisc)))
    coords=np.vstack((coords, linear_pts(np.array([1.0-r_1x, -1.0]), np.array([0.0, -1.0]), n_low, endpoint=True)))
    return coords

def Mtostream(a, b): #define coordinate transformation matrix to convert to streamwise coordinate system
    Mtost=np.zeros((3, 3))
    Mtost[0, :]=np.array([cos(a)*cos(b), -cos(a)*sin(b), sin(a)], dtype='double')
    Mtost[1, :]=np.cross(np.array([0.0, 0.0, 1.0]), Mtost[0, :])
    Mtost[1, :]/=lg.norm(Mtost[1, :])
    Mtost[2, :]=np.cross(Mtost[0, :], Mtost[1, :])
    return Mtost

def Mstreamtouni(a, b): #same but inverted
    return Mtostream(a, b).T

def PG_xmult(beta, a, b): #matrix to which multiply a point array to apply Prandtl-Glauert's correction
    return Mstreamtouni(a, b)@np.diag(np.array([1.0, beta, beta]))@Mtostream(a, b)

def PG_inv_xmult(beta, a, b): #inverse of the matrix before
    return lg.inv(PG_xmult(beta, a, b))

def PG_vtouni(beta, a, b): #returns matrix converting incompressible PG calculated velocities to compressible case
    return Mstreamtouni(a, b)@np.diag(np.array([1.0, beta, beta]))@Mtostream(a, b)