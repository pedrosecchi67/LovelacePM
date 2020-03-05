import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from scipy.interpolate import CubicSpline
import time as tm

import toolkit

def read_afl(afl, ext_append=False, header_lines=1, disc=0, strategy=lambda x: (np.sin(pi*x-pi/2)+1)/2, \
    remove_TE_gap=False, extra_intra=False, incidence=0.0, inverse=False): #read arfoil data points from Selig format file
    if ext_append:
        infile=open(afl+'.dat', 'r')
    else:
        infile=open(afl, 'r')
    aflpts=[]
    alltext=infile.read()
    lines=alltext.split('\n')
    for i in range(header_lines, len(lines)):
        linelist=lines[i].split()
        aflpts+=[[float(linelist[0]), float(linelist[1])]]
    aflpts=np.array(aflpts)
    if incidence!=0.0:
        R=np.array([[cos(incidence), -sin(incidence)], [sin(incidence), cos(incidence)]])
        aflpts=(R@(aflpts.T)).T/cos(incidence)
    if disc!=0:
        leading_edge_ind=np.argmin(aflpts[:, 0])
        extra=aflpts[0:leading_edge_ind, :]
        intra=aflpts[leading_edge_ind:np.size(aflpts, 0), :]
        xpts=strategy(np.linspace(1.0, 0.0, disc))
        extracs=CubicSpline(np.flip(extra[:, 0]), np.flip(extra[:, 1]))
        extra=np.vstack((xpts, extracs(xpts))).T
        xpts=strategy(np.linspace(0.0, 1.0, disc))
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
    if inverse:
        aflpts[:, 1]=-aflpts[:, 1]
    if not extra_intra:
        return aflpts
    else:
        tipind=np.argmin(aflpts[:, 0])
        extra=aflpts[0:tipind, :]
        intra=aflpts[tipind:np.size(aflpts, 0), :]
        extra=np.flip(extra, axis=0)
        return aflpts, extra, intra

def wing_afl_positprocess(afl, gamma=0.0, c=1.0, ypos=0.0, xpos=0.0, zpos=0.0):
    #position airfoil coordinate in 3D axis system
    R=np.array([[cos(gamma), sin(gamma), 0.0], [-sin(gamma), cos(gamma), 0.0], [0.0, 0.0, 1.0]])
    aflnew=(R@np.vstack((afl[:, 0]*c+xpos, np.zeros(np.size(afl, 0)), afl[:, 1]*c))).T
    aflnew[:, 1]+=ypos
    aflnew[:, 2]+=zpos
    aflnew[:, 0]+=xpos
    return aflnew

def trimlist(n, l): #trim list to defined length based on first element. For input handling
    while(len(l)<n):
        l+=[l[0]]
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