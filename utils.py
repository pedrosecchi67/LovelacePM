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
    remove_TE_gap=False, extra_intra=False): #read arfoil data points from Selig format file
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
    midpt=(aflpts[-1, :]+aflpts[0, :])/2
    aflpts[-1, :]=midpt
    aflpts[0, :]=midpt
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
        aflpts=np.vstack((extra, intra))
    if remove_TE_gap:
        midpoint=(aflpts[-1, :]+aflpts[0, :])/2
        aflpts[-1, :]=midpoint
        aflpts[0, :]=midpoint
    if not extra_intra:
        return aflpts
    else:
        tipind=np.argmin(aflpts[:, 0])
        extra=aflpts[0:tipind, :]
        intra=aflpts[tipind:np.size(aflpts, 0), :]
        extra=np.flip(extra, axis=0)
        return aflpts, extra, intra

def surf_grad(matdims, colmat, nvectmat, w):
    lambdax=np.linspace(0.0, 1.0, matdims[0])
    lambday=np.linspace(0.0, 1.0, matdims[1])
    dxdlx=np.gradient(colmat[:, :, 0], lambdax, axis=1)
    dydlx=np.gradient(colmat[:, :, 1], lambdax, axis=1)
    dzdlx=np.gradient(colmat[:, :, 2], lambdax, axis=1)
    dxdly=np.gradient(colmat[:, :, 0], lambday, axis=0)
    dydly=np.gradient(colmat[:, :, 1], lambday, axis=0)
    dzdly=np.gradient(colmat[:, :, 2], lambday, axis=0)
    dwdlx=np.gradient(w, lambdax, axis=1)
    dwdly=np.gradient(w, lambday, axis=0)
    l=[]
    for i in range(np.size(dxdlx, 0)):
        for j in range(np.size(dxdlx, 1)):
            l+=[lg.solve(np.vstack((np.array([[dxdlx[i, j], dydlx[i, j], dzdlx[i, j]], \
                [dxdly[i, j], dydly[i, j], dzdly[i, j]]]), nvectmat[i, j, :])), \
                    np.array([dwdlx[i, j], dwdly[i, j], 0.0]))]
    return l