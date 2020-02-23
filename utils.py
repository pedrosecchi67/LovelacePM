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
        aflpts=np.vstack((extra, intra[1:np.size(intra, 0), :]))
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