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
from body import *
from wing import *

def plot_Cps(sld, elems=[], xlim=[], ylim=[], zlim=[-2.5, 2.5]):
    fig=plt.figure()
    ax=plt.axes(projection='3d')
    xs=[]
    ys=[]
    cps=[]
    for elem in elems:
        for i in elem.paninds:
            if sld.panels[i].nvector[2]>=0:
                xs+=[sld.panels[i].colpoint[0]]
                ys+=[sld.panels[i].colpoint[1]]
                cps+=[sld.Cps[i]]
    ax.scatter3D(xs, ys, cps, 'blue')
    xs=[]
    ys=[]
    cps=[]
    for elem in elems:
        for i in elem.paninds:
            if sld.panels[i].nvector[2]<0:
                xs+=[sld.panels[i].colpoint[0]]
                ys+=[sld.panels[i].colpoint[1]]
                cps+=[sld.Cps[i]]
    ax.scatter3D(xs, ys, cps, 'red')
    if len(xlim)!=0:
        ax.set_xlim3d(xlim[0], xlim[1])
    if len(ylim)!=0:
        ax.set_ylim3d(ylim[0], ylim[1])
    if len(zlim)!=0:
        ax.set_zlim3d(zlim[0], zlim[1])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('$C_p$')
    plt.show()

def plot_Cls(sld, alpha=0.0, wings=[], axis=1):
    ypos=[0.0]
    plotcorr=False
    n=1
    for w in wings:
        if not w.coefavailable:
            w.calc_coefs(axis=axis, alpha=alpha)
        ypos+=[np.amax(w.Cls), np.amin(w.Cls)]
        plt.plot(w.ys, w.Cls, label='Inviscid w'+str(n))
        if w.plotcorrections:
            plt.plot(w.ys, w.Cls_corrected, label='Corrected w'+str(n))
            ypos+=[np.amax(w.Cls_corrected), np.amin(w.Cls_corrected)]
            plotcorr=True
        n+=1
    plt.title('Sectional lift coefficient')
    plt.xlabel('y')
    plt.ylabel('$C_l$')
    plt.ylim(min(ypos)-0.1*abs(min(ypos)), max(ypos)+0.1*abs(max(ypos)))
    plt.grid()
    if plotcorr:
        plt.legend()
    plt.show()

def plot_gammas(sld, alpha=0.0, Uinf=1.0, wings=[], axis=1):
    ypos=[0.0]
    for w in wings:
        if not w.coefavailable:
            w.calc_coefs(axis=axis, alpha=alpha)
        #gammas=w.Cls*w.cs*Uinf/2
        gammas=w.Gammas
        ypos+=[np.amax(gammas), np.amin(gammas)]
        plt.plot(w.ys, gammas, 'blue', label='Inviscid')
    plt.title('Sectional circulation')
    plt.xlabel('y')
    plt.ylabel('$\Gamma$')
    plt.ylim(min(ypos)-0.1*abs(min(ypos)), max(ypos)+0.1*abs(max(ypos)))
    plt.grid()
    plt.show()

def plot_Cds(sld, alpha=0.0, wings=[], axis=1):
    ypos=[0.0]
    plotcorr=False
    n=1
    for w in wings:
        if not w.coefavailable:
            w.calc_coefs(axis=axis, alpha=alpha)
        ypos+=[np.amax(w.Cds), np.amin(w.Cds)]
        plt.plot(w.ys, w.Cds, label='Inviscid w'+str(n))
        if w.plotcorrections:
            plt.plot(w.ys, w.Cds_corrected, label='Corrected w'+str(n))
            ypos+=[np.amax(w.Cds_corrected), np.amin(w.Cds_corrected)]
            plotcorr=True
        n+=1
    plt.title('Sectional drag coefficient')
    plt.xlabel('y')
    plt.ylabel('$C_d$')
    plt.ylim(min(ypos)-0.1*abs(min(ypos)), max(ypos)+0.1*abs(max(ypos)))
    plt.grid()
    if plotcorr:
        plt.legend()
    plt.show()

def plot_Cms(sld, alpha=0.0, wings=[], axis=1):
    ypos=[0.0]
    plotcorr=False
    n=1
    for w in wings:
        if not w.coefavailable:
            w.calc_coefs(axis=axis, alpha=alpha)
        ypos+=[np.amax(w.Cms), np.amin(w.Cms)]
        plt.plot(w.ys, w.Cms, label='Inviscid w'+str(n))
        if w.plotcorrections:
            plt.plot(w.ys, w.Cms_corrected, label='Corrected w'+str(n))
            ypos+=[np.amax(w.Cms_corrected), np.amin(w.Cms_corrected)]
            plotcorr=True
        n+=1
    plt.title('Sectional quarter-chord torque coefficient')
    plt.xlabel('y')
    plt.ylabel('$C_m1/4$')
    plt.ylim(min(ypos)-0.1*abs(min(ypos)), max(ypos)+0.1*abs(max(ypos)))
    plt.grid()
    if plotcorr:
        plt.legend()
    plt.show()