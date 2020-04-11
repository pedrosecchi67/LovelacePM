import numpy as np
from math import *
import numpy.linalg as lg

#file to replace toolkit.f90 functions with raw python functions

def aicm_lines_gen(lines, colpoints):
    nlin=np.size(lines, axis=0)
    npan=np.size(colpoints, axis=0)

    A=np.zeros((nlin, 3))
    B=np.zeros((nlin, 3))
    p=np.zeros(nlin)
    na=np.zeros(nlin)
    nb=np.zeros(nlin)
    M=np.zeros((3, npan, nlin))
    for i in range(npan):
        A=(lines[:, :, 0]-colpoints[i, :])
        B=(lines[:, :, 1]-colpoints[i, :])
        
        na=np.sqrt(A[:, 0]**2+A[:, 1]**2+A[:, 2]**2)
        nb=np.sqrt(B[:, 0]**2+B[:, 1]**2+B[:, 2]**2)

        p=A[:, 0]*B[:, 0]+A[:, 1]*B[:, 1]+A[:, 2]*B[:, 2]

        M[0, i, :]=(A[:, 1]*B[:, 2]-A[:, 2]*B[:, 1])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))
        M[1, i, :]=(A[:, 2]*B[:, 0]-A[:, 0]*B[:, 2])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))
        M[2, i, :]=(A[:, 0]*B[:, 1]-A[:, 1]*B[:, 0])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))

    return M

def aicm_norm_conv(aicm3, nvectmat):
    return (aicm3[0, :, :].T*nvectmat[:, 0]).T+(aicm3[1, :, :].T*nvectmat[:, 1]).T+(aicm3[2, :, :].T*nvectmat[:, 2]).T

def isitin(pts, pcont, centers, nvect, isvalid, tolerance, npan): #the function capable of ruining any man's pride!
    #isin=np.ones(npan, dtype='bool')
    isin=isvalid
    dot=np.zeros(npan)
    side=np.zeros((npan, 3))
    vec=np.zeros((npan, 3))
    prod=np.zeros((npan, 3))
    #vec=pcont-centers
    #dot=nvect[:, 0]*vec[:, 0]+nvect[:, 1]*vec[:, 1]+nvect[:, 2]*vec[:, 2]
    #isin=np.abs(dot)<tolerance
    for i in range(4):
        side[isin, :]=pts[isin, :, (i+1)%4]-pts[isin, :, i]
        vec[isin, :]=pcont[isin, :]-pts[isin, :, i]
        prod[isin, 0]=vec[isin, 1]*side[isin, 2]-vec[isin, 2]*side[isin, 1]
        prod[isin, 1]=vec[isin, 2]*side[isin, 0]-vec[isin, 0]*side[isin, 2]
        prod[isin, 2]=vec[isin, 0]*side[isin, 1]-vec[isin, 1]*side[isin, 0]
        dot[isin]=prod[isin, 0]*nvect[isin, 0]+prod[isin, 1]*nvect[isin, 1]+prod[isin, 2]*nvect[isin, 2]
        isin[isin]=np.logical_and(isin[isin], np.logical_or(dot[isin]<tolerance, lg.norm(side[isin, :], axis=1)<tolerance))
    return isin

def get_panel_contact(p, u, centers, nvects, pts, tolerance):
    npan=np.size(nvects, axis=0)
    #remember: shape(pts)==(npan, 3, 4)
    pnorm=(p[0]-centers[:, 0])*nvects[:, 0]+(p[1]-centers[:, 1])*nvects[:, 1]+(p[2]-centers[:, 2])*nvects[:, 2]
    unorm=(u@nvects.T).T
    isvalid=unorm<-tolerance
    lambdas=np.zeros(npan)
    lambdas[isvalid]=-pnorm[isvalid]/unorm[isvalid]
    pcont=np.zeros((npan, 3))
    pcont[:, 0]=p[0]+u[0]*lambdas
    pcont[:, 1]=p[1]+u[1]*lambdas
    pcont[:, 2]=p[2]+u[2]*lambdas
    isin=isitin(pts, pcont, centers, nvects, isvalid, tolerance, npan)
    error=not np.any(isin)
    if not error:
        pcont=pcont[isin, :]
        pcont=pcont[0, :]
    else:
        pcont=np.zeros(3)
    return pcont, error

def get_field_influence(lines, solution_lines, colpoint, tolerance):
    dv=np.zeros(3, dtype='double')
    nlin=np.size(lines, axis=0)

    A=np.zeros((nlin, 3))
    B=np.zeros((nlin, 3))
    p=np.zeros(nlin)
    na=np.zeros(nlin)
    nb=np.zeros(nlin)
    A=(lines[:, :, 0]-colpoint)
    B=(lines[:, :, 1]-colpoint)
    
    na=np.sqrt(A[:, 0]**2+A[:, 1]**2+A[:, 2]**2)
    nb=np.sqrt(B[:, 0]**2+B[:, 1]**2+B[:, 2]**2)

    p=A[:, 0]*B[:, 0]+A[:, 1]*B[:, 1]+A[:, 2]*B[:, 2]

    dv[0]=((A[:, 1]*B[:, 2]-A[:, 2]*B[:, 1])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4)))@solution_lines
    dv[1]=((A[:, 2]*B[:, 0]-A[:, 0]*B[:, 2])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4)))@solution_lines
    dv[2]=((A[:, 0]*B[:, 1]-A[:, 1]*B[:, 0])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4)))@solution_lines

    return dv

def aicm_lines_recalc(lineinds, lines, colpoints):
    nlin=len(lineinds)
    npan=np.size(colpoints, axis=0)

    A=np.zeros((nlin, 3))
    B=np.zeros((nlin, 3))
    p=np.zeros(nlin)
    na=np.zeros(nlin)
    nb=np.zeros(nlin)
    M=np.zeros((3, npan, nlin))
    for i in range(npan):
        A=(lines[lineinds, :, 0]-colpoints[i, :])
        B=(lines[lineinds, :, 1]-colpoints[i, :])
        
        na=np.sqrt(A[:, 0]**2+A[:, 1]**2+A[:, 2]**2)
        nb=np.sqrt(B[:, 0]**2+B[:, 1]**2+B[:, 2]**2)

        p=A[:, 0]*B[:, 0]+A[:, 1]*B[:, 1]+A[:, 2]*B[:, 2]

        M[0, i, :]=(A[:, 1]*B[:, 2]-A[:, 2]*B[:, 1])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))
        M[1, i, :]=(A[:, 2]*B[:, 0]-A[:, 0]*B[:, 2])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))
        M[2, i, :]=(A[:, 0]*B[:, 1]-A[:, 1]*B[:, 0])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))

    return M