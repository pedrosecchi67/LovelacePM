import numpy as np
from math import *

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
    return aicm3[0, :, :]*nvectmat[:, 0]+aicm3[1, :, :]*nvectmat[:, 1]+aicm3[2, :, :]*nvectmat[:, 2]