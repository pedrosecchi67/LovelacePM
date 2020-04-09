import numpy as np
from math import *

#file to replace toolkit.f90 functions with raw python functions

def aicm_lines_gen(lines, colpoints):
    nlin=np.size(lines, axis=0)
    npan=np.zise(colpoints, axis=0)

    A=np.zeros((3, npan, nlin))
    B=np.zeros((3, npan, nlin))
    for i in range(nlin):
        A[:, :, i]=lines[i, :, 0]-colpoints
        B[:, :, i]=lines[i, :, 1]-colpoints
    
    na=np.sqrt(A[:, :, 0]**2+A[:, :, 1]**2+A[:, :, 2]**2)
    nb=np.sqrt(B[:, :, 0]**2+B[:, :, 1]**2+B[:, :, 2]**2)

    p=A[:, :, 0]*B[:, :, 0]+A[:, :, 1]*B[:, :, 1]+A[:, :, 2]*B[:, :, 2]

    M=np.zeros((3, npan, nlin))
    M=(A[:, :, 1]*B[:, :, 2]-A[:, :, 2]*B[:, :, 1])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))
    M=(A[:, :, 2]*B[:, :, 0]-A[:, :, 0]*B[:, :, 2])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))
    M=(A[:, :, 0]*B[:, :, 1]-A[:, :, 1]*B[:, :, 0])*(1.0/na+1.0/nb)/((na*nb+p)*(np.pi*4))

    return M