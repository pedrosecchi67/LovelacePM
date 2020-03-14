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
from multiprocessing import *
import warnings as wn
import os
import os.path
from time import *
from func_timeout import *
import threading
import subprocess as sub
import platform as plat

def polar_data(name='n4412', ext_append=True, aseq=[-5.0, 20.0, 1.0], visc=True, Re=3e6, M=0.03, iter=300, flap=None):
    #flap variable: [x_hinge, y_hinge, deflection(angles)]. aseq: same input as required for xfoil command
    if os.path.exists('temppolar.plr'):
        os.remove('temppolar.plr')
    cmds=[]
    if ext_append:
        aflname=name+'.dat'
    else:
        aflname=name
    file_exists=os.path.exists(aflname)
    alphas=np.linspace(aseq[0], aseq[1], aseq[2])
    Cls=np.zeros(np.shape(alphas))
    Cds=np.zeros(np.shape(alphas))
    Cms=np.zeros(np.shape(alphas))
    if file_exists:
        cmds+=['load '+aflname]
        if flap!=None:
            cmds+=['gdes\nflap\n'+str(flap[0])+'\n'+str(flap[1])+'\n'+str(flap[2])+'\n']
        cmds+=['oper']
        if visc:
            cmds+=['visc '+str(Re)]
        cmds+=['Mach '+str(M)]
        cmds+=['iter '+str(iter)]
        cmds+=['pacc']
        cmds+=['temppolar.plr\n']
        cmds+=['aseq']
        cmds+=[str(aseq[0])]
        cmds+=[str(aseq[1])]
        cmds+=[str(aseq[2])]
        cmds+=['\nquit\n']
        xfoil=threading.Thread()
        if plat.system()=='Windows':
            xfoil.p=sub.Popen(['xfoil.exe'], stdin=sub.PIPE, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
        elif plat.system()=='Linux':
            xfoil.p=sub.Popen(['xfoil'], stdin=sub.PIPE, stdout=sub.DEVNULL, stderr=sub.DEVNULL)
        xfoil.start()
        xfoil.join(30)
        xfoil.p.communicate(str.encode('\n'.join(cmds)))
        if xfoil.is_alive():
            return -1
        pfile=open('temppolar.plr', 'r')
        for i in range(12):
            _=pfile.readline()
        totaltext=pfile.read()
        #totaltext=totaltext.split('------ -------- --------- --------- -------- -------- --------\n')[1]
        linelist=totaltext.split('\n')
        alphas=np.array([0.0]*(len(linelist)-1))
        Cls=np.array([0.0]*(len(linelist)-1))
        Cds=np.array([0.0]*(len(linelist)-1))
        Cms=np.array([0.0]*(len(linelist)-1))
        for i in range(len(linelist)-1):
            contentlist=linelist[i].split()
            alphas[i]=float(contentlist[0])
            Cls[i]=float(contentlist[1])
            Cds[i]=float(contentlist[2])
            Cms[i]=float(contentlist[4])
        pfile.close()
        if os.path.exists('temppolar.plr'):
            os.remove('temppolar.plr')
    return alphas, Cls, Cds, Cms

class polar_correction:
    def __init__(self, name='n4412', ext_append=True, aseq=[-10.0, 20.0, 2.0], Re_low=2e6, Re_high=3e6, Mach=0.03, flap=None, iter=300):
        self.Re_low=Re_low
        self.Re_high=Re_high
        self.alphas_Re_low, self.Cls_Re_low, self.Cds_Re_low, self.Cms_Re_low=polar_data(name=name, ext_append=ext_append, \
            aseq=aseq, Re=Re_low, M=Mach, visc=True, iter=iter)
        self.alphas_Re_high, self.Cls_Re_high, self.Cds_Re_high, self.Cms_Re_high=polar_data(name=name, ext_append=ext_append, \
            aseq=aseq, Re=Re_high, M=Mach, visc=True, iter=iter)
        self.alphas_inviscid, self.Cls_inviscid, self.Cds_inviscid, self.Cms_inviscid=polar_data(name=name, ext_append=ext_append, \
            aseq=aseq, M=Mach, visc=False, iter=iter)
    def dump(self, poldir='polars', polname='n4412', ext_append=True, echo=True): #dump to file
        fname=poldir+'/'+polname
        if ext_append:
            fname+='.plr'
        file=open(fname, 'w')
        file.write(str(len(self.alphas_inviscid))+'\n')
        if echo:
            print('Dumping '+polname+' polar')
            print('inviscid')
        for i in range(len(self.alphas_inviscid)):
            file.write(str(self.alphas_inviscid[i])+' '+str(self.Cls_inviscid[i])+' '+\
                str(self.Cds_inviscid[i])+' '+str(self.Cms_inviscid[i])+'\n')
            if echo:
                print(str(self.alphas_inviscid[i])+' '+str(self.Cls_inviscid[i])+' '+\
                    str(self.Cds_inviscid[i])+' '+str(self.Cms_inviscid[i])+'\n')
        file.write(str(self.Re_low)+' '+str(len(self.alphas_Re_low))+'\n')
        if echo:
            print('Re '+str(self.Re_low))
        for i in range(len(self.alphas_Re_low)):
            file.write(str(self.alphas_Re_low[i])+' '+str(self.Cls_Re_low[i])+' '+\
                str(self.Cds_Re_low[i])+' '+str(self.Cms_Re_low[i])+'\n')
            if echo:
                print(str(self.alphas_Re_low[i])+' '+str(self.Cls_Re_low[i])+' '+\
                    str(self.Cds_Re_low[i])+' '+str(self.Cms_Re_low[i])+'\n')
        if echo:
            print('Re '+str(self.Re_high))
        file.write(str(self.Re_high)+' '+str(len(self.alphas_Re_high))+'\n')
        for i in range(len(self.alphas_Re_high)):
            file.write(str(self.alphas_Re_high[i])+' '+str(self.Cls_Re_high[i])+' '+\
                str(self.Cds_Re_high[i])+' '+str(self.Cms_Re_high[i])+'\n')
            if echo:
                print(str(self.alphas_Re_high[i])+' '+str(self.Cls_Re_high[i])+' '+\
                    str(self.Cds_Re_high[i])+' '+str(self.Cms_Re_high[i])+'\n')
        file.close()
    def call(self, Re=2e6, cubic=True):
        eta=np.interp(Re, np.array([self.Re_low, self.Re_high]), np.array([0.0, 1.0]))
        if cubic:
            alphas=sinterp.CubicSpline(self.Cls_inviscid, self.alphas_inviscid)
            Cls_inviscid=np.interp(self.alphas_Re_low, self.alphas_inviscid, self.Cls_inviscid)
            Cds_inviscid=np.interp(self.alphas_Re_low, self.alphas_inviscid, self.Cds_inviscid)
            Cms_inviscid=np.interp(self.alphas_Re_low, self.alphas_inviscid, self.Cms_inviscid)
            Cls=sinterp.CubicSpline(Cls_inviscid, self.Cls_Re_low*(1.0-eta)+self.Cls_Re_high*eta)
            Cds=sinterp.CubicSpline(Cls_inviscid, (self.Cds_Re_low-Cds_inviscid)*(1.0-eta)+(self.Cds_Re_high-Cds_inviscid)*eta)
            Cms=sinterp.CubicSpline(Cls_inviscid, (self.Cms_Re_low-Cms_inviscid)*(1.0-eta)+(self.Cms_Re_high-Cms_inviscid)*eta)
        else:
            alphas=sinterp.interp1d(self.Cls_inviscid, self.alphas_inviscid)
            Cds_inviscid=np.interp(self.alphas_Re_high, self.alphas_inviscid, self.Cds_inviscid)
            Cls_inviscid=np.interp(self.alphas_Re_high, self.alphas_inviscid, self.Cls_inviscid)
            Cms_inviscid=np.interp(self.alphas_Re_high, self.alphas_inviscid, self.Cms_inviscid)
            Cls=sinterp.interp1d(Cls_inviscid, self.Cls_Re_low*(1.0-eta)+self.Cls_Re_high*eta)
            Cds=sinterp.interp1d(Cls_inviscid, (self.Cds_Re_low-Cds_inviscid)*(1.0-eta)+(self.Cds_Re_high-Cds_inviscid)*eta)
            Cms=sinterp.interp1d(Cls_inviscid, (self.Cms_Re_low-Cms_inviscid)*(1.0-eta)+(self.Cms_Re_high-Cms_inviscid)*eta)
        return alphas, Cls, Cds, Cms

def read_polar(poldir='polars', polname='n4412', ext_append=True, echo=True): #read polars from dumped file
    newpolar=polar_correction(name='')
    fname=poldir+'/'+polname
    if ext_append:
        fname+='.plr'
    if not os.path.exists(fname):
        print('WARNING: inexisting polar file provided to read_polar()')
    file=open(fname, 'r')
    alltext=file.read()
    all_lines=alltext.split('\n')
    if echo:
        print('Reading '+polname+' polar')
        print('inviscid')
    ninvisc=int(all_lines[0])
    already_read=1
    newpolar.alphas_inviscid=np.zeros(ninvisc)
    newpolar.Cls_inviscid=np.zeros(ninvisc)
    newpolar.Cds_inviscid=np.zeros(ninvisc)
    newpolar.Cms_inviscid=np.zeros(ninvisc)
    for i in range(already_read, ninvisc+already_read):
        linelist=all_lines[i].split()
        newpolar.alphas_inviscid[i-already_read]=float(linelist[0])
        newpolar.Cls_inviscid[i-already_read]=float(linelist[1])
        newpolar.Cds_inviscid[i-already_read]=float(linelist[2])
        newpolar.Cms_inviscid[i-already_read]=float(linelist[3])
        if echo:
            print(all_lines[i])
    already_read+=ninvisc
    linelist=all_lines[already_read].split()
    nrelow=int(linelist[1])
    newpolar.Re_low=float(linelist[0])
    print('Re '+str(newpolar.Re_low))
    newpolar.alphas_Re_low=np.zeros(nrelow)
    newpolar.Cls_Re_low=np.zeros(nrelow)
    newpolar.Cds_Re_low=np.zeros(nrelow)
    newpolar.Cms_Re_low=np.zeros(nrelow)
    already_read+=1
    for i in range(already_read, nrelow+already_read):
        linelist=all_lines[i].split()
        newpolar.alphas_Re_low[i-already_read]=float(linelist[0])
        newpolar.Cls_Re_low[i-already_read]=float(linelist[1])
        newpolar.Cds_Re_low[i-already_read]=float(linelist[2])
        newpolar.Cms_Re_low[i-already_read]=float(linelist[3])
        if echo:
            print(all_lines[i])
    already_read+=nrelow
    linelist=all_lines[already_read].split()
    nrehigh=int(linelist[1])
    newpolar.Re_high=float(linelist[0])
    print('Re '+str(newpolar.Re_high))
    newpolar.alphas_Re_high=np.zeros(nrehigh)
    newpolar.Cls_Re_high=np.zeros(nrehigh)
    newpolar.Cds_Re_high=np.zeros(nrehigh)
    newpolar.Cms_Re_high=np.zeros(nrehigh)
    already_read+=1
    for i in range(already_read, nrehigh+already_read):
        linelist=all_lines[i].split()
        newpolar.alphas_Re_high[i-already_read]=float(linelist[0])
        newpolar.Cls_Re_high[i-already_read]=float(linelist[1])
        newpolar.Cds_Re_high[i-already_read]=float(linelist[2])
        newpolar.Cms_Re_high[i-already_read]=float(linelist[3])
        if echo:
            print(all_lines[i])
    already_read+=nrehigh
    file.close()
    return newpolar

plr=read_polar()
print(vars(plr))