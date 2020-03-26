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

''' Xfoil automation to temporarily supply the absense of a parallel, estimative boundary layer solver '''

def polar_data(name='n4412', afldir='', ext_append=True, aseq=[-5.0, 20.0, 1.0], visc=True, Re=3e6, M=0.03, iter=300, flap=None, npan=300, LE_con=0.4, inverse=False):
    #flap variable: [x_hinge, y_hinge, deflection(angles)]. aseq: same input as required for xfoil command
    ordir=os.getcwd()
    if len(afldir)==0:
        os.chdir(ordir)
    else:
        os.chdir(afldir)
    if os.path.exists('temppolar.plr'):
        os.remove('temppolar.plr')
    cmds=[]
    if ext_append:
        aflname=name+'.dat'
    else:
        aflname=name
    file_exists=os.path.exists(aflname)
    alphas=np.arange(aseq[0], aseq[1]+aseq[2], aseq[2])
    Cls=2*pi*np.radians(alphas)
    Cds=np.zeros(np.shape(alphas))
    Cms=np.zeros(np.shape(alphas))
    if file_exists:
        cmds+=['load '+aflname]
        if flap!=None:
            cmds+=['gdes\nflap\n'+str(flap[0])+'\n'+str(flap[1])+'\n'+str(flap[2])+'\n']
        cmds+=['ppar']
        cmds+=['n '+str(npan)]
        cmds+=['r '+str(LE_con)]
        cmds+=[' ']*2
        cmds+=['oper']
        if visc:
            cmds+=['visc '+str(Re)]
        cmds+=['Mach '+str(M)]
        cmds+=['iter '+str(iter)]
        cmds+=['pacc']
        cmds+=['temppolar.plr\n']
        for a in np.arange(aseq[0], aseq[1], aseq[2]):
            cmds+=['a']
            cmds+=[str((-1 if inverse else 1)*a)]
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
            if visc:
                Cds[i]=float(contentlist[2])
            else:
                Cds[i]=float(contentlist[3])
            Cms[i]=float(contentlist[4])
        pfile.close()
        if os.path.exists('temppolar.plr'):
            os.remove('temppolar.plr')
        if os.path.exists('*.bl'):
            os.remove('*.bl')
    os.chdir(ordir)
    return alphas, Cls, Cds, Cms

class polar_correction:
    def __init__(self, name='n4412', afldir='', ext_append=True, aseq=[-5.0, 20.0, 2.0], Re_low=2e6, Re_high=3e6, Mach=0.03, flap=None, iter=300, cubic=True):
        self.Re_low=Re_low
        self.Re_high=Re_high
        amed=(aseq[0]+aseq[1])/2
        alphas_Re_low1, Cls_Re_low1, Cds_Re_low1, Cms_Re_low1=polar_data(name=name, afldir=afldir, ext_append=ext_append, \
            aseq=[amed, aseq[1], aseq[2]], Re=Re_low, M=Mach, visc=True, iter=iter)
        alphas_Re_high1, Cls_Re_high1, Cds_Re_high1, Cms_Re_high1=polar_data(name=name, afldir=afldir, ext_append=ext_append, \
            aseq=[amed, aseq[1], aseq[2]], Re=Re_high, M=Mach, visc=True, iter=iter)
        alphas_Re_low2, Cls_Re_low2, Cds_Re_low2, Cms_Re_low2=polar_data(name=name, afldir=afldir, ext_append=ext_append, \
            aseq=[amed-aseq[2], aseq[0], -aseq[2]], Re=Re_low, M=Mach, visc=True, iter=iter)
        alphas_Re_high2, Cls_Re_high2, Cds_Re_high2, Cms_Re_high2=polar_data(name=name, afldir=afldir, ext_append=ext_append, \
            aseq=[amed-aseq[2], aseq[0], -aseq[2]], Re=Re_high, M=Mach, visc=True, iter=iter)
        self.alphas_Re_low=np.hstack((alphas_Re_low1, alphas_Re_low2))
        self.Cls_Re_low=np.hstack((Cls_Re_low1, Cls_Re_low2))
        self.Cds_Re_low=np.hstack((Cds_Re_low1, Cds_Re_low2))
        self.Cms_Re_low=np.hstack((Cms_Re_low1, Cms_Re_low2))
        order=np.argsort(self.alphas_Re_low)
        self.alphas_Re_low=self.alphas_Re_low[order]
        self.Cls_Re_low=self.Cls_Re_low[order]
        self.Cds_Re_low=self.Cds_Re_low[order]
        self.Cms_Re_low=self.Cms_Re_low[order]
        self.alphas_Re_high=np.hstack((alphas_Re_high1, alphas_Re_high2))
        self.Cls_Re_high=np.hstack((Cls_Re_high1, Cls_Re_high2))
        self.Cds_Re_high=np.hstack((Cds_Re_high1, Cds_Re_high2))
        self.Cms_Re_high=np.hstack((Cms_Re_high1, Cms_Re_high2))
        order=np.argsort(self.alphas_Re_high)
        self.alphas_Re_high=self.alphas_Re_high[order]
        self.Cls_Re_high=self.Cls_Re_high[order]
        self.Cds_Re_high=self.Cds_Re_high[order]
        self.Cms_Re_high=self.Cms_Re_high[order]
        self.alphas_inviscid, self.Cls_inviscid, self.Cds_inviscid, self.Cms_inviscid=polar_data(name=name, afldir=afldir, ext_append=ext_append, \
            aseq=aseq, M=Mach, visc=False, iter=iter)
        if len(name)!=0:
            self.create_functions(cubic=cubic)
    def dump(self, poldir='', polname='n4412', ext_append=True, echo=True): #dump to file
        ordir=os.getcwd()
        if len(poldir)==0:
            os.chdir(ordir)
        else:
            os.chdir(poldir)
        #fname=poldir+'/'+polname
        fname=polname
        if ext_append:
            fname+='.plr'
        file=open(fname, 'r')
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
        os.chdir(ordir)
    def create_functions(self, cubic=True):
        Cls_inviscid_Re_low=np.interp(self.alphas_Re_low, self.alphas_inviscid, self.Cls_inviscid)
        Cls_inviscid_Re_high=np.interp(self.alphas_Re_high, self.alphas_inviscid, self.Cls_inviscid)
        Cds_inviscid_Re_low=np.interp(self.alphas_Re_low, self.alphas_inviscid, self.Cds_inviscid)
        Cds_inviscid_Re_high=np.interp(self.alphas_Re_high, self.alphas_inviscid, self.Cds_inviscid)
        Cms_inviscid_Re_low=np.interp(self.alphas_Re_low, self.alphas_inviscid, self.Cms_inviscid)
        Cms_inviscid_Re_high=np.interp(self.alphas_Re_high, self.alphas_inviscid, self.Cms_inviscid)
        if cubic:
            self.alphas_fun=sinterp.CubicSpline(self.Cls_inviscid, self.alphas_inviscid)
            self.Cls_Re_low_fun=sinterp.CubicSpline(Cls_inviscid_Re_low, self.Cls_Re_low)
            self.Cls_Re_high_fun=sinterp.CubicSpline(Cls_inviscid_Re_high, self.Cls_Re_high)
            self.Cds_Re_low_fun=sinterp.CubicSpline(Cls_inviscid_Re_low, self.Cds_Re_low-Cds_inviscid_Re_low)
            self.Cds_Re_high_fun=sinterp.CubicSpline(Cls_inviscid_Re_high, self.Cds_Re_high-Cds_inviscid_Re_high)
            self.Cms_Re_low_fun=sinterp.CubicSpline(Cls_inviscid_Re_low, self.Cms_Re_low-Cms_inviscid_Re_low)
            self.Cms_Re_high_fun=sinterp.CubicSpline(Cls_inviscid_Re_high, self.Cms_Re_high-Cms_inviscid_Re_high)
        else:
            self.alphas_fun=sinterp.interp1d(self.Cls_inviscid, self.alphas_inviscid)
            self.Cls_Re_low_fun=sinterp.interp1d(Cls_inviscid_Re_low, self.Cls_Re_low)
            self.Cls_Re_high_fun=sinterp.interp1d(Cls_inviscid_Re_high, self.Cls_Re_high)
            self.Cds_Re_low_fun=sinterp.interp1d(Cls_inviscid_Re_low, self.Cds_Re_low-Cds_inviscid_Re_low)
            self.Cds_Re_high_fun=sinterp.interp1d(Cls_inviscid_Re_high, self.Cds_Re_high-Cds_inviscid_Re_high)
            self.Cms_Re_low_fun=sinterp.interp1d(Cls_inviscid_Re_low, self.Cms_Re_low-Cms_inviscid_Re_low)
            self.Cms_Re_high_fun=sinterp.interp1d(Cls_inviscid_Re_high, self.Cms_Re_high-Cms_inviscid_Re_high)
    def __call__(self, Re=2e6, inverse=False): #return fitting functions for a given Reynolds number
        eta=np.interp(Re, np.array([self.Re_low, self.Re_high]), np.array([0.0, 1.0]))
        if not inverse:
            alphas=lambda Cl: self.alphas_fun(-Cl)
            Cls=lambda Cl: self.Cls_Re_low_fun(Cl)*(1.0-eta)+eta*self.Cls_Re_high_fun(Cl)
            Cds=lambda Cl: self.Cds_Re_low_fun(Cl)*(1.0-eta)+eta*self.Cds_Re_high_fun(Cl)
            Cms=lambda Cl: self.Cms_Re_low_fun(Cl)*(1.0-eta)+eta*self.Cms_Re_high_fun(Cl)
        else:
            alphas=lambda Cl: -self.alphas_fun(-Cl)
            Cls=lambda Cl: -(self.Cls_Re_low_fun(-Cl)*(1.0-eta)+eta*self.Cls_Re_high_fun(-Cl))
            Cds=lambda Cl: self.Cds_Re_low_fun(-Cl)*(1.0-eta)+eta*self.Cds_Re_high_fun(-Cl)
            Cms=lambda Cl: -(self.Cms_Re_low_fun(-Cl)*(1.0-eta)+eta*self.Cms_Re_high_fun(-Cl))
        return alphas, Cls, Cds, Cms

def read_polar(poldir='', polname='n4412', ext_append=True, echo=True, cubic=True): #read polars from dumped file
    ordir=os.getcwd()
    if len(poldir)==0:
        os.chdir(ordir)
    else:
        os.chdir(poldir)
    newpolar=polar_correction(name='')
    #fname=poldir+'/'+polname
    fname=polname
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
    newpolar.create_functions(cubic=cubic)
    os.chdir(ordir)
    return newpolar