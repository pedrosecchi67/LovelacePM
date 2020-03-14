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

class aircraft: #class to ease certain case studies for full aircraft
    #generate aircraft BEFORE patchcompose functions
    def parameter_report(self):
        print('Freestream parameters:')
        print('%10s %10s %10s %10s %10s %10s' % ('a', 'b', 'p', 'q', 'r', 'Uinf'))
        print('%10f %10f %10f %10f %10f %10f' % (degrees(self.a), degrees(self.b), self.p, self.q, self.r, self.Uinf))
        print('Geometry parameters:')
        print('%10s %10s %10s' % ('Sref', 'cref', 'bref'))
        print('%10f %10f %10f' % (self.Sref, self.cref, self.bref))
        if self.hascontrol():
            print('Control parameters:')
            controlnames_str=''
            controlvals_str=''
            for cname in self.controlset:
                controlnames_str+='%10s' % cname
                controlvals_str+='%10f' % degrees(self.controlset[cname].state)
            print(controlnames_str)
            print(controlvals_str)
        '''if self.massavailable:
            print('Mass parameters not yet present')
        else:
            print('%10s %10s %10s %10s %10s %10s %10s %10s' % ('m', 'Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Ixz', 'Iyz'))
            print('%10f %10f %10f %10f %10f %10f %10f %10f' % (self.m, self.Ixx, self.Iyy, self.Izz, self.Ixy, self.Ixz, self.Iyz))'''
    def __init__(self, sld, elems=[], Sref=0.0, cref=0.0, bref=0.0, echo=True):
        self.wings=[]
        self.bodies=[]
        for e in elems:
            e.set_aircraft(self)
            if type(e)==wing:
                self.wings+=[e]
            elif type(e)==body:
                self.bodies+=[e]

        #BE VARY AWARE OF THESE WARNINGS WHEN USINGG FUSELAGE-SEPARATED SEMI-WINGS AS WING OBJECTS
        S, mac, b=self.wings[0].calc_reference()
        if Sref==0.0:
            print('WARNING: no reference surface provided. Using first wing input as reference. Address to aircraft constructor\'s Sref argument for better definitions')
            self.Sref=S
        else:
            self.Sref=Sref
        if cref==0.0:
            print('WARNING: no reference chord provided. Using first wing input as reference. Address to aircraft constructor\'s cref argument for better definitions')
            self.cref=mac
        else:
            self.cref=cref
        if bref==0.0:
            print('WARNING: no reference span provided. Using first wing input as reference. Address to aircraft constructor\'s bref argument for better definitions')
            self.bref=b
        else:
            self.bref=bref

        self.a=0.0
        self.b=0.0
        self.p=0.0
        self.q=0.0
        self.r=0.0
        self.Uinf=1.0
        
        #defining available controls
        self.controlset={}
        for wng in self.wings:
            for wngqd in wng.wingquads:
                if wngqd.hascontrol():
                    for cname in wngqd.controls:
                        if not cname in self.controlset:
                            self.controlset[cname]=control_DOF()
                            wngqd.controls[cname].DOF=self.controlset[cname]

        if echo:
            self.parameter_report()

        self.sld=sld
        self.forcesavailable=False
        self.stabavailable=False
        self.massavailable=False

        #defining plotting limits
        xkeypts=[]
        ykeypts=[]
        zkeypts=[]
        for wng in self.wings:
            for wngqd in wng.wingquads:
                xkeypts+=[wngqd.sect1.CA_position[0]-wngqd.sect1.c/4, wngqd.sect1.CA_position[0]+3*wngqd.sect1.c/4, \
                    wngqd.sect2.CA_position[0]-wngqd.sect2.c/4, wngqd.sect2.CA_position[0]+3*wngqd.sect2.c/4]
                ykeypts+=[wngqd.sect1.CA_position[1], wngqd.sect2.CA_position[1]]
                zkeypts+=[wngqd.sect1.CA_position[2], wngqd.sect2.CA_position[2]]
        for bdy in self.bodies:
            xkeypts+=[bdy.sections[i].center[0] for i in range(len(bdy.sections))]
            ykeypts+=[bdy.sections[i].center[1]-bdy.sections[i].Rs[i] for i in range(len(bdy.sections))]
            zkeypts+=[bdy.sections[i].center[2]-bdy.sections[i].Rs[i] for i in range(len(bdy.sections))]
            ykeypts+=[bdy.sections[i].center[1]+bdy.sections[i].Rs[i] for i in range(len(bdy.sections))]
            zkeypts+=[bdy.sections[i].center[2]+bdy.sections[i].Rs[i] for i in range(len(bdy.sections))]
        xmax=max(xkeypts)
        xmin=min(xkeypts)
        ymax=max(ykeypts)
        ymin=min(ykeypts)
        zmax=max(zkeypts)
        zmin=min(zkeypts)
        self.plotlim=max([abs(xmax), abs(ymax), abs(zmax), abs(xmin), abs(ymin), abs(zmin)])
    def hascontrol(self):
        return len(self.controlset)!=0
    def addwake(self, offset=1000.0):
        #once the wake has been added, one can end geometry pre-processing of panels in question:
        self.sld.end_preprocess()
        for wng in self.wings:
            wng.genwakepanels(offset=offset, a=self.a, b=self.b)
    def calcforces(self, echo=True, custCp=np.array([])):
        self.CX=sum([self.sld.forces[i][0] for i in range(self.sld.npanels)])/self.Sref
        self.CY=sum([self.sld.forces[i][1] for i in range(self.sld.npanels)])/self.Sref
        self.CZ=sum([self.sld.forces[i][2] for i in range(self.sld.npanels)])/self.Sref
        u=np.array([cos(self.a)*cos(self.b), sin(self.b)*cos(self.a), sin(self.a)])
        v=np.array([-sin(self.a)*cos(self.b), -sin(self.b)*sin(self.a), cos(self.a)])
        adforces=np.array([self.CX, self.CY, self.CZ])
        self.CL=v@adforces
        self.CD=u@adforces
        self.Cl=sum([self.sld.moments[i][0] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
        self.Cm=sum([self.sld.moments[i][1] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
        self.Cn=sum([self.sld.moments[i][2] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
        if echo:
            self.forces_report()
    def calcstab(self, echo=True):
        #calculate variation on Cps in respect to several parameters
        pars=['Uinf', 'a', 'b', 'p', 'q', 'r']
        self.stabderivative_dict={}
        for p in pars:
            dCps=self.sld.calc_derivative(self.Uinf, a=self.a, b=self.b, p=self.p, q=self.q, r=self.r, par=p)
            forces=[-dCps[i]*self.sld.panels[i].S*self.sld.panels[i].nvector for i in range(self.sld.npanels)]
            moments=[np.cross(self.sld.panels[i].colpoint, forces[i]) for i in range(self.sld.npanels)]
            dCX=sum([forces[i][0] for i in range(self.sld.npanels)])/self.Sref
            dCY=sum([forces[i][1] for i in range(self.sld.npanels)])/self.Sref
            dCZ=sum([forces[i][2] for i in range(self.sld.npanels)])/self.Sref
            u=np.array([cos(self.a)*cos(self.b), sin(self.b)*cos(self.a), sin(self.a)])
            v=np.array([-sin(self.a)*cos(self.b), -sin(self.b)*sin(self.a), cos(self.a)])
            adforces=np.array([dCX, dCY, dCZ])
            dCL=v@adforces
            dCD=u@adforces
            dCl=sum([moments[i][0] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
            dCm=sum([moments[i][1] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
            dCn=sum([moments[i][2] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
            self.stabderivative_dict[p]={'dCX':dCX, 'dCY':dCY, 'dCZ':dCZ, 'dCL':dCL, 'dCD':dCD, 'dCl':dCl, 'dCm':dCm, 'dCn':dCn}
        self.stabavailable=True
        if echo:
            self.stabreport()
    def stabreport(self): #report calculated stability derivatives, or calculate them in place
        if not self.stabavailable:
            self.calcstab(echo=False)
        print('%5s | %10s %10s %10s %10s %10s %10s %10s %10s' % ('par', 'dCX', 'dCY', 'dCZ', 'dCL', 'dCD', 'dCl', 'dCm', 'dCn'))
        for p in self.stabderivative_dict:
            print('%5s | %10f %10f %10f %10f %10f %10f %10f %10f' % (p, self.stabderivative_dict[p]['dCX'], self.stabderivative_dict[p]['dCY'], \
                self.stabderivative_dict[p]['dCZ'], self.stabderivative_dict[p]['dCL'], self.stabderivative_dict[p]['dCD'], self.stabderivative_dict[p]['dCl'], \
                    self.stabderivative_dict[p]['dCm'], self.stabderivative_dict[p]['dCn']))
    def edit_parameters(self, pardict):
        for par in pardict:
            val=pardict[par]
            if par=='a':
                self.a=radians(val)
            elif par=='b':
                self.b=radians(val)
            elif par=='p':
                self.p=val
            elif par=='q':
                self.q=val
            elif par=='r':
                self.r=val
            elif par=='Uinf':
                self.Uinf=val
            else:
                self.controlset[par].state=radians(val)
        #reset result readiness
        self.stabavailable=False
        self.forcesavailable=False
        self.parameter_report()
    def eulersolve(self, echo=True, damper=0.0):
        self.sld.eulersolve(target=np.array([]), a=self.a, b=self.b, p=self.p, q=self.q, r=self.r, damper=damper, Uinf=self.Uinf, echo=True)
    def forces_report(self):
        print('========Total Forces Report========')
        if not self.forcesavailable:
            self.calcforces(echo=False)
        print('%10s %10s %10s %10s %10s %10s | %10s %10s' % ('CX', 'CY', 'CZ', 'Cl', 'Cm', 'Cn', 'CL', 'CD'))
        print('%10f %10f %10f %10f %10f %10f | %10f %10f' % (self.CX, self.CY, self.CZ, self.Cl, self.Cm, self.Cn, self.CL, self.CD))
        print('===================================')