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

import fdyn
import toolkit
from utils import *

from paneller import *
from body import *
from wing import *

class aircraft: #class to ease certain case studies for full aircraft
    #generate aircraft BEFORE patchcompose functions
    def parameter_report(self):
        print('Freestream parameters:')
        print('%10s %10s %10s %10s %10s %10s %10s' % ('a', 'b', 'p', 'q', 'r', 'Uinf', 'M'))
        print('%10f %10f %10f %10f %10f %10f %10f' % (degrees(self.a), degrees(self.b), self.p, self.q, self.r, self.Uinf, self.M))
        print('Environment parameters:')
        print('%10s %10s %10s' % ('rho', 'mu', 'g'))
        print('%10f %10f %10f' % (self.rho, self.mu, self.g))
        print('Geometry parameters:')
        print('%10s %10s %10s %10s' % ('Sref', 'cref', 'bref', 'AR'))
        print('%10f %10f %10f %10f' % (self.Sref, self.cref, self.bref, self.AR))
        if self.hascontrol():
            print('Control parameters:')
            controlnames_str=''
            controlvals_str=''
            for cname in self.controlset:
                controlnames_str+='%10s' % cname
                controlvals_str+='%10f' % degrees(self.controlset[cname].state)
            print(controlnames_str)
            print(controlvals_str)
        if not self.massavailable:
            print('Mass parameters not yet present')
        else:
            print('Mass parameters:')
            print('CG=(%10f %10f %10f)' % (self.CG[0], self.CG[1], self.CG[2]))
            print('%10s %10s %10s %10s %10s %10s %10s' % ('m', 'Ixx', 'Iyy', 'Izz', 'Ixy', 'Ixz', 'Iyz'))
            print('%10f %10f %10f %10f %10f %10f %10f' % (self.m, self.Inertia[0, 0], self.Inertia[1, 1], self.Inertia[2, 2], \
                self.Inertia[0, 1], self.Inertia[0, 2], self.Inertia[1, 2]))
    def addmass(self, m=0.0, Ixx=0.0, Iyy=0.0, Izz=0.0, Ixy=0.0, Ixz=0.0, Iyz=0.0, echo=True): #sets mass variables for the aircraft
        self.m=m
        self.Inertia=np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
        self.massavailable=True
        if echo:
            self.parameter_report()
    def dynamic_simulation(self, nstep=2e5, dt=1e-4, start_point=np.array([0.0, 0.0, 0.0]), \
        perturbations={}, trim=True, visc=True, onboard=np.zeros(3)):
        if not self.massavailable:
            print('WARNING: dynamic analysis requested without any mass data. Returning empty arrays')
            return np.array([]), np.array([]), np.array([]), np.array([])
        perturbations=np.array([start_point[0], start_point[1], start_point[2], \
            perturbations['phi'], perturbations['theta'], perturbations['psi'], \
                perturbations['u'], perturbations['v'], perturbations['w'], \
                    perturbations['p'], perturbations['q'], perturbations['q']])
        coeffs=np.zeros(6)
        derivs=np.zeros((5, 6))
        if not trim:
            coeffs=np.array([-self.CX, self.CY, -self.CZ, -self.Cl, self.Cm, -self.Cn]) #note that we're using stability coordinates here
            if visc:
                coeffs+=np.array([-self.dCX, self.dCY, -self.dCZ, -self.dCl, self.dCm, -self.dCn])
        i=0
        for par in self.stabderivative_dict:
            derivs[i, :]=np.array([-self.stabderivative_dict[par]['dCX'], self.stabderivative_dict[par]['dCY'], -self.stabderivative_dict[par]['dCZ'], \
                -self.stabderivative_dict[par]['dCl'], self.stabderivative_dict[par]['dCm'], -self.stabderivative_dict[par]['dCn']])
            i+=1
        external_history, alpha_history, beta_history, euler_history=fdyn.tstep_solve(nstep, dt, self.rho, self.Uinf, self.g, \
            self.Sref, self.cref, self.bref, perturbations, onboard, self.Inertia, self.m, \
                coeffs, derivs)
        return external_history, alpha_history, beta_history, euler_history
    def balance(self, SM=0.1, echo=True): #balancing the aircraft for given static margin
        if not self.stabavailable:
            print('WARNING: stability derivatives not available for balancing. Performing in place')
            self.calcstab(echo=echo)
        loc_SM=-self.stabderivative_dict['a']['dCm']/self.stabderivative_dict['a']['dCL']
        self.transp_byvec(-(SM-loc_SM)*np.array([1.0, 0.0, 0.0])*self.cref)
        new_SM=-self.stabderivative_dict['a']['dCm']/self.stabderivative_dict['a']['dCL']
        if echo:
            print('==========Balancing========')
            print('%10s %10s' % ('SM init', 'SM final'))
            print('%10f %10f' % (loc_SM, new_SM))
            self.forces_report()
            if self.stabavailable:
                self.stabreport()
            print('===========================')
    def transp_byvec(self, vec): #transport pre-calculated forces by a 'vec' deallocation in CG
        mom_transp=np.cross(vec, np.array([self.CX, self.CY, self.CZ]))
        self.Cl-=mom_transp[0]/self.bref
        self.Cm-=mom_transp[1]/self.cref
        self.Cn-=mom_transp[2]/self.bref
        if self.stabavailable:
            pars=['a', 'b', 'p', 'q', 'r']
            for par in pars:
                dmom_transp=np.cross(vec, np.array([self.stabderivative_dict[par]['dCX'], self.stabderivative_dict[par]['dCY'], \
                    self.stabderivative_dict[par]['dCZ']]))
                self.stabderivative_dict[par]['dCl']-=dmom_transp[0]/self.bref
                self.stabderivative_dict[par]['dCm']-=dmom_transp[1]/self.cref
                self.stabderivative_dict[par]['dCn']-=dmom_transp[2]/self.bref
        if self.hascorrections:
            mom_transp=np.cross(vec, np.array([self.dCX, self.dCY, self.dCZ]))
            self.dCl-=mom_transp[0]/self.bref
            self.dCm-=mom_transp[1]/self.cref
            self.dCn-=mom_transp[2]/self.bref
    def transp_to_cg(self, CG):
        self.transp_byvec(-self.CG)
        self.transp_byvec(CG)
        self.CG=CG
    def __init__(self, sld, elems=[], Sref=0.0, cref=0.0, bref=0.0, echo=True, CG=np.array([0.0, 0.0, 0.0])):
        self.wings=[]
        self.bodies=[]
        for e in elems:
            e.set_aircraft(self)
            if type(e)==wing:
                self.wings+=[e]
            elif type(e)==body:
                self.bodies+=[e]
        #BE VERY AWARE OF THESE WARNINGS WHEN USINGG FUSELAGE-SEPARATED SEMI-WINGS AS WING OBJECTS
        if len(self.wings)!=0:
            S, mac, b=self.wings[0].calc_reference()
        else:
            S=1.0; mac=1.0; b=1.0
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
        self.AR=(self.bref**2)/self.Sref

        self.a=0.0
        self.b=0.0
        self.p=0.0
        self.q=0.0
        self.r=0.0
        self.Uinf=1.0
        self.M=0.0
        self.rho=1.225
        self.mu=1.72e-5
        self.g=9.81

        self.CG=CG
        
        #defining available controls
        self.controlset={}
        for wng in self.wings:
            for wngqd in wng.wingquads:
                if wngqd.hascontrol():
                    for cname in wngqd.controls:
                        if not cname in self.controlset:
                            self.controlset[cname]=control_DOF()
                            wngqd.controls[cname].DOF=self.controlset[cname]

        self.sld=sld
        self.forcesavailable=False
        self.stabavailable=False
        self.massavailable=False
        self.haseqflatplate=False
        self.hascorrections=any([wng.hascorrection() for wng in self.wings])

        if echo:
            self.parameter_report()

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
    def addwake(self, offset=1000.0, wakedisc=1, strategy=lambda x: ((np.exp(x)-1.0)/(exp(1)-1.0))**2):
        #once the wake has been added, one can end geometry pre-processing of panels in question:
        self.sld.end_preprocess()
        for wng in self.wings:
            wng.genwakepanels(offset=offset, a=self.a, b=self.b, wakedisc=wakedisc, strategy=strategy)
    def calcforces(self, echo=True):
        self.CX=sum([self.sld.forces[i][0] for i in range(self.sld.npanels)])/self.Sref
        self.CY=sum([self.sld.forces[i][1] for i in range(self.sld.npanels)])/self.Sref
        self.CZ=sum([self.sld.forces[i][2] for i in range(self.sld.npanels)])/self.Sref
        u=np.array([cos(self.a)*cos(self.b), sin(self.b)*cos(self.a), sin(self.a)])
        v=np.array([-sin(self.a)*cos(self.b), -sin(self.b)*sin(self.a), cos(self.a)])
        adforces=np.array([self.CX, self.CY, self.CZ])
        self.CL=v@adforces
        self.CD=u@adforces
        self.Cl=sum([self.sld.moments[i][0] for i in range(self.sld.npanels)])/(self.Sref*self.bref)
        self.Cm=sum([self.sld.moments[i][1] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
        self.Cn=sum([self.sld.moments[i][2] for i in range(self.sld.npanels)])/(self.Sref*self.bref)

        if self.hascorrections:
            self.dCX=0.0; self.dCY=0.0; self.dCZ=0.0; self.dCL=0.0; self.dCD=0.0; self.dCl=0.0; self.dCm=0.0; self.dCn=0.0
            for wng in self.wings:
                wng.calc_coefs(alpha=self.a, axis=wng.axis)
                if wng.hascorrection():
                    dCX_wng, dCY_wng, dCZ_wng, dCl_wng, dCm_wng, dCn_wng=wng.calc_corrected_forces()
                    self.dCX+=dCX_wng
                    self.dCY+=dCY_wng
                    self.dCZ+=dCZ_wng
                    self.dCl+=dCl_wng
                    self.dCm+=dCm_wng
                    self.dCn+=dCn_wng
            adforces=np.array([self.dCX, self.dCY, self.dCZ])
            self.dCL=adforces@v
            self.dCD=adforces@u
        
        self.transp_byvec(self.CG)

        if echo:
            self.forces_report()
        self.forcesavailable=True
    def freestream_derivative(self, par='a'): #definition of derivative of freestream velocity vectors u and v (normal vectors)
        if par=='a':
            du=np.array([-sin(self.a)*cos(self.b), -sin(self.b)*sin(self.a), cos(self.a)])
            dv=np.array([-cos(self.a)*cos(self.b), -sin(self.b)*cos(self.a), -sin(self.a)])
        elif par=='b':
            du=np.array([-cos(self.a)*sin(self.b), cos(self.b)*cos(self.a), 0.0])
            dv=np.array([sin(self.a)*sin(self.b), -cos(self.b)*sin(self.a), 0.0])
        else:
            du=np.zeros(3)
            dv=np.zeros(3)
        return du, dv
    def calcstab(self, echo=True):
        #calculate variation on Cps in respect to several parameters
        pars=['a', 'b', 'p', 'q', 'r']
        if self.forcesavailable:
            self.transp_byvec(-self.CG)
        self.stabderivative_dict={}
        orforces=np.array([self.CX, self.CY, self.CZ])
        for p in pars:
            dCps=self.sld.calc_derivative(self.Uinf, a=self.a, b=self.b, p=self.p*2*self.Uinf/self.bref, q=self.q*2*self.Uinf/self.cref, r=self.r*2*self.Uinf/self.bref, par=p)
            forces=[-dCps[i]*self.sld.panels[i].S*self.sld.panels[i].nvector for i in range(self.sld.npanels)]
            moments=[np.cross(self.sld.panels[i].colpoint, forces[i]) for i in range(self.sld.npanels)]
            dCX=sum([forces[i][0] for i in range(self.sld.npanels)])/self.Sref
            dCY=sum([forces[i][1] for i in range(self.sld.npanels)])/self.Sref
            dCZ=sum([forces[i][2] for i in range(self.sld.npanels)])/self.Sref
            u=np.array([cos(self.a)*cos(self.b), sin(self.b)*cos(self.a), sin(self.a)])
            v=np.array([-sin(self.a)*cos(self.b), -sin(self.b)*sin(self.a), cos(self.a)])
            adforces=np.array([dCX, dCY, dCZ])
            du, dv=self.freestream_derivative(par=p)
            dCL=v@adforces+orforces@dv
            dCD=u@adforces+orforces@du
            dCl=sum([moments[i][0] for i in range(self.sld.npanels)])/(self.Sref*self.bref)
            dCm=sum([moments[i][1] for i in range(self.sld.npanels)])/(self.Sref*self.cref)
            dCn=sum([moments[i][2] for i in range(self.sld.npanels)])/(self.Sref*self.bref)
            self.stabderivative_dict[p]={'dCX':dCX, 'dCY':dCY, 'dCZ':dCZ, 'dCL':dCL, 'dCD':dCD, 'dCl':dCl, 'dCm':dCm, 'dCn':dCn}
        self.stabavailable=True
        self.transp_byvec(self.CG)
        if echo:
            self.stabreport()
    def stabreport(self): #report calculated stability derivatives, or calculate them in place
        if not self.stabavailable:
            self.calcstab(echo=False)
        print('Stability derivatives')
        print('%5s | %10s %10s %10s %10s %10s %10s %10s %10s' % ('par', 'dCX', 'dCY', 'dCZ', 'dCL', 'dCD', 'dCl', 'dCm', 'dCn'))
        for p in self.stabderivative_dict:
            print('%5s | %10f %10f %10f %10f %10f %10f %10f %10f' % (p, self.stabderivative_dict[p]['dCX'], self.stabderivative_dict[p]['dCY'], \
                self.stabderivative_dict[p]['dCZ'], self.stabderivative_dict[p]['dCL'], self.stabderivative_dict[p]['dCD'], -self.stabderivative_dict[p]['dCl'], \
                    self.stabderivative_dict[p]['dCm'], self.stabderivative_dict[p]['dCn']))
    def edit_parameters(self, pardict, echo=True):
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
            elif par=='M':
                self.M=val
            elif par=='rho':
                self.rho=val
            elif par=='mu':
                self.mu=val
            elif par=='g':
                self.g=val
            else:
                self.controlset[par].state=radians(val)
        #reset result readiness
        self.stabavailable=False
        self.forcesavailable=False
        self.haseqflatplate=False
        if echo:
            self.parameter_report()
    def bodies_eqflatplate_apply(self, rho=1.225, mu=1.72*10e-5, turbulent_criterion=Re2e5, Cf_l_rule=Blausius_Cf_l, Cf_t_rule=Prandtl_1_7th):
        self.haseqflatplate=True
        for bdy in self.bodies:
            bdy.apply_eqflatplate(rho=rho, mu=mu, turbulent_criterion=turbulent_criterion, Cf_l_rule=Cf_l_rule, Cf_t_rule=Cf_t_rule, Uinf=self.Uinf)
    def eulersolve(self, echo=True, damper=0.0, tolerance=1e-5, wakeiter=0):
        self.sld.eulersolve(target=np.array([]), a=self.a, b=self.b, p=self.p*2*self.Uinf/self.bref, q=self.q*2*self.Uinf/self.cref, \
            r=self.r*2*self.Uinf/self.bref, damper=damper, Uinf=self.Uinf, echo=True, tolerance=tolerance, wakeiter=wakeiter, beta=sqrt(1.0-self.M**2))
    def forces_report(self):
        print('========Total Forces Report========')
        if self.hascorrections:
            if self.haseqflatplate:
                if self.haseqflatplate:
                    str_Cf='+local Cf estimation'
            else:
                str_Cf=''
            print('Inviscid'+str_Cf+' case: ')
        if not self.forcesavailable:
            self.calcforces(echo=False)
        print('%10s %10s %10s %10s %10s %10s | %10s %10s' % ('CX', 'CY', 'CZ', 'Cl', 'Cm', 'Cn', 'CL', 'CD'))
        print('%10f %10f %10f %10f %10f %10f | %10f %10f' % (self.CX, self.CY, self.CZ, -self.Cl, self.Cm, -self.Cn, self.CL, self.CD))
        if self.hascorrections:
            print('After viscous external polar corrections: ')
            print('%10s %10s %10s %10s %10s %10s | %10s %10s' % ('CX', 'CY', 'CZ', 'Cl', 'Cm', 'Cn', 'CL', 'CD'))
            print('%10f %10f %10f %10f %10f %10f | %10f %10f' % (self.CX+self.dCX, self.CY+self.dCY, self.CZ+self.dCZ, \
                -(self.Cl+self.dCl), self.Cm+self.dCm, -(self.Cn+self.dCn), self.CL+self.dCL, self.CD+self.dCD))
        print('===================================')
    def plot_input(self, xlim=[], ylim=[], zlim=[]):
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        for elem in self.bodies+self.wings:
            elem.plot_input(ax=ax, fig=fig)
        if len(xlim)!=0:
            ax.set_xlim3d(xlim[0], xlim[1])
        else:
            ax.set_xlim3d(-self.plotlim, self.plotlim)
        if len(ylim)!=0:
            ax.set_ylim3d(ylim[0], ylim[1])
        else:
            ax.set_ylim3d(-self.plotlim, self.plotlim)
        if len(zlim)!=0:
            ax.set_zlim3d(zlim[0], zlim[1])
        else:
            ax.set_zlim3d(-self.plotlim, self.plotlim)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()
    def plotgeometry(self, xlim=[], ylim=[], zlim=[], factor=1.0, velfield=True):
        if len(xlim)==0:
            xlim=[-self.plotlim, self.plotlim]
        if len(ylim)==0:
            ylim=[-self.plotlim, self.plotlim]
        if len(zlim)==0:
            zlim=[-self.plotlim, self.plotlim]
        self.sld.plotgeometry(xlim=xlim, ylim=ylim, zlim=zlim, factor=factor, velfield=velfield)
