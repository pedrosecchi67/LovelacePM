import numpy as np
from math import *
import numpy.linalg as lg

#module to replace fdyn.f90

#solve aircraft positions as they vary in time
def tstep_solve(nstep, dt, rho, U0, g, Sref, cref, bref, perturbation, onboard, inertia, m, \
    coeffs, derivs):
    external_history=np.zeros((nstep, 3))
    euler_history=np.zeros((nstep, 3))
    alpha_history=np.zeros(nstep)
    beta_history=np.zeros(nstep)
    time_history=np.linspace(0.0, dt*nstep, nstep)
    a=0.0
    b=0.0
    Uinf=0.0
    forces=np.zeros(3)
    moments=np.zeros(3)
    mtemp=np.zeros(3)
    #for pertutbation, here we use the notation employed by Drela in Flight Viechle Aerodynamics
    #equations are also directly obtained from there
    #notation for coefs: (CX, CY, CZ, Cl, Cm, Cn)
    #notation for derivs: (CX, CY, CZ, Cl, Cm, Cn), first row: a derivatives; second, b derivatives; third, p; 4, q; 5, r
    xe=perturbation[0]
    ye=perturbation[1]
    ze=perturbation[2]
    phi=perturbation[3]
    theta=perturbation[4]
    psi=perturbation[5]
    u=perturbation[6]+U0
    v=perturbation[7]
    w=perturbation[8]
    p=perturbation[9]*2*U0/bref
    q=perturbation[10]*2*U0/cref
    r=perturbation[11]*2*U0/bref #denormalizing angular velocities
    invI=lg.inv(inertia)
    for i in range(nstep):
        Uinf=sqrt(u**2+v**2+w**2)
        a=atan(w/u)
        b=atan(v/sqrt(u**2+w**2))
        for j in range(3):
            forces[j]=(coeffs[j]+derivs[0, j]*a+derivs[1, j]*b+derivs[2, j]*p*bref/(2*Uinf)+\
                derivs[3, j]*q*cref/(2*Uinf)+derivs[4, j]*r*bref/(2*Uinf))*Sref*Uinf**2*rho/2
            moments[j]=(coeffs[j+3]+derivs[0, j+3]*a+derivs[1, j+3]*b+derivs[2, j+3]*p*bref/(2*Uinf)+\
                derivs[3, j+3]*q*cref/(2*Uinf)+derivs[4, j+3]*r*bref/(2*Uinf))*Sref*(cref if j==1 else bref)*Uinf**2*rho/2
        u=u+dt*(forces[0]/m-g*sin(theta)-q*w+r*v)
        v=v+dt*(forces[1]/m+g*sin(phi)*cos(theta)-r*u+p*w)
        w=w+dt*(forces[2]/m+g*cos(phi)*cos(theta)-p*v+q*u)
        mtemp[0]=moments[0]-(inertia[2, 2]-inertia[1, 1])*q*r-inertia[1, 2]*(q**2-r**2)-\
            inertia[0, 2]*p*q+inertia[0, 1]*p*r-onboard[2]*q+onboard[1]*r
        mtemp[1]=moments[1]-(inertia[0, 0]-inertia[2, 2])*r*p-inertia[0, 2]*(r**2-p**2)-\
            inertia[0, 1]*q*r+inertia[1, 2]*q*p-onboard[0]*r+onboard[2]*p
        mtemp[2]=moments[2]-(inertia[2, 2]-inertia[0, 0])*p*q-inertia[0, 1]*(p**2-q**2)-\
            inertia[1, 2]*r*p+inertia[0, 2]*r*q-onboard[1]*p+onboard[0]*q
        mtemp=invI@mtemp
        p+=dt*mtemp[0]
        q+=dt*mtemp[1]
        r+=dt*mtemp[2]
        phi=phi+dt*(p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta))
        theta=theta+dt*(q*cos(phi)-r*sin(phi))
        psi=psi+dt*(q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta))
        alpha_history[i]=a
        beta_history[i]=b
        euler_history[i, 0]=phi
        euler_history[i, 1]=theta
        euler_history[i, 2]=psi
        external_history[i, 0]=xe
        external_history[i, 1]=ye
        external_history[i, 2]=ze
    return external_history, alpha_history, beta_history, euler_history, time_history