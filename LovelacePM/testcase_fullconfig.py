import numpy as np
import numpy.linalg as lg
import scipy.linalg as slg
import scipy.linalg.lapack as slapack
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import time as tm

from paneller import *
from utils import *
from wing import *
from body import *
from aircraft import *
from xfoil_visc import *
from aerodynamic_output import *

ordir=os.getcwd()
os.chdir(os.path.dirname(os.path.abspath(__file__)))

a=0.0
Uinf=15
rho=1.225
mu=1.72*10e-5
b=1.2
taper=0.5
croot=0.2
fuselage_aftlen=0.4
fuselage_rearlen=0.6
fuselage_width=0.1

sld=Solid()

#generating xfoil corrected polars
n4412_polar=polar_correction(name='n4412')

fuselage=standard_body(sld, defsect=circdefsect, body_thdisc=50, body_width=fuselage_width, nose_loc=np.array([-fuselage_aftlen, 0.0, -0.01]))

sect1=wing_section(afl='n4412', CA_position=np.array([0.0, -b/2, 0.0]), c=croot*taper, xdisc=20, correction=n4412_polar, Re=Uinf*rho*croot*taper/mu, closed=True)
sect2=wing_section(afl='n4412', CA_position=np.array([0.0, -fuselage_width/2, 0.0]), c=croot, xdisc=20, correction=n4412_polar, Re=Uinf*rho*croot/mu)

wingquad_left=wing_quadrant(sld, sect1=sect1, sect2=sect2)

sect3=wing_section(afl='n4412', CA_position=np.array([0.0, fuselage_width/2, 0.0]), c=croot, xdisc=20, correction=n4412_polar, Re=Uinf*rho*croot/mu)
sect4=wing_section(afl='n4412', CA_position=np.array([0.0, b/2, 0.0]), c=croot*taper, xdisc=20, correction=n4412_polar, Re=Uinf*rho*croot*taper/mu, closed=True)

wingquad_right=wing_quadrant(sld, sect1=sect3, sect2=sect4)

wing_left=wing(sld, wingquads=[wingquad_left])
wing_right=wing(sld, wingquads=[wingquad_right])

wing_right.trim_bybody(fuselage, sectside=1)
wing_left.trim_bybody(fuselage, sectside=2)

stabsect1=wing_section(CA_position=np.array([fuselage_rearlen, -0.3, 0.3]), inverse=True, c=0.05, incidence=-radians(0.0), xdisc=30, \
    correction=n4412_polar, Re=Uinf*rho*0.05/mu, closed=True)
stabsect2=wing_section(CA_position=np.array([fuselage_rearlen, 0.0, 0.3]), inverse=True, c=0.1, incidence=-radians(0.0), xdisc=30, \
    correction=n4412_polar, Re=Uinf*rho*0.1/mu, closed=True)
stabsect3=wing_section(CA_position=np.array([fuselage_rearlen, 0.3, 0.3]), inverse=True, c=0.05, incidence=-radians(0.0), xdisc=30, \
    correction=n4412_polar, Re=Uinf*rho*0.05/mu, closed=True)

horz_emp_left=wing_quadrant(sld, sect1=stabsect1, sect2=stabsect2)
horz_emp_right=wing_quadrant(sld, sect1=stabsect2, sect2=stabsect3)
horz_emp=wing(sld, wingquads=[horz_emp_left, horz_emp_right])

Sref=croot*(1+taper)*(b+fuselage_width)/2
acft=aircraft(sld, elems=[wing_left, wing_right, fuselage, horz_emp], Sref=Sref, bref=b)
acft.plot_input()

horz_emp.patchcompose(ydisc=24)
wing_right.patchcompose(ydisc=20, ystrategy=lambda x: x)#lambda x: (np.sin(pi*x-pi/2)+1)/2)
wing_left.patchcompose(ydisc=20, ystrategy=lambda x: x)#lambda x: (np.sin(pi*x-pi/2)+1)/2)

fuselage.patchcompose(leftqueue=[wing_left], rightqueue=[wing_right], xdisc=60, \
    thdisc_upleft=5, thdisc_downleft=5, thdisc_downright=5, thdisc_upright=5)


#sld.plotnormals(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], factor=0.1)
#sld.plotnormals(xlim=[-0.1, 0.3], ylim=[0.5, 0.7], zlim=[-0.2, 0.2], factor=0.1)
#sld.genwakepanels(wakecombs=wingquad_left.wakecombs+wingquad_right.wakecombs, wakeinds=[[0, 0]], a=a)
acft.edit_parameters({'a':a, 'Uinf':Uinf})
acft.addwake(offset=10.0, wakedisc=30, strategy=lambda x: x)
sld.plotgeometry(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0])
acft.bodies_eqflatplate_apply(rho=rho, mu=mu)
acft.eulersolve(wakeiter=1)
acft.forces_report()
acft.stabreport()
acft.balance()
#    return acft.CD, acft.dCD, acft.CL, acft.dCL

'''CLs=[]
CDs=[]
CLs_visc=[]
CDs_visc=[]
for a in np.arange(0.0, 10.0, 1.0):
    CD, dCD, CL, dCL=polar_data(a)
    CLs+=[CL]
    CLs_visc+=[CL+dCL]
    CDs+=[CD]
    CDs_visc+=[CD+dCD]
plt.plot(CLs_visc, CDs_visc, label='viscous')
plt.plot(CLs, CDs, label='inviscid')
plt.legend()
plt.show()'''
acft.plotgeometry(velfield=False)
#sld.plotgeometry(xlim=[-0.1, 0.3], ylim=[0.5, 0.7], zlim=[-0.2, 0.2])
print('npanels: ', sld.npanels)
plot_Cps(sld, elems=[wing_left, wing_right])
plot_Cls(sld, wings=[wing_left, wing_right])
plot_Cds(sld, wings=[wing_left, wing_right])
plot_Cms(sld, wings=[wing_left, wing_right])
plot_gammas(sld, wings=[wing_left, wing_right])
plot_Cps(sld, elems=[horz_emp])
plot_Cls(sld, wings=[horz_emp])
plot_Cds(sld, wings=[horz_emp])
plot_Cms(sld, wings=[horz_emp])
plot_gammas(sld, wings=[horz_emp])

os.chdir(ordir)