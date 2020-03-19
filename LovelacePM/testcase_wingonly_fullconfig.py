''' EQUIVALENT TO FULLCONFIG TEST CASE, BUT INCLUDING ONLY ITS LIFTING BODIES'''
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
from aerodynamic_output import *

a=-radians(0.0)
Uinf=0.01
b=1.2
taper=0.5
croot=0.2
fuselage_aftlen=0.4
fuselage_rearlen=0.6
fuselage_width=0.1

sld=Solid()

sect1=wing_section(CA_position=np.array([0.0, -b/2, 0.0]), c=croot*taper, xdisc=20)
sect2=wing_section(CA_position=np.array([0.0, 0.0, 0.0]), c=croot, xdisc=20)

wingquad_left=wing_quadrant(sld, sect1=sect1, sect2=sect2)

sect3=wing_section(CA_position=np.array([0.0, b/2, 0.0]), c=croot*taper, xdisc=20)

wingquad_right=wing_quadrant(sld, sect1=sect2, sect2=sect3)

wng=wing(sld, wingquads=[wingquad_left, wingquad_right])

wng.patchcompose(ydisc=40, ystrategy=lambda x: x)#lambda x: (np.sin(pi*x-pi/2)+1)/2)

stabsect1=wing_section(CA_position=np.array([fuselage_rearlen, -0.3, 0.3]), inverse=True, c=0.05, incidence=-radians(0.0), xdisc=30)
stabsect2=wing_section(CA_position=np.array([fuselage_rearlen, 0.0, 0.3]), inverse=True, c=0.1, incidence=-radians(0.0), xdisc=30)
stabsect3=wing_section(CA_position=np.array([fuselage_rearlen, 0.3, 0.3]), inverse=True, c=0.05, incidence=-radians(0.0), xdisc=30)

horz_emp_left=wing_quadrant(sld, sect1=stabsect1, sect2=stabsect2)
horz_emp_right=wing_quadrant(sld, sect1=stabsect2, sect2=stabsect3)
horz_emp=wing(sld, wingquads=[horz_emp_left, horz_emp_right])

horz_emp.patchcompose(ydisc=24)

Sref=croot*(1+taper)*(b+fuselage_width)/2
acft=aircraft(sld, elems=[wng, horz_emp], Sref=Sref, bref=b)

sld.plotgeometry(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0])
#sld.plotnormals(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], factor=0.1)
#sld.plotnormals(xlim=[-0.1, 0.3], ylim=[0.5, 0.7], zlim=[-0.2, 0.2], factor=0.1)
#sld.genwakepanels(wakecombs=wingquad_left.wakecombs+wingquad_right.wakecombs, wakeinds=[[0, 0]], a=a)
acft.edit_parameters(par='a', val=a)
acft.addwake()
acft.eulersolve()
acft.forces_report()
acft.stabreport()
sld.plotgeometry(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0])
#sld.plotgeometry(xlim=[-0.1, 0.3], ylim=[0.5, 0.7], zlim=[-0.2, 0.2])
print(sld.npanels)
plot_Cps(sld, elems=[wng])
plot_Cls(sld, wings=[wng])
plot_Cds(sld, wings=[wng])
plot_Cms(sld, wings=[wng])
plot_gammas(sld, wings=[wng])
plot_Cps(sld, elems=[horz_emp])
plot_Cls(sld, wings=[horz_emp])
plot_Cds(sld, wings=[horz_emp])
plot_Cms(sld, wings=[horz_emp])
plot_gammas(sld, wings=[horz_emp])