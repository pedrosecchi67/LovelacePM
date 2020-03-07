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
from aerodynamic_output import *

a=radians(0.0)
Uinf=0.01

sld=Solid()

fig=plt.figure()
ax=plt.axes(projection='3d')

fuselage=standard_body(sld, body_thdisc=100, nose_loc=np.array([-0.4, 0.0, 0.0]))
fuselage.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

sect1=wing_section(CA_position=np.array([0.0, -0.6, 0.0]), c=0.1, xdisc=20)
sect2=wing_section(CA_position=np.array([0.0, -0.05, 0.0]), c=0.2, xdisc=20)

wingquad_left=wing_quadrant(sld, sect1=sect1, sect2=sect2)

sect3=wing_section(CA_position=np.array([0.0, 0.05, 0.0]), c=0.2, xdisc=20)
sect4=wing_section(CA_position=np.array([0.0, 0.6, 0.0]), c=0.1, xdisc=20)

wingquad_right=wing_quadrant(sld, sect1=sect3, sect2=sect4)

wing_left=wing(wingquads=[wingquad_left])
wing_right=wing(wingquads=[wingquad_right])

wing_right.patchcompose(ydisc=20, ystrategy=lambda x: x)#lambda x: (np.sin(pi*x-pi/2)+1)/2)
wing_left.patchcompose(ydisc=20, ystrategy=lambda x: x)#lambda x: (np.sin(pi*x-pi/2)+1)/2)

wing_right.trim_bybody(fuselage, sectside=1)
wing_left.trim_bybody(fuselage, sectside=2)

stabsect1=wing_section(CA_position=np.array([0.6, -0.3, 0.3]), inverse=True, c=0.05, incidence=-radians(2.0), xdisc=30)
stabsect2=wing_section(CA_position=np.array([0.6, 0.0, 0.3]), inverse=True, c=0.1, incidence=-radians(2.0), xdisc=30)
stabsect3=wing_section(CA_position=np.array([0.6, 0.3, 0.3]), inverse=True, c=0.05, incidence=-radians(2.0), xdisc=30)

horz_emp_left=wing_quadrant(sld, sect1=stabsect1, sect2=stabsect2)
horz_emp_right=wing_quadrant(sld, sect1=stabsect2, sect2=stabsect3)
horz_emp=wing(wingquads=[horz_emp_left, horz_emp_right])

horz_emp.patchcompose(ydisc=24)

wing_left.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

wing_right.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

fuselage.patchcompose(leftqueue=[wing_left], rightqueue=[wing_right], xdisc=60, \
    thdisc_upleft=5, thdisc_downleft=5, thdisc_downright=5, thdisc_upright=5)
    
horz_emp.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

plt.show()

sld.end_preprocess()

wing_left.genwakepanels(a=a)
wing_right.genwakepanels(a=a)
horz_emp.genwakepanels(a=a)

sld.plotgeometry(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0])
#sld.plotnormals(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], factor=0.1)
#sld.plotnormals(xlim=[-0.1, 0.3], ylim=[0.5, 0.7], zlim=[-0.2, 0.2], factor=0.1)
#sld.genwakepanels(wakecombs=wingquad_left.wakecombs+wingquad_right.wakecombs, wakeinds=[[0, 0]], a=a)
sld.eulersolve(Uinf=Uinf, a=a)
sld.plotgeometry(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0])
#sld.plotgeometry(xlim=[-0.1, 0.3], ylim=[0.5, 0.7], zlim=[-0.2, 0.2])
print(sld.npanels)
'''plot_Cps(sld, elems=[wing_left, wing_right])
plot_Cls(sld, wings=[wing_left, wing_right])
plot_Cds(sld, wings=[wing_left, wing_right])
plot_Cms(sld, wings=[wing_left, wing_right])
plot_gammas(sld, wings=[wing_left, wing_right])
plot_Cps(sld, elems=[horz_emp])
plot_Cls(sld, wings=[horz_emp])
plot_Cds(sld, wings=[horz_emp])
plot_Cms(sld, wings=[horz_emp])
plot_gammas(sld, wings=[horz_emp])'''