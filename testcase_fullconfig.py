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

a=radians(5.0)
Uinf=0.01

sld=Solid()

fig=plt.figure()
ax=plt.axes(projection='3d')

fuselage=standard_body(sld, body_thdisc=40, nose_loc=np.array([-0.4, 0.0, 0.0]))
fuselage.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

sect1=wing_section(CA_position=np.array([0.0, -0.6, 0.0]), c=0.1)
sect2=wing_section(CA_position=np.array([0.0, -0.05, 0.0]), c=0.2)

wingquad_left=wing_quadrant(sld, sect1=sect1, sect2=sect2)
wingquad_left.trim_bybody(fuselage, sectside=2)
wingquad_left.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

wingquad_left.patchcompose()

sect3=wing_section(CA_position=np.array([0.0, 0.05, 0.0]), c=0.2)
sect4=wing_section(CA_position=np.array([0.0, 0.6, 0.0]), c=0.1)

wingquad_right=wing_quadrant(sld, sect1=sect3, sect2=sect4)
wingquad_right.trim_bybody(fuselage, sectside=1)
wingquad_right.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

plt.show()

wingquad_right.patchcompose()

fuselage.patchcompose(leftqueue=[wingquad_left], rightqueue=[wingquad_right], xdisc=30, \
    thdisc_upleft=5, thdisc_downleft=5, thdisc_downright=5, thdisc_upright=5)

sld.end_preprocess()
#sld.plotnormals(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], factor=0.01)
sld.genwakepanels(wakecombs=wingquad_left.wakecombs+wingquad_right.wakecombs, wakeinds=[[0, 0]], a=a)
sld.eulersolve(Uinf=Uinf, a=a)
#sld.plotgeometry(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0])
sld.plotgeometry(xlim=[-0.1, 0.3], ylim=[-0.2, 0.2], zlim=[-0.2, 0.2])
print(sld.npanels)