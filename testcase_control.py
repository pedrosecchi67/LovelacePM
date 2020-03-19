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
from control import *
from aircraft import *
from aerodynamic_output import *

sld=Solid()

b=1.0
c=0.1
Uinf=0.01

sect1=wing_section(afl='n4412', c=c, CA_position=np.array([0.0, -b/2, 0.0]), xdisc=30)
sect2=wing_section(afl='n4412', c=c, CA_position=np.array([0.0, 0.0, 0.0]), xdisc=30)
sect3=wing_section(afl='n4412', c=c, CA_position=np.array([0.0, b/2, 0.0]), xdisc=30)
wngqd1=wing_quadrant(sld, sect1=sect1, sect2=sect2)
wngqd2=wing_quadrant(sld, sect1=sect2, sect2=sect3, control_names=['aileron'], control_axpercs_x=[[0.75, 0.75]], control_axpercs_thickness=[[0.5, 0.5]])
wng=wing(sld, wingquads=[wngqd1, wngqd2])
acft=aircraft(sld, elems=[wng])
acft.edit_parameters({'Uinf':0.01, 'aileron':1})
wng.patchcompose(ydisc=70)
sld.plotgeometry(xlim=[-c/2, c*3/4], ylim=[-c/2, c/2], zlim=[-c/2, c/2])
acft.addwake()
acft.eulersolve()
plot_Cps(sld, elems=[wng])
plot_Cls(sld, wings=[wng])
sld.plotgeometry(xlim=[-b/2, b/2], ylim=[-b/2, b/2], zlim=[-b/2, b/2])
plot_Cps(sld, elems=[wng])
plot_Cls(sld, wings=[wng])
plot_Cds(sld, wings=[wng])
plot_Cms(sld, wings=[wng])