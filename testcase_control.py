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
flap=control(p0=np.array([c/2, -b/2, 0.0]), p1=np.array([c/2, b/2, 0.0]))

sect1=wing_section(afl='n4412', c=c, CA_position=np.array([0.0, -b/2, 0.0]), xdisc=30)
sect2=wing_section(afl='n4412', c=c, CA_position=np.array([0.0, b/2, 0.0]), xdisc=30)
wngqd=wing_quadrant(sld, sect1=sect1, sect2=sect2, control_names=['aileron'], control_axpercs_x=[[0.75, 0.75]], control_axpercs_thickness=[[0.5, 0.5]])
wng=wing(wingquads=[wngqd])
acft=aircraft(sld, elems=[wng])
wng.patchcompose(ydisc=20)
acft.addwake()
acft.edit_parameters(par='Uinf', val=Uinf)
acft.eulersolve()
plot_Cps(sld, elems=[wng])
plot_Cls(sld, wings=[wng])
sld.plotgeometry(xlim=[-b/2, b/2], ylim=[-b/2, b/2], zlim=[-b/2, b/2])
