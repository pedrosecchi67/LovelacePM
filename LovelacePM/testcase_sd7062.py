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

a=radians(0.0)
Uinf=0.05
b=1.2
taper=1.0
croot=0.3

sld=Solid()

sect1=wing_section(afl='sd7062', CA_position=np.array([0.0, -b/2, 0.0]), c=croot*taper, xdisc=50, closed=True)#, xstrategy=lambda x: x)
sect2=wing_section(afl='sd7062', CA_position=np.array([0.0, b/2, 0.0]), c=croot, xdisc=50, closed=True)#, xstrategy=lambda x: x)

wng1=wing_quadrant(sld, sect1=sect1, sect2=sect2)
wng=wing(sld, wingquads=[wng1])
acft=aircraft(sld, elems=[wng], Sref=b*croot*(1+taper)/2)
wng.patchcompose(ydisc=50)
acft.edit_parameters({'a':a, 'Uinf':Uinf})
acft.addwake()

'''sld.plotnormals(xlim=[-0.6, 0.6], ylim=[-0.6, 0.6], zlim=[-0.6, 0.6], factor=0.1)
sld.plotnormals(xlim=[-0.2, 0.2], ylim=[-0.8, -0.4], zlim=[-0.2, 0.2], factor=0.1)'''
acft.eulersolve()
acft.forces_report()
acft.stabreport()
'''plot_Cps(sld, elems=[wng])
plot_Cls(sld, wings=[wng])
plot_Cds(sld, wings=[wng])
plot_Cms(sld, wings=[wng])
plot_gammas(sld, wings=[wng])
sld.plotgeometry(xlim=[-0.6, 0.6], ylim=[-0.6, 0.6], zlim=[-0.6, 0.6])'''