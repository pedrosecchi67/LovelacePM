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

sld=Solid()

fig=plt.figure()
ax=plt.axes(projection='3d')

fuselage=standard_body(sld, body_thdisc=40, nose_loc=np.array([-0.4, 0.0, 0.0]))
fuselage.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

sect1=wing_section(CA_position=np.array([0.0, -0.6, 0.0]), c=0.1)
sect2=wing_section(CA_position=np.array([0.0, -0.05, 0.0]), c=0.2)

wingquad=wing_quadrant(sld, sect1=sect1, sect2=sect2)
wingquad.trim_bybody(fuselage)
wingquad.plot_input(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], \
    ax=ax, fig=fig)

plt.show()

wingquad.patchcompose()
sld.end_preprocess()
sld.plotnormals(xlim=[-0.5, 1.5], ylim=[-1.0, 1.0], zlim=[-1.0, 1.0], factor=0.05)