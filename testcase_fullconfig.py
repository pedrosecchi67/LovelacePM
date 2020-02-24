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

'''script to test full configuration calculations based on high-end functions'''

#fuselage only test
sld=Solid()
pts=lin_fus_surf()
sld.addpatch(pts, wraps=[0])
sld.end_preprocess()
sld.plotgeometry()

'''wingafl='n4412'
c=0.3
b=2.0
fuselage_width=0.1
AR=b/c
wake_offset=1000.0
xdisc=20
ydisc=30

alpha=radians(5.0)
beta=0.0
Uinf=0.01

wakecombs=[]
sld=Solid([])
lwing=wing(sld)
lw_horzlines, lw_vertlines, lw_paninds, lw_sldpts, lw_wakecombs=lwing.patchcompose(afls=[wingafl], chords=[c], ypos=[-b/1, -fuselage_width/2], \
    zbase=0.0, xpos=[0.0], gammas=[0.0], yspacing=np.linspace(-b/2, -fuselage_width/2, ydisc), xdisc=xdisc)
wakecombs+=lw_wakecombs
rwing=wing(sld)
rw_horzlines, rw_vertlines, rw_paninds, rw_sldpts, rw_wakecombs=rwing.patchcompose(afls=[wingafl], chords=[c], ypos=[-fuselage_width/2, b/1], \
    zbase=0.0, xpos=[0.0], gammas=[0.0], prevlines={'left':[lw_vertlines[i][0] for i in range(len(lw_vertlines))]}, \
        yspacing=np.linspace(-fuselage_width/2, b/2, ydisc), xdisc=xdisc)
wakecombs+=rw_wakecombs

sld.end_preprocess()
sld.genwakepanels(wakecombs=wakecombs, offset=wake_offset, a=alpha, b=beta)
#rwing.plotpatch(xlim=[-b/2, b/2], ylim=[-b/2, b/2], zlim=[-b/2, b/2])
sld.eulersolve(a=alpha, b=beta, Uinf=Uinf)
sld.plotgeometry(xlim=[-b/2, b/2], ylim=[-b/2, b/2], zlim=[-b/2, b/2])
lwing.plotpressure3D()
rwing.plotpressure3D()'''