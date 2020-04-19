from LovelacePM import *
import numpy as np

ordir=os.getcwd()
os.chdir(os.path.dirname(os.path.abspath(__file__)))

c=0.5
AR=20.0
b=c*AR
nu=1.72e-5
rho=1.225
Uinf=10.0
Re=Uinf*rho*c

a=0.0

n0012_visc=polar_correction('n0012', aseq=[-5.0, 5.0, 0.5])

sld=Solid()

sl=wing_section(afl='n0012', c=c, CA_position=np.array([0.0, -b/2, 0.0]), xdisc=30, closed=True, correction=n0012_visc, Re=Re)
sr=wing_section(afl='n0012', c=c, CA_position=np.array([0.0, b/2, 0.0]), xdisc=30, closed=True, correction=n0012_visc, Re=Re)
wngqd=wing_quadrant(sld, sect1=sl, sect2=sr)
wng=wing(sld, wingquads=[wngqd])
acft=aircraft(sld, elems=[wng])

acft.edit_parameters({'a':a, 'Uinf':Uinf})

wng.patchcompose(ydisc=30)
acft.addwake()

acft.eulersolve()
acft.forces_report()

plot_Cps(sld, elems=[wng])

os.chdir(ordir)