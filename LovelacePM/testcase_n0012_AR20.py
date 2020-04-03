from LovelacePM import *
import numpy as np

ordir=os.getcwd()
os.chdir(os.path.dirname(os.path.abspath(__file__)))

c=1.0
AR=20.0
b=c*AR

a=10.0

sld=Solid()

sl=wing_section(afl='n0012', c=c, CA_position=np.array([0.0, -b/2, 0.0]), xdisc=30, closed=True)
sr=wing_section(afl='n0012', c=c, CA_position=np.array([0.0, b/2, 0.0]), xdisc=30, closed=True)
wngqd=wing_quadrant(sld, sect1=sl, sect2=sr)
wng=wing(sld, wingquads=[wngqd])
acft=aircraft(sld, elems=[wng])

acft.edit_parameters({'a':a})

wng.patchcompose(ydisc=30)
acft.addwake()

acft.eulersolve()

plot_Cps(sld, elems=[wng])

os.chdir(ordir)