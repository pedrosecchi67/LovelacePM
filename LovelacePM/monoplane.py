import numpy as np
import os
from math import *
from LovelacePM import *
import matplotlib.pyplot as plt

ordir=os.getcwd()
os.chdir(os.path.dirname(os.path.abspath(__file__)))

Uinf=200/3.6 #200kph

fuselage_height=2.2
cabin_screen_height=0.8
motor_spinner_height=0.9
motor_width=1.2
radiator_height=0.7
fuselage_width=1.6
motor_length=2.5
cabin_aft_projection=0.5
screen_aft_projection=1.0
cabin_length=2.4
tailcone_length=5.0
tailcone_height=1.1
spinner_length=0.3

wing_height=0.3
b=13.0; c=2.0; taper=0.5; dihedron=5.0

emph_b=4.0
emph_c=1.0
emph_taper=0.6
emph_height=1.1
emph_offset=0.6

empv_b=2.0
empv_c=1.0
empv_taper=0.6
empv_height=1.1
empv_offset=0.6

sld=Solid()

fusdefsect=smooth_angle_defsect_function(r_1x=0.5, r_2x=0.5, r_2y=0.5, r_1y=0.5, ldisc=20, thdisc=10)

sects=[fusdefsect(center=np.array([-spinner_length-motor_length, 0.0, motor_spinner_height]), R=0.0), \
    fusdefsect(center=np.array([-motor_length, 0.0, motor_spinner_height]), R=motor_width/2, z_expand=radiator_height/motor_width), \
        fusdefsect(center=np.array([-screen_aft_projection, 0.0, (fuselage_height-cabin_screen_height)/2]), R=(fuselage_height-cabin_screen_height)/2, \
            y_expand=motor_width/(fuselage_height-cabin_screen_height)), \
            fusdefsect(center=np.array([-cabin_aft_projection, 0.0, fuselage_height/2]), R=fuselage_width/2, z_expand=fuselage_height/fuselage_width), \
                fusdefsect(center=np.array([cabin_length-cabin_aft_projection, 0.0, fuselage_height/2]), R=fuselage_width/2, z_expand=fuselage_height/fuselage_width), \
                    fusdefsect(center=np.array([cabin_length-cabin_aft_projection+tailcone_length, 0.0, tailcone_height]), R=0.0)]

fuselage=body(sld, sections=sects)

root_1=wing_section(afl='atr_foil', afldir=os.getcwd(), c=c, CA_position=np.array([0.0, -fuselage_width/2, wing_height]), xdisc=20)
tip_1=wing_section(afl='atr_foil', afldir=os.getcwd(), c=c*taper, CA_position=np.array([0.0, -b/2, wing_height+tan(radians(dihedron))*b/2]), xdisc=20, closed=True)
wing1quad=wing_quadrant(sld, sect1=tip_1, sect2=root_1)
wing1=wing(sld, wingquads=[wing1quad])
root_2=wing_section(afl='atr_foil', afldir=os.getcwd(), c=c, CA_position=np.array([0.0, fuselage_width/2, wing_height]), xdisc=20)
tip_2=wing_section(afl='atr_foil', afldir=os.getcwd(), c=c*taper, CA_position=np.array([0.0, b/2, wing_height+tan(radians(dihedron))*b/2]), xdisc=20, closed=True)
wing2quad=wing_quadrant(sld, sect1=root_2, sect2=tip_2)
wing2=wing(sld, wingquads=[wing2quad])

wing1.trim_bybody(fuselage, sectside=2)
wing2.trim_bybody(fuselage, sectside=1)

horzemp_tip1=wing_section(afl='atr_emp', afldir=os.getcwd(), CA_position=np.array([cabin_length-cabin_aft_projection+tailcone_length-0.75*emph_c-emph_offset, -emph_b/2, emph_height]), \
    c=emph_c*emph_taper, inverse=True, xdisc=20, closed=True)
horzemp_root1=wing_section(afl='atr_emp', afldir=os.getcwd(), CA_position=np.array([cabin_length-cabin_aft_projection+tailcone_length-0.75*emph_c-emph_offset, 0.0, emph_height]), \
    c=emph_c, inverse=True, xdisc=20)
horzemp1_quad=wing_quadrant(sld, sect1=horzemp_tip1, sect2=horzemp_root1)
horzemp1=wing(sld, wingquads=[horzemp1_quad])

horzemp1.trim_bybody(fuselage, sectside=2)

horzemp_tip2=wing_section(afl='atr_emp', afldir=os.getcwd(), CA_position=np.array([cabin_length-cabin_aft_projection+tailcone_length-0.75*emph_c-emph_offset, emph_b/2, emph_height]), \
    c=emph_c*emph_taper, inverse=True, xdisc=20, closed=True)
horzemp_root2=wing_section(afl='atr_emp', afldir=os.getcwd(), CA_position=np.array([cabin_length-cabin_aft_projection+tailcone_length-0.75*emph_c-emph_offset, 0.0, emph_height]), \
    c=emph_c, inverse=True, xdisc=20)
horzemp2_quad=wing_quadrant(sld, sect2=horzemp_tip2, sect1=horzemp_root2)
horzemp2=wing(sld, wingquads=[horzemp2_quad])

horzemp2.trim_bybody(fuselage, sectside=1)

vertemp_root=wing_section(c=empv_c, afl='atr_rud', afldir=os.getcwd(), CA_position=np.array([cabin_length-cabin_aft_projection+tailcone_length-0.75*empv_c-empv_offset, 0.0, empv_height]), \
    gamma=-90.0, xdisc=15)
vertemp_tip=wing_section(c=empv_c*empv_taper, afl='atr_rud', afldir=os.getcwd(), CA_position=np.array([cabin_length-cabin_aft_projection+tailcone_length-0.75*empv_c-empv_offset, 0.0, empv_height+empv_b]), \
    gamma=-90.0, xdisc=15, closed=True)
vertemp_quad=wing_quadrant(sld, sect1=vertemp_tip, sect2=vertemp_root)
vertemp=wing(sld, wingquads=[vertemp_quad])

vertemp.trim_bybody(fuselage, sectside=2)

S, mac, _=wing1.calc_reference()
acft=aircraft(sld, elems=[fuselage, wing1, wing2, horzemp1, horzemp2, vertemp], cref=mac, Sref=2*S, bref=b, CG=np.array([0.0, 0.0, fuselage_height/2]))

acft.plot_input()

wing1.patchcompose(ydisc=20)
wing2.patchcompose(ydisc=20)
horzemp1.patchcompose(ydisc=10)
horzemp2.patchcompose(ydisc=10)
vertemp.patchcompose(ydisc=10)
fuselage.patchcompose(xdisc=100, thdisc_upright=5, thdisc_downright=5, thdisc_upleft=5, thdisc_downleft=5, leftqueue=[wing1, horzemp1], rightqueue=[wing2, horzemp2], upqueue=[vertemp])

closex=cabin_length-cabin_aft_projection+tailcone_length-0.75*emph_c-emph_offset

#inertial data from https://www.researchgate.net/publication/317371325_Design_Methodology_and_Flight_Test_Protocols_for_a_Dynamically-Scaled_General_Aviation_Aircraft/figures?lo=1&utm_source=google&utm_medium=organic

#conversion to metric units:
m=2260*0.453592
Ixx=14.59390*948*0.3048**2
Iyy=14.59390*1.346*0.3048**2
Izz=14.59390*1.967*0.3048**2
acft.edit_parameters({'Uinf':Uinf})

acft.addwake()
acft.plotgeometry()#(xlim=[closex-0.5, closex+0.5], zlim=[0.6, 1.6], ylim=[-0.5+tailcone_height, 0.5+tailcone_height])
#sld.plotnormals(xlim=[closex-0.5, closex+0.5], zlim=[0.6, 1.6], ylim=[-0.5, 0.5], factor=0.05)
acft.eulersolve()
acft.forces_report()
acft.stabreport()
acft.balance(SM=0.1)
acft.addmass(m=m, Ixx=Ixx, Iyy=Iyy, Izz=Izz)
external_history, alpha_history, beta_history, euler_history, time_history=acft.dynamic_simulation(nstep=2.5e3, dt=1e-4, perturbations={'w':10.0})
plt.plot(time_history, alpha_history)
plt.xlabel('t [s]')
plt.ylabel('Var. in angle of attack [rad]')
plt.show()
external_history, alpha_history, beta_history, euler_history, time_history=acft.dynamic_simulation(nstep=2.5e3, dt=1e-4, perturbations={'v':10.0})
plt.plot(time_history, beta_history)
plt.xlabel('t [s]')
plt.ylabel('Var. in sideslip angle [rad]')
plt.show()

os.chdir(ordir)