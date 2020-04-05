# LovelacePM

LovelacePM is an open source, Python-interpretable 3D vortex panel method code for optimized full aircraft configuration analysis.

It provides optimized access to three-dimentional, viscous-corrected potential flow calculations in the easy-to-access fashion preferred for optimization purposes in academia.

Named after Ada Lovelace and her fascination for flying machines!

## Quick start

### Defining a wing

To install LovelacePM, use:

```pip3 install LovelacePM```

For further installation instructions, check out the final section of this README.

To quickly analysie an ONERA M6 wing, start by importing the package and defining its dimensions. You can then run a simulation as the one below. File 'onerad.dat' should be located in the simulation's directory.

```
from LovelacePM import *
import numpy as np
from math import tan, radians

b=1.1963; croot=0.806; taper=0.56; sweep=26.7
alpha=5.0; Uinf=10.0

root_sect=wing_section(afl='onerad', c=croot, xdisc=30, sweep=sweep) 
#CA_position set to origin as default
left_tip_sect=wing_section(afl='onerad', c=croot*taper, \
CA_position=np.array([b*tan(radians(sweep)), -b, 0.0]), closed=True, xdisc=30, sweep=sweep)
right_tip_sect=wing_section(afl='onerad', c=croot*taper, \
CA_position=np.array([b*tan(radians(sweep)), b, 0.0]), closed=True, xdisc=30, sweep=sweep)

sld=Solid()

left_wingquad=wing_quadrant(sld, sect1=left_tip_sect, sect2=root_sect)
right_wingquad=wing_quadrant(sld, sect1=root_sect, sect2=right_tip_sect)
wng=wing(sld, wingquads=[left_wingquad, right_wingquad])

acft=aircraft(sld, elems=[wng])
acft.edit_parameters({'a':alpha, 'Uinf':Uinf, 'M':0.5})

wng.patchcompose(ydisc=50)
acft.plotgeometry()

acft.addwake()
acft.eulersolve()
acft.forces_report()
plot_Cls(sld, wings=[wng])
```

You can also run the case above without any console output using keyword argument

```echo=False```

in every "aircraft" class method execution. Check out 

```
from LovelacePM import *; from LovelacePM.documentation import *
print(aircraft.__doc__)
```

The code above follows the steps:

* Instantiating a Solid (class containing info about geometry and aerodynamic singularities);
* Instantiating wing sections (with airfoil x-discretization, in number of panels, specified by kwarg "xdisc");
* Gathering pairs of wing sections into wing quadrants (thus specifying an order);
* Gathering wing quadrants into wings (from leftmost/highest quadrant to rightmost/lowest quadrant);
* Gathering all wings and non-lifting bodies into a wing
* Using "patchcompose" methods to generate panel networks;
* Adding a wake;
* Generating an euler solution;
* Calculate forces and stability derivatives;
* Plot results, if desired.

### Adding a non-lifting body

You can add a "standard body" (with ellipsoidal nose and conical tailcone) with code as in the following example:

```
fuselage=standard_body(sld, defsect=circdefsect, body_thdisc=50, \
body_width=fuselage_width, nose_loc=np.array([-fuselage_aftlen, 0.0, -0.01]))

# ... create wings wing_left, wing_right as you want them ...

#trim right wing's leftmost/highest side so that none of its points 
#penetrate the fuselage
wing_right.trim_bybody(fuselage, sectside=1)
#trim right wing's rightmost/lowest side so that none of its points 
#penetrate the fuselage
wing_left.trim_bybody(fuselage, sectside=2)

acft=aircraft(sld, elems=[wing_left, wing_right, fuselage], Sref=Sref, bref=b)

wing_right.patchcompose(ydisc=20, ystrategy=lambda x: x)
wing_left.patchcompose(ydisc=20, ystrategy=lambda x: x)

fuselage.patchcompose(leftqueue=[wing_left], rightqueue=[wing_right], xdisc=60, \
    thdisc_upleft=5, thdisc_downleft=5, thdisc_downright=5, thdisc_upright=5)
```

Note the arguments "leftqueue", "rightqueue", etc. in the "fuselage.patchcompose" method. They are lists identifying the lifting surfaces which are connected to the body a given side of it, ordered in x-axis\'s positive direction.

Arguments "thdisc_downleft", "thdisc_upright", etc. identify the angular discretization (in number of panels) specified for the body between queues "lowqueue" and "leftqueue", "upqueue" and "rightqueue", etc.

### Customizing a body\'s geometry

You can also define a body instantiating "body" class and provide it with sections instantiated x-position by x-position, as in:

```
fusdefsect=smooth_angle_defsect_function(r_1x=0.5, r_2x=0.5, r_2y=0.5, r_1y=0.5, ldisc=20, thdisc=10)

sects=[
    fusdefsect(center=np.array([-spinner_length-motor_length, 0.0, motor_spinner_height]), R=0.0), \
    \
    fusdefsect(center=np.array([-motor_length, 0.0, motor_spinner_height]), \
    R=motor_width/2, z_expand=radiator_height/motor_width), \
    \
    fusdefsect(center=np.array([-screen_aft_projection, 0.0, (fuselage_height-cabin_screen_height)/2]), \
    R=(fuselage_height-cabin_screen_height)/2, \
    y_expand=motor_width/(fuselage_height-cabin_screen_height)), \
    \
    fusdefsect(center=np.array([-cabin_aft_projection, 0.0, fuselage_height/2]), \
    R=fuselage_width/2, z_expand=fuselage_height/fuselage_width), \
    \
    fusdefsect(center=np.array([cabin_length-cabin_aft_projection, 0.0, fuselage_height/2]), \
    R=fuselage_width/2, z_expand=fuselage_height/fuselage_width), \
    \
    fusdefsect(center=np.array([cabin_length-cabin_aft_projection+tailcone_length, 0.0, tailcone_height]), R=0.0)
]

fuselage=body(sld, sections=sects)
```

To check out the fuselage composed above, use

```
from LovelacePM import *; from LovelacePM.monoplane import *
```

Note that a function with suffix "defsect" is used in the code above. A defsect is a lambda function that defines a body section based on data simpler than its point-by-point geometry - speciffically, its maximum dimension, center position in 3D space and y and z axis scale factors (to make it easier to transform a tubular into an egg-shape fuselage, if you so desire =D). If you want further info on what is a defsect and how to customize it, check out:

```
from LovelacePM import *
from LovelacePM.documentation import *
print(body.__doc__)
```

You can also use pre-programmed defsects "circdefsect" and "squaredefsect", which can be instantiated individually as in the example above, or provided as an argument to function "standard_body", as in:

```
#circular section fuselage
fuselage=standard_body(sld, defsect=circdefsect, body_thdisc=50, \
body_width=fuselage_width, nose_loc=np.array([-fuselage_aftlen, 0.0, -0.01]))

#square section fuselage
fuselage=standard_body(sld, defsect=squaredefsect, body_thdisc=50, \
body_width=fuselage_width, nose_loc=np.array([-fuselage_aftlen, 0.0, -0.01]))
```

For further info, check:

```
from LovelacePM import *
from LovelacePM.documentation import *

help(body)
help(body_section)
help(circdefsect)
help(squaredefsect)
help(smooth_angle_defsect_function)
```
### Adding viscous corrections

LovelacePM comes with an automation of Mark Drela's viscous-inviscid Xfoil panel method code to ease strip theory viscous corrections and thus make the package complete in itself in its capability of assisting aircraft design.

To get instructions on how to use this automation for viscous corrections, check out:

```
from LovelacePM import *
from LovelacePM.documentation import *

help(xfoil_visc)
help(polar_correction)
```

Or begin by following this example:

```
#generating xfoil corrected polars
n4412_polar=polar_correction(name='n4412')
#notice that n4412.dat Selig format airfoil file must be included in your script's directory

#adding them to a wing section
sect1=wing_section(CA_position=np.array([0.0, -b/2, 0.0]), c=croot*taper, xdisc=20, correction=n4412_polar, Re=Uinf*rho*croot*taper/mu, closed=True)
```

Note that, for the automation to work, **Xfoil must be located within the user's PATH environment variable.**

## Contributing to LovelacePM

You can contribute to LovelacePM by assisting us in meeting the system requirements and code conventions listed in the sections below.

### Introduction to the project

This README refers to capabilities and usage conditions referrent to version 0.1.0. Version 0.1.0 is a beta testing version and has not yet met all system requirements listed below.

This project is subject to GNU GPL v3.0.

The system requirements defined for the final product are labelled below, along with their objective, numbered as:
(1) - performance optimization in highly iterative optimization problems;
(2) - flexibility of configuration and calculation method definition by the user;
(3) - higher accuracy and completeness of the analysis performed.

* Performing full configuration analysis of subsonic aircraft in an amount of time that doesn't exceed 50 s per AOA for a 5000-panel mesh; (1) [completed]
* Performing viscous effect estimations based on viscid-inviscid coupling exclusively with geometry, Mach and Reynolds specifications, accepting custom initial guesses for boundary layer thickness; (2, 3) [not yet implemented]
* Modelling fuselages, nacelles and engines accurately in terms of Euler solution and viscous flow; (3) [completed for non-lifting bodies]
* Approximating rotor behaviour through means of geometry and flow information, provided either by previous BEMT or LLT calculations (I. E., recieving local inflow characteristics along the rotor); (3) [not yet implemented]
* Providing user control over every step of calculations through Python functions so as to make the applied method customized for the user's own research purpose, recurring to precompiled FORTRAN backend for optimal performance only for strictly mathematical steps. (2)

### Test cases

This code has been tested according to the test cases reported in Ashok Srivastava's paper "Quadrilateral Vortex Panel Code Numerical Testing" (National Aerospace Laboratory Project Document CF 9123).

They consist in the scrips the user can summon with:

```
from LovelacePM import *

#test case for AR 20 NACA 0012 wing at alpha=10 deg
from LovelacePM.testcase_n0012_AR20 import *

#test case for AR 5 NACA 0012 wing rotating at 3 deg/s and advancing at 1 ft/s, at AOA 0
from LovelacePM.testcase_n0012_rotating import *

#test case for sphere
from LovelacePM.testcase_sphere import *

#test case for cilynder
from LovelacePM.testcase_cilynder import *
```

And compare with references in:

1. For Euler flow over a sphere with stream flowing from pole to pole equal to 3/2 times its azimuthal polar coordinate;
2. For Euler flow over a cilynder: tangential velocity equal to 2 times the sine of its polar coordinate;
3. For Euler flow over an $$AR=20$$ NACA-0012 wing: a pressure distribution as in Maskew, B.: "Prediction of Subsonic Aerodynamic Characteristics: A Case for Low Order Panel Methods". Journal of Aircraft, Feb. 1982, pp 157-153;
4. For Euler flow on a rotating NACA 0012 wing with free wake model: NASA Technical Memorandum 101024: Ashby, D., Dudley, M.: "Development and Validation of an Advanced Low Order Panel Method". Ames Research Center, Moffett Field, California, Oct. 1988.

### Code conventions

* Python's suggested naming conventions are adopted. Initials may also be referred to in capital letters;
* All functions related to Euler solutions should be either defined in paneller.py or called from it, so as to ease access. utils.py may also contain independent geometry-related functions;
* All FORTRAN backend subroutines should be defined in toolkit.f90 (for aerodynamics related functions) or fdyn.f90 (for flight dynamics related functions), obeying the FORTRAN 90 standard.

### Dependencies

The program's viscid correction module, ```xfoil_visc.py```, currently depends only on the presence of an xfoil binary either on the user's PATH environment variable or on the LovelacePM folder.

### Further instalation instructions

To clone and install directly from git, use:

```
git clone https://github.com/pedrosecchi67/LovelacePM.git
cd LovelacePM
pip install -e . #installing using pip and setup.py, thus adding the package to PYTHONPATH variable and making the repository folder the source folder for calls in any python shell
```

For custom compilation of mathematical modules, use:

```
cd LovelacePM
make
```

Yep, that's it =)

If using custom-compiled BLAS and LAPACK libraries, run:

```make LIBFLAGS="-L$(LIBDIR) -llapack -lblas"```

With ```LIBDIR``` indicating your custom directory.

If compiling for Windows, use:
```f2py -c toolkit.f90 -llapack -lblas -m toolkit```
And you should be good to go.

Pedro de Almeida Secchi (https://github.com/pedrosecchi67), 14/02/2020
