# On how to contribute

## Current projects

You can line up contribution intentions with current projects using our [kanbans](https://github.com/pedrosecchi67/LovelacePM/projects) on GitHub.

## File structure

You can also propose file addition (mind modularity and organization) and reorganization with pull requests.

```
LovelacePM
├── aerodynamic_output.py: functions for plotting aerodynamic calculation results
├── aircraft.py: classes containing high-end function definitions to be directly summoned 
by the user when running simulations
├── body.py: classes defining non-lifting body geometry and viscous corrections
├── control.py: classes defining control surfaces (which are controlled from within aircraft.py)
├── documentation.py: module containing documentation. Imported separately from __init__ 
for memory saving
├── fdyn.f90: legacy Fortran (f2py) routines for flight dynamics simulations
├── __init__.py: autoupdate control with LoveUpdate and importation of essential modules
├── Makefile: legacy f2py modules compilation Makefile
├── monoplane.py: monoplane example case for Tutorial 3
├── paneller.py: Solid and Panel class definitions and methods directly related 
to panel method
├── pyfdyn.py: new, pure-Python replacement for fdyn.f90
├── pytoolkit.py: new, pure-Python replacement for toolkit.f90
├── testcase_airfoil.py: test case with NACA-0012 nearly elliptic wing
├── testcase_cilynder.py: test case for Euler solution around cilynder
├── testcase_control.py: test case for control deflection on NACA-0012 straight wing
├── testcase_fullconfig.py: test case for full configuration analysis with tubular body
├── testcase_n0012_AR20.py: test case for Srivastava's AR20 NACA-0012 wing (see README)
├── testcase_n0012_rotating.py: test case for PMARC's AR5, NACA-0012 wing (see README)
├── testcase_n65415.py: test case for APAME's swept, NACA-65415 wing
├── testcase_oneram6.py: test case for ONERA M6 wing (displayed on README and tut. 1)
├── testcase_plate.py: test case for Euler solution around flat plate (obsolete)
├── testcase_random.py: test case for AIC generation (obsolete)
├── testcase_sd7062.py: test case for free wake model around straight SD-7062 wing
├── testcase_sphere.py: test case for Euler flow around sphere
├── testcase_wingonly_fullconfig.py: test case for testcase_fullconfig.py configuration
without its fuselage (obsolete)
├── toolkit.f90: legacy Fortran AIC matrix calculator (also with other lengthy mathematical 
task subroutines)
├── utils.py: small batch of short aerodynamics and airfoil manipulation functions
├── wing.py: definition for wing geometry defining classes and viscous correction application
├── xfoil_visc.py: definition of viscous correction functions and Xfoil automation
└── *.dat: airfoil files for test cases
```
