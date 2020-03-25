# LovelacePM
## Open source, Python-interpretable 3D vortex panel method code for optimized full aircraft configuration analysis

This README refers to capabilities and usage conditions referrent to version 0.0.3. Version 0.0.3 is an MVP and is therefore not bug free, and has not yet met all system requirements listed below.

### Introduction

This program is meant to provide optimized access to three-dimentional, viscous-corrected potential flow calculations in the easy-to-access fashion preferred for optimization purposes in academia.

The project is subject to GNU GPL v3.0.

Named after Ada Lovelace and her fascination for flying machines.

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

This product will be tested with all subsonic test cases displayed in the PAN AIR Case Manual v1.0 in its transition from MVP/alpha version (v0.X.Y) to beta version (v1.X.Y). Results will be available in annex pdf document and LaTeX source as the program is completed and tests are performed.

### Code conventions

* Python's suggested naming conventions are adopted. Initials may also be referred to in capital letters;
* All functions related to Euler solutions should be either defined in paneller.py or called from it, so as to ease access. utils.py may also contain independent geometry-related functions;
* All FORTRAN backend subroutines should be defined in toolkit.so, obeying the FORTRAN 90 standard.

### Dependencies

The program's viscid correction module, ```xfoil_visc.py```, currently depends only on the presence of an xfoil binary either on the user's PATH environment variable or on the LovelacePM folder.

### Instalation instructions

To install via pip, use:

```pip3 install LovelacePM==0.0.3```

To clone and install directly from git, use:

```git clone https://github.com/pedrosecchi67/LovelacePM.git
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
