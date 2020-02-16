# LovelacePM
## Open source, Python-interpretable 3D vortex panel method code for optimized full aircraft configuration analysis

### Introduction

This program is meant to provide optimized access to three-dimentional, viscous-corrected potential flow calculations in the easy-to-access fashion preferred for optimization purposes in academia.

The project is subject to GNU GPL v3.0.

Named after Ada Lovelace and her fascination for flying machines.

The system requirements defined for the final product are labelled below, along with their purpose, numbered as:
(1) - performance optimization in highly iterative optimization problems;
(2) - flexibility of configuration and calculation method definition by the user;
(3) - higher accuracy and completeness of the analysis performed, in spite of the low-fidelity Euler equation solver that serves as its kernel.

* Performing full configuration analysis of subsonic aircraft in an amount of time that doesn't exceed 50 s per AOA for a 5000-panel mesh; (1) [not yet completed]
* Performing viscous effect estimations based on viscid-inviscid coupling exclusively with geometry, Mach and Reynolds specifications, accepting custom initial guesses for boundary layer thickness; (2, 3) [not yet implemented]
* Modelling fuselages, nacelles and engines accurately in terms of euler solution and viscous flow; (3) [not yet completed]
* Approximating rotor behaviour through means of geometry and flow information, provided either by previous BEMT or LLT calculations (I. E., recieving local inflow characteristics along the rotor); (3) [not yet implemented]
* Providing user control over every step of calculations through Python functions so as to make the applied method customized for the user's own research purpose, recurring to precompiled FORTRAN backend for optimal performance only for strictly mathematical steps. (2)

### Test cases

This product will be tested with all subsonic test cases displayed in the PAN AIR Case Manual v1.0. Results will be available in annex pdf document and LaTeX source as the program is completed and tests are performed.

### Code conventions

* Python's suggested naming conventions are adopted. Initials may also be referred to in capital letters;
* All functions related to aerodynamics should be either defined in paneller.py or called from it, so as to ease access. utils.py may also contain independent geometry-related functions;
* All FORTRAN backend subroutines should be defined in toolkit.so, obeying the FORTRAN 90 standard;
* If you don't have fun with your design, you're probably doing it wrong ;)

### Dependencies

The program's dependencies currently include:
* Python 3.7;
* numpy;
* scipy;
* matplotlib;
* mpl_toolkits (for 3D plotting);
* f2py (currently using v2 for tested versions);
* gfortran (currently using gcc 9.2.1 for tested versions);
* make (currently using 4.2.1);
* The appropriate libblas.so and liblapack.so files necessary for previous dependencies.

### Compiling instructions

```make compile```
Yep, that's it =)

If using custom-compiled BLAS and LAPACK libraries, run:
```make LIBFLAGS="-L$(LIBDIR) -llapack -lblas"```
With ```LIBDIR``` indicating your custom directory.

If compiling for Windows, use:
```f2py -c toolkit.f90 -llapack -lblas -m toolkit```
And you should be good to go.

Pedro de Almeida Secchi (https://github.com/pedrosecchi67), 14/02/2020
