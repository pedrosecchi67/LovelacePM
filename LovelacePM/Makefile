MODULE=toolkit
SRC=toolkit.f90
LIBFLAGS=-llapack -lblas

CC=f2py
CFLAGS=-c $(SRC) $(LIBFLAGS) -m $(MODULE)

PIP=pip3
PIPFLGS=install
PIPPCKS=numpy scipy matplotlib

all: dependencies compile

.PHONY: dependencies
dependencies:
	echo "Installing python dependencies"
	$(PIP) $(PIPFLGS) $(PIPPCKS)

.PHONY: compile
compile:
	echo "Compiling FORTRAN wrapped routines"
	$(CC) $(CFLAGS)