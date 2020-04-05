MODULE1=toolkit
SRC1=toolkit.f90
MODULE2=fdyn
SRC2=fdyn.f90
LIBFLAGS=-llapack -lblas

CC=f2py
CFLAGS=$(LIBFLAGS)

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
	$(CC) -c $(SRC1) -m $(MODULE1) $(CFLAGS)
	$(CC) -c $(SRC2) -m $(MODULE2) $(CFLAGS)