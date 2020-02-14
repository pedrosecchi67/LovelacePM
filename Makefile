MODULE=toolkit
SRC=toolkit.f90
LIBFLAGS=-llapack -lblas

CC=f2py
CFLAGS=-c $(SRC) $(LIBFLAGS) -m $(MODULE)

all: compile

.PHONY: compile
compile:
	$(CC) $(CFLAGS)