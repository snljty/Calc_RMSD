# Makefile for calc_rmsd

SHELL = cmd
FC := gfortran
FLINKER = $(FC)
OPTS = -s
LAPACKROOT := C:\lapack-3.10.1
LIBPATH := -L $(LAPACKROOT)\lib
LIB := -l lapack -l blas

.PHONY: all
all: calc_rmsd.exe

%.exe: %.obj
	@echo Linking $@ against $^ ...
	$(FLINKER) -o $@ $^ $(LIBPATH) $(LIB) $(OPTS) -static

%.obj: %.f90
	@echo Compiling $@ from $< ...
	$(FC) -o $@ -c $< $(OPTS) -cpp

.PHONY: clean
clean:
	-del /q calc_rmsd.exe 2> NUL

.PHONY: test
test: all ref.xyz traj.xyz
	-del /q RMSD.txt RMSD.txt.bak 2> NUL
	.\calc_rmsd.exe ref.xyz traj.xyz

