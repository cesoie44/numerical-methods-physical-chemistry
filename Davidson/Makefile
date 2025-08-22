#
#   Makefile
#
FC = gfortran
#FFLAGS = -fopenmp -g -fbacktrace -fcheck=all
FFLAGS = -g -fbacktrace -fcheck=all  -O3
LIBS = -llapack -lblas
MODS   = utils.o ortho.o davidson_solver.o
OBJS   = main.o
#
all:    $(MODS) $(OBJS)
	$(FC) $(FFLAGS) -o main.exe $(OBJS) $(MODS) $(LIBS)
#
%.o: %.f
	$(FC) $(FFLAGS) -c $*.f
%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90
#
clean:
	rm -fr $(MODS) $(OBJS) *.exe *.mod
