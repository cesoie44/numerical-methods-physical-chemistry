#
#   Makefile
#

FC = gfortran
FFLAGS = -fopenmp -O3
LIBS = -lblas -llapack

MODS = utils.o jacobi_solver.o
OBJS = main.o

#
all: $(MODS) $(OBJS)
	$(FC) $(FFLAGS) -o main.exe $(OBJS) $(MODS) $(LIBS)

#
%.o: %.f
	$(FC) $(FFLAGS) -c $*.f

%.o: %.f90
	$(FC) $(FFLAGS) -c $*.f90

#
clean:
	rm -fr $(MODS) $(OBJS) *.exe *.mod

