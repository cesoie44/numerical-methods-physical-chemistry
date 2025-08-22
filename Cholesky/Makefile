#
#   Makefile
#
FC = gfortran
FFLAGS = -fopenmp -O3
LIBS = -lblas -llapack

# Moduli e oggetti
MODS = utils.o cholesky_solver.o
OBJS = main.o

#
all: $(MODS) $(OBJS)
	$(FC) $(FFLAGS) -o main.exe $(OBJS) $(MODS) $(LIBS)

#
%.o: %.f
	$(FC) $(FFLAGS) -c $<

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#
clean:
	rm -f $(MODS) $(OBJS) *.exe *.mod

