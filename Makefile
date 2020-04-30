#LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_gf_lp64 -lmkl_sequential -lmkl_lapack95_lp64
FFLAGS=-O3
FC=gfortran

OBJ =     \
modules   \
main      \
step      \

OBJs=$(patsubst %, .%.o, $(OBJ))

all: $(OBJs)
	$(FC) $(ADDFL) $(LIBS) $(FFLAGS) $(OBJs) -o swlin

.%.o: %.f90
	$(FC) $(ADDFL) $(FFLAGS) -c -o $@ $<
.%.o: %.f
	$(FC) $(ADDFL) $(FFLAGS) -c -o $@ $<

.main.o: .modules.o
.step.o: .modules.o

clean:
	rm .*.o *.mod
