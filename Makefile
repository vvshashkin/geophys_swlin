LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_gf_lp64 -lmkl_sequential -lmkl_lapack95_lp64
FFLAGS=-O3 -traceback

OBJ =     \
modules   \
main      \
step      \

OBJs=$(patsubst %, .%.o, $(OBJ))

all: $(OBJs)
	ifort $(ADDFL) $(LIBS) $(FFLAGS) $(OBJs) -o swlin

.%.o: %.f90
	ifort $(ADDFL) $(FFLAGS) -c -o $@ $<
.%.o: %.f
	ifort $(ADDFL) $(FFLAGS) -c -o $@ $<

.main.o: .modules.o
.step.o: .modules.o

clean:
	rm .*.o *.mod
