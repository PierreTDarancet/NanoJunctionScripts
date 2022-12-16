#
#===============================
# package
#===============================
#
include ../make.sys


OBJS= \
alphagrid.o \
biasgrid.o  \
control_m.o \
egrid.o \
green_functions.o \
green_workspace_m.o \
hamiltonian_construction_m.o \
hamiltonian_m.o \
input.o \
occupation.o \
print_m.o \
photongrid.o \
photonvariables.o \
transport_formula.o 

EXTERNAL_LIBS=  ../libs/libwant.a
BIN_DEPS=      $(OBJS)  $(EXTERNAL_LIBS)

#
# main targets
#
all: conductor_model_beta.x 


conductor_model_beta.x: conductor.o $(BIN_DEPS)
	$(LD) -o $@ conductor.o $(OBJS)  $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../transport_model_molecular_heterojunction/$@ . )

clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

