#
#===============================
# package
#===============================
#
include ./make.sys


OBJS= \
control_m.o \
egrid.o \
green_functions.o \
green_workspace_m.o \
hamiltonian_m.o \
input.o \
kgrid_module.o \
occupation.o \
read_module.o \
se_functions.o \
transport_formula.o 

#EXTERNAL_LIBS=  ../libs/libwant.a
BIN_DEPS=      $(OBJS)  $(EXTERNAL_LIBS)

#
# main targets
#
all: conductor_MW_photon.x 


conductor_MW_photon.x : conductor.o $(BIN_DEPS)
	$(LD) -o $@ conductor.o $(OBJS)  $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../transport_model_molecular_heterojunction/$@ . )

clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

