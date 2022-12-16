#
#===============================
# WanT package
#===============================
#
include ../make.sys


OBJS= \
cleanup.o \
control.o \
egrid.o \
gamma.o \
green_C_inv.o \
green_C_tridiag.o \
in_hamiltonian_m.o \
in_matrix_m.o \
input.o \
input_parameters.o \
summary.o \

WANNIER_OBJS= \
../wannier/io_module.o \
../wannier/startup.o \
../wannier/version.o


EXTERNAL_LIBS= ../iotk/libiotk.a ../libs/libwant.a
BIN_DEPS=      $(OBJS) $(WANNIER_OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all: conductor_qudot_imp_gamma.x 


conductor_qudot_imp_gamma.x: conductor.o $(BIN_DEPS)
	$(LD) -o $@ conductor.o $(OBJS) $(WANNIER_OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../transport_qudot/$@ . )

clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

