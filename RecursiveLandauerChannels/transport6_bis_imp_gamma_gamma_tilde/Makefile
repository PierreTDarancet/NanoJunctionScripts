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
green_C_tridiag.o \
in_matrix_m.o \
input.o \
input_parameters.o \
summary.o \
transmittance_min.o \

WANNIER_OBJS= \
../wannier/io_module.o \
../wannier/startup.o \
../wannier/version.o


EXTERNAL_LIBS= ../iotk/libiotk.a ../libs/libwant.a
BIN_DEPS=      $(OBJS) $(WANNIER_OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all: conductor6_imp_gamma.x 


conductor6_imp_gamma_gamma_tilde.x: conductor.o $(BIN_DEPS)
	$(LD) -o $@ conductor.o $(OBJS) $(WANNIER_OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../transport6_bis_imp_gamma_gamma_tilde/$@ . )

clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

