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
gamma_m.o \
green_analytique.o \
green_tridiag.o \
green_want.o \
in_matrix_m.o \
input.o \
input_parameters.o \
matrix_m.o \
summary.o \
transfer.o \

WANNIER_OBJS= \
../wannier/io_module.o \
../wannier/startup.o \
../wannier/version.o


EXTERNAL_LIBS= ../iotk/libiotk.a ../libs/libwant.a
BIN_DEPS=      $(OBJS) $(WANNIER_OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all: gamma_c.x 


gamma_c.x: gamma_c.o $(BIN_DEPS)
	$(LD) -o $@ gamma_c.o $(OBJS) $(WANNIER_OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../gamma/$@ . )

clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

