#
#===============================
# WanT package
#===============================
#
include ../make.sys


OBJS= \
cleanup.o \
coeff_m.o \
construct.o \
control_variable.o \
data_m.o \
dim_variable.o \
distance_m.o \
hamiltonian_m.o \
identite_m.o \
input_base.o \
input.o \
input_parameters.o \
io_module.o \
orbitale_m.o \
output.o \
startup.o \
summary.o

BIN= \
hamiltonian.x \

EXTERNAL_LIBS= ../iotk/libiotk.a ../libs/libwant.a 
BIN_DEPS=      $(OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all : $(BIN)


hamiltonian.x: hamiltonian.o $(BIN_DEPS)
	$(LD) -o $@ hamiltonian.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../hamiltonian/$@ . )

#
# cleaning
#
clean:
	- /bin/rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

