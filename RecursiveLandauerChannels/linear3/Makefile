#
#===============================
# WanT package
#===============================
#
include ../make.sys


OBJS= \
cleanup.o \
coeff_m.o \
control_variable.o \
data_m.o \
dim_variable.o \
distance_m.o \
identite_m.o \
input_base.o \
input.o \
input_parameters.o \
io_module.o \
orbitale_m.o \
output.o \
recursion.o \
recursion_function_m.o \
recursion_init.o \
recursion_m.o \
startup.o \
subspace_variable.o \
summary.o

BIN= \
linear.x \

EXTERNAL_LIBS= ../iotk/libiotk.a ../libs/libwant.a 
BIN_DEPS=      $(OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all : $(BIN)

linear.x: linear.o $(BIN_DEPS)
	$(LD) -o $@ linear.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear3/$@ . )

#
# cleaning
#
clean:
	- /bin/rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

