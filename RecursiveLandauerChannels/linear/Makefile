#
#===============================
# WanT package
#===============================
#
include ../make.sys


OBJS= \
cleanup.o \
control_variable.o \
dim_variable.o \
input.o \
input_parameters.o \
io_module.o \
output.o \
recursion.o \
recursion_function_m.o \
recursion_init.o \
recursion_m.o \
recursion_read_matrix.o \
startup.o \
subspace_variable.o \
summary.o

BIN= \
linear.x \
util_input.x \

EXTERNAL_LIBS= ../iotk/libiotk.a ../libs/libwant.a 
BIN_DEPS=      $(OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all : $(BIN)


linear.x: linear.o $(BIN_DEPS)
	$(LD) -o $@ linear.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear/$@ . )
util_input.x:     util_input.o $(BIN_DEPS)
	$(LD) -o $@ util_input.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear/$@ . )


#
# cleaning
#
clean:
	- /bin/rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

