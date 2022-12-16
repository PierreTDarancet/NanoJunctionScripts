#
#===============================
# WanT package
#===============================
#
include ../make.sys


OBJS= \
cleanup.o \
convert2_control_variable.o \
convert2_input.o \
convert2_input_parameters.o \
convert2_m.o \
convert3_control_variable.o \
convert3_input.o \
convert3_input_parameters.o \
convert3_m.o \
convert4_control_variable.o \
convert4_input.o \
convert4_input_parameters.o \
convert4_m.o \
convert_control_variable.o \
convert_input.o \
convert_input_parameters.o \
convert_m.o \
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
convert4.x \
convert3.x \
convert2.x \
convert.x \
linear.x \

EXTERNAL_LIBS= ../iotk/libiotk.a ../libs/libwant.a 
BIN_DEPS=      $(OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all : $(BIN)
convert4.x: convert4.o $(BIN_DEPS)
	$(LD) -o $@ convert4.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear2/$@ . )

convert3.x: convert3.o $(BIN_DEPS)
	$(LD) -o $@ convert3.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear2/$@ . )

convert2.x: convert2.o $(BIN_DEPS)
	$(LD) -o $@ convert2.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear2/$@ . )

convert.x: convert.o $(BIN_DEPS)
	$(LD) -o $@ convert.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear2/$@ . )

linear.x: linear.o $(BIN_DEPS)
	$(LD) -o $@ linear.o $(OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../linear2/$@ . )

#
# cleaning
#
clean:
	- /bin/rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

