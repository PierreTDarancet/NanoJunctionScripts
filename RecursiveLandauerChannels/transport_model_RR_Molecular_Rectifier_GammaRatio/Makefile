#
#===============================
# package
#===============================
#
include ../make.sys


OBJS= \
control_m.o \
biasgrid.o \
gapgrid.o \
gammagrid.o  \
input.o 

EXTERNAL_LIBS=  ../libs/libwant.a
BIN_DEPS=      $(OBJS)  $(EXTERNAL_LIBS)

#
# main targets
#
all: conductor_Paper.x 


conductor_Paper.x: conductor.o $(BIN_DEPS)
	$(LD) -o $@ conductor.o $(OBJS)  $(LDFLAGS)
	
clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

