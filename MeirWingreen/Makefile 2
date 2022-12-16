#
#===============================
# WanT package
#===============================
#
include ../make.sys


OBJS= \

WANNIER_OBJS= \


EXTERNAL_LIBS=  ../libs/libwant.a
BIN_DEPS=      $(OBJS) $(WANNIER_OBJS) $(EXTERNAL_LIBS)

#
# main targets
#
all: conductor_model_molecular_HJ.x 


conductor_model_molecular_HJ.x: conductor.o $(BIN_DEPS)
	$(LD) -o $@ conductor.o $(OBJS) $(WANNIER_OBJS) $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../transport_model_molecular_heterojunction/$@ . )

clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

