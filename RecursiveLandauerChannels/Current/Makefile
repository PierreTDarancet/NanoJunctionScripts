#
#===============================
# package
#===============================
#
include ../make.sys


OBJS= \


EXTERNAL_LIBS=  ../libs/libwant.a
BIN_DEPS=      $(OBJS)  $(EXTERNAL_LIBS)

#
# main targets
#
all: current.x currentSmall.x


current.x : current.o  $(BIN_DEPS)
	$(LD) -o $@ current.o $(OBJS)  $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../Current/$@ . )
currentSmall.x : currentSmall.o  $(BIN_DEPS)
	$(LD) -o $@ current.o $(OBJS)  $(LDFLAGS)
	- ( cd ../bin ; ln -sf ../Current/$@ . )

clean:
	rm -rf *.o *.a *core *.F90 work.pc work.pcl *.mod *.x

# dependencies
include make.depend

