#
#===============================
# package
#===============================
#
include ./make.sys


OBJS= \
IOModule.o \
Calculate.o \
kind.o \
constants.o 

#EXTERNAL_LIBS=  ./libwant.a
BIN_DEPS=      $(OBJS)  $(EXTERNAL_LIBS)

#
# main targets
#
all: ScriptPPV.x ScriptDensity.x ScriptPlotDensity.x ScriptPlotPotential.x ScriptPotential.x ScriptDipole.x ScriptWavefunctionInPotential.x Plot.x

ScriptPPV.x: ScriptPPV.o $(BIN_DEPS)
	$(LD) -o $@ ScriptPPV.o $(OBJS)  $(LDFLAGS)

ScriptDensity.x: ScriptDensity.o $(BIN_DEPS)
	$(LD) -o $@ ScriptDensity.o $(OBJS)  $(LDFLAGS)

ScriptPlotDensity.x: ScriptPlotDensity.o $(BIN_DEPS)
	$(LD) -o $@ ScriptPlotDensity.o $(OBJS)  $(LDFLAGS)

ScriptPlotPotential.x: ScriptPlotPotential.o $(BIN_DEPS)
	$(LD) -o $@ ScriptPlotPotential.o $(OBJS)  $(LDFLAGS)

ScriptPotential.x: ScriptPotential.o $(BIN_DEPS)
	$(LD) -o $@ ScriptPotential.o $(OBJS)  $(LDFLAGS)

ScriptDipole.x: ScriptDipole.o $(BIN_DEPS)
	$(LD) -o $@ ScriptDipole.o $(OBJS)  $(LDFLAGS)

ScriptWavefunctionInPotential.x: ScriptWavefunctionInPotential.o $(BIN_DEPS)
	$(LD) -o $@ ScriptWavefunctionInPotential.o $(OBJS)  $(LDFLAGS)

Plot.x: Plot.o $(BIN_DEPS)
	$(LD) -o $@ Plot.o $(OBJS)  $(LDFLAGS)
	 
clean:
	rm -rf *.x

# dependencies
include make.depend

