BUILDDIR = ../build
BINDIR = ../bin
LIBDIR = ../lib

EXEC = ${BINDIR}/genec
EXEC_MAKEINI = ${BINDIR}/genec_makeini
LIB = ${LIBDIR}/libgenec.a

FFLAGS = -Wall -fimplicit-none -J$(BUILDDIR)
#FFLAGS = $(FFLAGS) -fallow-invalid-boz -pedantic -ffpe-trap=invalid,zero,overflow,underflow

PGPLOT = TRUE

###########
## IFORT ##
###########
#BEWARE: if you compile with ifort,
#you need to uncomment the corresponding line in the compilation rules !
#FC = ifort
#LIBS = -lm
#FFLAGS = -debug -check -warn -r8 -zero -g
#FFLAGS = -r8 -O3

###########
### PGI ###
###########
#LIBS = -lm
#FC = pgf95
#FFLAGS = -r8 -O2

##############
## GFORTRAN ##
##############
FC = gfortran

#on the Mac (MacPorts):
ifeq ($(PGPLOT),TRUE)
        LIBS = -lm -L/opt/local/lib -lpgplot
else
        LIBS = -lm
endif

# Modules for GENEC - the order is important
SRC_BASE = \
io_definitions.f90 \
evol.f90 \
storage.f90 \
const.f90 \
caramodele.f90 \
inputparam.f90 \
interpolation.f90 \
nagmod.f90 \
SmallFunc.f90 \
abundmod.f90 \
safestop.f90 \
diffadvmod.f90 \
equadiffmod.f90 \
geomod.f90 \
magmod.f90 \
rotmod.f90 \
strucmod.f90 \
ionisation.f90 \
convection.f90 \
timestep.f90 \
PrintAll.f90 \
omegamod.f90 \
EOS.f90 \
nablas.f90 \
energy.f90 \
opacity.f90 \
envelope.f90 \
LayersShift.f90 \
advection.f90 \
diffusion.f90 \
chemicals.f90 \
bintidemod.f90 \
winds.f90

SRC_PLOT = $(SRC_BASE) \
	   PGPlotModule.f90 \
	   henyey_solver.f90 \
	   WriteSaveClose.f90 \
	   genec.f90
SRC_NOPLOT = $(SRC_BASE) \
	     PGPlotModule_Null.f90 \
	     henyey_solver.f90 \
	     WriteSaveClose.f90 \
	     genec.f90
ifeq ($(PGPLOT),TRUE)
	SRC = $(SRC_PLOT) main.f90
else
	SRC = $(SRC_NOPLOT) main.f90
endif

OBJ = $(patsubst %.f90, $(BUILDDIR)/%.o, $(SRC))

SRC_MAKEINI = \
io_definitions.f90 \
evol.f90 \
storage.f90 \
const.f90 \
caramodele.f90 \
inputparam.f90 \
abundmod.f90 \
inichemmod.f90 \
interpolmod.f90 \
modinimod.f90 \
makeini.f90 \
main_inifile.f90

OBJ_MAKEINI = $(patsubst %.f90, $(BUILDDIR)/%.o, $(SRC_MAKEINI))

SRC_LIB = $(SRC_NOPLOT) \
	  helpers.f90 \
	  inichemmod.f90 \
	  interpolmod.f90 \
	  modinimod.f90 \
	  makeini.f90
OBJ_LIB = $(patsubst %.f90, $(BUILDDIR)/%.o, $(SRC_LIB))

# Targets
all: dirs genec genec_makeini

genec: dirs $(OBJ)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJ)

genec_makeini: dirs $(OBJ_MAKEINI)
	$(FC) $(FFLAGS) -o $(EXEC_MAKEINI) $(OBJ_MAKEINI)

libgenec: dirs $(OBJ_LIB)
	ar crs ${LIB} $(OBJ_LIB)

$(BUILDDIR)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

dirs:
	mkdir -p $(BUILDDIR) $(BINDIR) $(LIBDIR)

clean:
	rm -f $(BUILDDIR)/*.o $(BUILDDIR)/*.mod 

distclean: clean
	rm -f $(EXEC) $(EXEC_MAKEINI) $(LIB)
	rmdir $(BUILDDIR) $(BINDIR) $(LIBDIR)

.PHONY: clean distclean dirs
