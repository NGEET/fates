PROG = test_norman

SRCS = FatesConstantsMod.F90 FatesGlobals.F90 log_mod.f90                      \
			 FatesHydraulicsMemMod.F90 FatesHydroWTFMod.F90 FatesLitterMod.F90       \
			 CohortModule.f90 FatesRunningMeanMod.F90 PatchModule.f90                \
			 PRTGenericMod.F90 PRTParametersMod.F90 norman_constants.f90             \
			 IOMod.f90 norman_radiation.f90 test_norman.f90 OrbitalValsMod.f90
			 

OBJS = FatesConstantsMod.o FatesGlobals.o log_mod.o                            \
			 FatesHydraulicsMemMod.o FatesHydroWTFMod.o FatesLitterMod.o             \
			 CohortModule.o FatesRunningMeanMod.o PatchModule.o                      \
			 PRTGenericMod.o PRTParametersMod.o norman_constants.o                   \
			 IOMod.o norman_radiation.o test_norman.o OrbitalValsMod.o

## need both
LIBS = -L/usr/local/Cellar/netcdf/4.9.2/lib -lnetcdf -L/usr/local/Cellar/netcdf-fortran/4.6.0/lib -lnetcdff -lnetcdf 

INCLUDES = -I/usr/local/Cellar/netcdf-fortran/4.6.0/include

F90 =  gfortran

##DBG = -g -fbacktrace -ffpe-trap=invalid,zero,overflow,underflow,denormal

##OPT = -O2 -axcore-avx2 -xavx -align array64byte -ipo

F90FLAGS = $(DBG) $(OPT)

all: $(PROG)

## link objects into a library (don't need includes?)
$(PROG): $(OBJS)
	$(F90) -o $@ $(OBJS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90 .F90 .f95 .f

# percent can stand for anything XX.o depends on XX.mod
%.o: %.mod

# for each .o file, make it from a 
.f90.o .f95.o .F90.o:
	$(F90) $(F90FLAGS) $(INCLUDES) -c $<

# dependencies (:) - if norman_constants.o gets updated, remake parameter_mod, etc.
FatesGlobals.o: FatesConstantsMod.o                         
FatesHydraulicsMemMod.o: FatesConstantsMod.o FatesHydroWTFMod.o norman_constants.o
FatesHydroWTFMod.o: FatesConstantsMod.o FatesGlobals.o log_mod.o
FatesLitterMod.o: FatesConstantsMod.o FatesGlobals.o log_mod.o       
CohortModule.o: FatesConstantsMod.o PRTGenericMod.o FatesHydraulicsMemMod.o
FatesRunningMeanMod.o: FatesConstantsMod.o log_mod.o FatesGlobals.o norman_constants.o
PatchModule.o: CohortModule.o FatesLitterMod.o FatesRunningMeanMod.o FatesConstantsMod.o PRTGenericMod.o norman_constants.o                    
PRTGenericMod.o: FatesConstantsMod.o FatesGlobals.o log_mod.o PRTParametersMod.o
PRTParametersMod.o: FatesConstantsMod.o              
IOMod.o: norman_constants.o
OrbitalValsMod.o: IOMod.o
norman_radiation.o: norman_constants.o
norman_constants.o: FatesConstantsMod.o
test_norman.o: norman_radiation.o norman_constants.o IOmod.o PatchModule.o OrbitalValsMod.o


