#!/bin/bash

# Path to FATES src

FATES_SRC=../../

CNP_SRC=/home/rgknox/SyncLRC/PARTEH/FModules/

F_OPTS="-shared -fPIC -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check"

MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.mod

gfortran $F_OPTS $MOD_FLAG bld/ -o bld/FatesConstants.o ${FATES_SRC}/main/FatesConstantsMod.F90

# Generic Integration routines (all native types except defined constants)
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/  -o bld/FatesIntegratorsMod.o ${FATES_SRC}/main/FatesIntegratorsMod.F90

# Support Modules, fairly trivial contents
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ -o bld/FatesWrapMod.o f_wrapper_modules/FatesWrapMod.F90

# This defines and fills the global pft parameter structures (stripped down from fates version)
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ -o bld/FatesPFTWrapMod.o f_wrapper_modules/FatesPFTWrapMod.F90

# Allometry Module, take this from FATES directly
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ -o bld/FatesAllometryMod.o ${FATES_SRC}/biogeochem/FatesAllometryMod.F90 

# The Generic (parent) PARTEH module
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ -o bld/PRTGenericMod.o ${FATES_SRC}/parteh/PRTGenericMod.F90

# Loss Fluxes and phenology
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ -o bld/PRTLossFluxesMod.o ${FATES_SRC}/parteh/PRTLossFluxesMod.F90

# The carbon-only PARTEH module
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ -o bld/PRTAllometricCarbonMod.o ${FATES_SRC}/parteh/PRTAllometricCarbonMod.F90

# The CNP allometric target model 
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ -o bld/PRTAllometricCNPMod.o ${CNP_SRC}/PRTAllometricCNPMod.F90

# Initialize PARTEH instance and mapping functions
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/ bld/PRTGenericMod.o bld/PRTAllometricCarbonMod.o -o bld/FatesPARTEHWrapMod.o f_wrapper_modules/FatesPARTEHWrapMod.F90

# The cohort instances and initialization
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/  -o bld/FatesCohortWrapMod.o f_wrapper_modules/FatesCohortWrapMod.F90



