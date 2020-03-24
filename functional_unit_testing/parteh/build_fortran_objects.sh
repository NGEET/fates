#!/bin/bash

# Path to FATES src

FC='gfortran'

F_OPTS="-shared -fPIC -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check"
#F_OPTS="-shared -fPIC -O"


MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.mod


# First copy over the FatesConstants file, but change the types of the fates_r8 and fates_int

old_fates_r8_str=`grep -e integer ../../main/FatesConstantsMod.F90 | grep fates_r8 | sed 's/^[ \t]*//;s/[ \t]*$//'`
new_fates_r8_str='use iso_c_binding, only: fates_r8  => c_double'

old_fates_int_str=`grep -e integer ../../main/FatesConstantsMod.F90 | grep fates_int | sed 's/^[ \t]*//;s/[ \t]*$//'`
new_fates_int_str='use iso_c_binding, only: fates_int => c_int'

# Add the new lines (need position change, don't swap)

sed "/implicit none/i $new_fates_r8_str" ../../main/FatesConstantsMod.F90 > f90src/FatesConstantsMod.F90
sed -i "/implicit none/i $new_fates_int_str" f90src/FatesConstantsMod.F90
sed -i "/private /i public :: fates_r8" f90src/FatesConstantsMod.F90
sed -i "/private /i public :: fates_int" f90src/FatesConstantsMod.F90

# Delete the old lines

sed -i "/$old_fates_r8_str/d" f90src/FatesConstantsMod.F90
sed -i "/$old_fates_int_str/d" f90src/FatesConstantsMod.F90


# This re-writes the wrapper so that it uses all the correct parameters
# in FatesAllometryMod.F90
python AutoGenVarCon.py


# Procedure for auto-generating AllomUnitWrap
# 1) scan FatesAllometry and create list of EDPftVarcon_inst variables
# 2) scan EDpftVarcon and get the name of the in-file parameter names associated
#    with these variables



# Build the new file with constants

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesConstantsMod.o  f90src/FatesConstantsMod.F90


# The Parameter Definitions
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/ -o bld/PRTParametersMod.o ../../parteh/PRTParametersMod.F90



${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/UnitWrapMod.o f90src/UnitWrapMod.F90




# Generic Integration routines (all native types except defined constants)
gfortran $F_OPTS $MOD_FLAG bld/ -I bld/  -o bld/FatesIntegratorsMod.o ../../main/FatesIntegratorsMod.F90

# Allometry Module, take this from FATES directly
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/ -o bld/FatesAllometryMod.o ../../biogeochem/FatesAllometryMod.F90 

# The Generic (parent) PARTEH module
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/ -o bld/PRTGenericMod.o ../../parteh/PRTGenericMod.F90

# Loss Fluxes and phenology
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/ -o bld/PRTLossFluxesMod.o ../../parteh/PRTLossFluxesMod.F90

# The carbon-only PARTEH module
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/ -o bld/PRTAllometricCarbonMod.o ../../parteh/PRTAllometricCarbonMod.F90

# The CNP allometric target model 
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/ -o bld/PRTAllometricCNPMod.o ../../parteh/PRTAllometricCNPMod.F90

# Initialize PARTEH instance and mapping functions
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/ -o bld/FatesPARTEHWrapMod.o f90src/FatesPARTEHWrapMod.F90

# The cohort instances and initialization
${FC} ${F_OPTS} $MOD_FLAG bld/ -I bld/  -o bld/FatesCohortWrapMod.o f90src/FatesCohortWrapMod.F90



