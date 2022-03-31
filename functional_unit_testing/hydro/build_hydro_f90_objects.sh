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

sed "/implicit none/i $new_fates_r8_str" ../../main/FatesConstantsMod.F90 > f90_src/FatesConstantsMod.F90
sed -i "/implicit none/i $new_fates_int_str" f90_src/FatesConstantsMod.F90
sed -i "/private /i public :: fates_r8" f90_src/FatesConstantsMod.F90
sed -i "/private /i public :: fates_int" f90_src/FatesConstantsMod.F90

# Delete the old lines

sed -i "/$old_fates_r8_str/d" f90_src/FatesConstantsMod.F90
sed -i "/$old_fates_int_str/d" f90_src/FatesConstantsMod.F90

# Build the new file with constants

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesConstantsMod.o  f90_src/FatesConstantsMod.F90

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/UnitWrapMod.o f90_src/UnitWrapMod.F90

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesHydroWTFMod.o ../../biogeophys/FatesHydroWTFMod.F90

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/HydroUnitWrapMod.o f90_src/HydroUnitWrapMod.F90




