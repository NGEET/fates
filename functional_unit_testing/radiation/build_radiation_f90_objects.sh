#!/bin/bash

# Path to FATES src

FC='gfortran'

F_OPTS="-shared -fPIC -g -O0 -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wall"
#F_OPTS="-shared -fPIC -O"


MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.mod

# Build the new file with constants

${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/WrapShrMod.o f90_src/WrapShrMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/FatesRadiationMemMod.o ../../radiation/FatesRadiationMemMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/TwoStreamMLPEMod.o  ../../radiation/TwoStreamMLPEMod.F90
${FC} ${F_OPTS} -I bld/ ${MOD_FLAG} bld/ -o bld/RadiationWrapMod.o f90_src/RadiationWrapMod.F90




