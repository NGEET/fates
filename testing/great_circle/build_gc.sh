#!/bin/bash

# Path to FATES src

FC='gfortran'

#F_OPTS="-fPIC -O3 -llapack"
F_OPTS="-g -fPIC"
F_OBJ_OPTS="-shared"

FATES_PATH='../../'

#F_OPTS="-fPIC -O0 -g -ffpe-trap=zero,overflow,underflow -fbacktrace -fbounds-check -Wall"

MOD_FLAG="-J"

rm -f bld/*.o
rm -f bld/*.so
rm -f bld/*.mod
rm -f bld/*.a

# Build dgesv from lapack
${FC} ${F_OPTS} -c -I bld/ -J./bld/ -o bld/libFatesConstantsMod.so ${FATES_PATH}/main/FatesConstantsMod.F90
${FC} ${F_OPTS} -c -I bld/ -J./bld/ -o bld/libWrapGreatCircleMod.so WrapGreatCircleMod.F90
${FC} ${F_OPTS} -c -I bld/ -J./bld/ -o bld/libFatesUtilsMod.so ${FATES_PATH}/main/FatesUtilsMod.F90
${FC} ${F_OPTS} -I bld/ -J./bld/ -L./bld/ -lFatesConstantsMod -lFatesUtilsMod -lWrapGreatCircleMod -o test_gc TestGreatCircle.F90

