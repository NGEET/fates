cmake_minimum_required(VERSION 2.8)

#list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
#include(CIME_initial_setup)

project(radiation_tests Fortran C)

#include(CIME_utils)

set(CLM_ROOT "../../")

find_package (NetCDF COMPONENTS F90)
include_directories(${NETCDF_F90_INCLUDE_DIRS})

message(" inc =  ${NETCDF_F90_INCLUDE_DIRS}  ")

#add_subdirectory(${CLM_ROOT}/src/fates/unit_testing/utils)
#add_subdirectory(${CLM_ROOT}/src/fates/main)
#add_subdirectory(${CLM_ROOT}/src/fates/unit_testing/test)

#add_library(utils ${utils_sources})

#target_link_libraries(utils, ${NETCDF_F90_INCLUDE_DIRS})
#include_directories(${CMAKE_CURRENT_BINARY_DIR})
#link_directories(${CMAKE_CURRENT_BINARY_DIR})