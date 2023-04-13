cmake_minimum_required(VERSION 3.4)

project(radiation_tests Fortran C)

enable_language(Fortran)

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
include(CIME_initial_setup)
include(CIME_utils)

set(CLM_ROOT "../../")
add_subdirectory(${CLM_ROOT}/src/fates/main)
add_subdirectory(${CLM_ROOT}/src/fates/unit_testing/utils)
add_subdirectory(${CLM_ROOT}/src/fates/unit_testing/test)

set(NETCDF_C_DIR "/usr/local/Cellar/netcdf/4.9.2")
set(NETCDF_FORTRAN_DIR "/usr/local/Cellar/netcdf-fortran/4.6.0")
MESSAGE(">> CISM_NETCDF_C_DIR set to : ${NETCDF_C_DIR}")
MESSAGE(">> CISM_NETCDF_FORTRAN_DIR set to : ${NETCDF_FORTRAN_DIR}")

FIND_PATH(NETCDFC_FOUND libnetcdf.a ${NETCDF_C_DIR}/lib)
FIND_PATH(NETCDFF_FOUND libnetcdff.a ${NETCDF_FORTRAN_DIR}/lib)
MESSAGE(" NETCDFC_FOUND = ${NETCDFC_FOUND}")
MESSAGE(" NETCDFF_FOUND = ${NETCDFF_FOUND}")

IF (${NETCDFC_FOUND} STREQUAL ${NETCDF_C_DIR}/lib )
   SET(NETCDF_C_LIBS "netcdff;netcdf" CACHE STRING "Netcdf Library Names(s)")
ELSE()
   SET(NETCDF_C_LIBS "netcdf" CACHE STRING "Netcdf Library Names(s)") 
ENDIF()
MESSAGE(">> NETCDF_LIBS Library(s) set to : ${NETCDF_C_LIBS}")


include_directories(${NETCDF_C_DIR}/include
                    ${NETCDF_FORTRAN_DIR}/include)
link_directories(${NETCDF_C_DIR}/lib
                ${NETCDF_FORTRAN_DIR}/lib)

add_library(utils ${utils_sources})

include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})