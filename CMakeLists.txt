cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../share/cmake")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../components/cmeps/cmake")

FIND_PATH(NETCDFC_FOUND libnetcdf.a ${NETCDF_C_DIR}/lib)
FIND_PATH(NETCDFF_FOUND libnetcdff.a ${NETCDF_FORTRAN_DIR}/lib)
MESSAGE(" NETCDFC_FOUND = ${NETCDFC_FOUND}")
MESSAGE(" NETCDFF_FOUND = ${NETCDFF_FOUND}")

string(APPEND LDFLAGS " -L${NETCDF_FORTRAN_DIR}/lib -lnetcdff")
string(APPEND LDFLAGS " -L${NETCDF_C_DIR}/lib -lnetcdf")

include(CIME_initial_setup)

project(FATES_tests Fortran C)

include(CIME_utils)

set(HLM_ROOT "../../")

if (DEFINED ENV{ESMF_ROOT})
  list(APPEND CMAKE_MODULE_PATH $ENV{ESMF_ROOT}/cmake)
endif()
find_package(ESMF REQUIRED)

# This adds include directories needed for ESMF
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${ESMF_F90COMPILEPATHS} ")
# This (which is *not* done in the share CMakeLists.txt) adds all directories and
# libraries needed when linking ESMF, including any dependencies of ESMF. (But note that
# this does *not* include the "-lesmf" itself). In particular, note that this includes any
# link flags needed to link against PIO, which is needed on some systems (including
# derecho); bringing in these PIO-related link flags via this ESMF mechanism allows us to
# avoid explicitly including PIO as a link library, which wouldn't work on systems where
# there is no separate PIO library and instead ESMF is built with its internal PIO
# library.
link_libraries(${ESMF_INTERFACE_LINK_LIBRARIES})

# Add source directories from other share code (csm_share, etc.)
add_subdirectory(${HLM_ROOT}/share/src csm_share)
add_subdirectory(${HLM_ROOT}/share/unit_test_stubs/util csm_share_stubs)

# Add FATES source directories
add_subdirectory(${HLM_ROOT}/src/fates/main fates_main)
add_subdirectory(${HLM_ROOT}/src/fates/biogeochem fates_biogeochem)
add_subdirectory(${HLM_ROOT}/src/fates/biogeophys fates_biogeophys)
add_subdirectory(${HLM_ROOT}/src/fates/parteh fates_parteh)
add_subdirectory(${HLM_ROOT}/src/fates/fire fates_fire)
add_subdirectory(${HLM_ROOT}/src/fates/radiation fates_radiation)

# Testing directories
add_subdirectory(${HLM_ROOT}/src/fates/testing/tests/fortran_shr test_share)
add_subdirectory(${HLM_ROOT}/src/fates/testing/tests/functional/fire/shr fire_share)

# Remove shr_mpi_mod from share_sources. We are avoiding anything that depends on mpi in the unit test build.
remove_source_file("shr_mpi_mod.F90" share_sources)

# Remove some unneeded share sources that depend on shr_mpi_mod.
#
# CTSM handles this differently, by using a stub shr_mpi_mod in place of the real one. It
# needs this different approach because it needs these share sources in its unit tests.
# That would be another solution, which could be achieved via:
# add_subdirectory(${HLM_ROOT}/src/unit_test_stubs/csm_share csm_share_ctsm_stubs)
remove_source_file("shr_dust_emis_mod.F90" share_sources)
remove_source_file("shr_fire_emis_mod.F90" share_sources)

# Remove some other things from share_sources.
remove_source_file("shr_cal_mod.F90" share_sources)
remove_source_file("shr_pio_mod.F90" share_sources)

# Build libraries containing stuff needed for the unit tests.
# Eventually, these add_library calls should probably be distributed into the correct location, rather than being in this top-level CMakeLists.txt file.
add_library(csm_share ${share_sources})
declare_generated_dependencies(csm_share "${share_genf90_sources}")
add_library(fates ${fates_sources})
add_dependencies(fates csm_share ESMF)

# We need to look for header files here, in order to pick up shr_assert.h
include_directories(${HLM_ROOT}/share/include)

# This needs to be something we add dynamically
# via some calls using cime
set(NETCDF_C_DIR ${NETCDF_C_PATH})
set(NETCDF_FORTRAN_DIR ${NETCDF_F_PATH})

include_directories(${NETCDF_C_DIR}/include
                    ${NETCDF_FORTRAN_DIR}/include)
link_directories(${NETCDF_C_DIR}/lib
                ${NETCDF_FORTRAN_DIR}/lib)
                

# Tell cmake to look for libraries & mod files here, because this is where we built libraries
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Directories and libraries to include in the link step
link_directories(${CMAKE_CURRENT_BINARY_DIR})
link_libraries(esmf)
              
# Add the main test directory
add_subdirectory(${HLM_ROOT}/src/fates/testing)
