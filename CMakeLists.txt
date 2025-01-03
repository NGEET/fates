cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../share/cmake")

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
add_subdirectory(${HLM_ROOT}/src/fates/testing/testing_shr test_share)
add_subdirectory(${HLM_ROOT}/src/fates/testing/functional_testing/fire/shr fire_share)

# Remove shr_mpi_mod from share_sources.
# This is needed because we want to use the mock shr_mpi_mod in place of the real one
#
# TODO: this should be moved into a general-purpose function in Sourcelist_utils.
# Then this block of code could be replaced with a single call, like:
# remove_source_file(${share_sources} "shr_mpi_mod.F90")
foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_mpi_mod.F90" match_found ${sourcefile})
  if(match_found)
    list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()

# Remove shr_cal_mod from share_sources.
foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_cal_mod.F90" match_found ${sourcefile})
  if(match_found)
    list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()

# Remove shr_pio_mod from share_sources.
foreach (sourcefile ${share_sources})
  string(REGEX MATCH "shr_pio_mod.F90" match_found ${sourcefile})
  if(match_found)
    list(REMOVE_ITEM share_sources ${sourcefile})
  endif()
endforeach()

# Build libraries containing stuff needed for the unit tests.
# Eventually, these add_library calls should probably be distributed into the correct location, rather than being in this top-level CMakeLists.txt file.
add_library(csm_share ${share_sources})
declare_generated_dependencies(csm_share "${share_genf90_sources}")
add_library(fates ${fates_sources})
add_dependencies(fates csm_share)

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
              
# Add the main test directory
add_subdirectory(${HLM_ROOT}/src/fates/testing)
