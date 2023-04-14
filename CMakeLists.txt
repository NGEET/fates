cmake_minimum_required(VERSION 3.4)

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
include(CIME_initial_setup)

project(FATES_tests Fortran C)

enable_language(Fortran)

include(CIME_utils)

set(CLM_ROOT "../../")

# Add source directories from other share code (csm_share, etc.)
add_subdirectory(${CLM_ROOT}/share/src csm_share)

# Add FATES source directories
add_subdirectory(${CLM_ROOT}/src/fates/main fates_main)
add_subdirectory(${CLM_ROOT}/src/fates/parteh fates_parteh)
add_subdirectory(${CLM_ROOT}/src/fates/unit_testing/utils fates_test_utils)

# Build libraries containing stuff needed for the unit tests.
# Eventually, these add_library calls should probably be distributed into the correct location, rather than being in this top-level CMakeLists.txt file.
add_library(csm_share ${share_sources})
declare_generated_dependencies(csm_share "${share_genf90_sources}")

add_library(fates ${fates_sources})
add_library(utils ${utils_sources})
add_dependencies(fates csm_share)

# We need to look for header files here, in order to pick up shr_assert.h
include_directories(${CLM_ROOT}/share/include)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

# Add the test directories
add_subdirectory(${CLM_ROOT}/src/fates/unit_testing/radiation_test fates_testing)
