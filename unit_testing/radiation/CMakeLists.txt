# Local files
set(test_sources FatesUnitTestRadiation.F90)

# Add sources we need from other directories
set(sources_needed 
    FatesUnitTestIOMod.F90 
    FatesUnitTestOrbitalMod.F90)
extract_sources("${sources_needed}" "${utils_sources}" test_sources)

add_executable(FATES_rad_test ${test_sources})
