# Local files
set(test_sources FatesUnitTestRadiation.F90)

#Add sources we need from other directories
set(sources_needed 
   FatesUnitTestIOMod.F90 
   FatesUnitTestOrbitalMod.F90
   FatesConstantsMod.F90)

extract_sources("${sources_needed}" "${utils_sources}" test_sources)

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

add_executable(FATES_rad_test ${test_sources})

add_test(rad_test FATES_rad_test)

# Tell CTest how to figure out that "STOP 1" fails for the current
define_Fortran_stop_failure(rad_test)
