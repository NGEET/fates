set(fire_sources FatesUnitTestSF.F90 FatesUnitTestSFMod.F90)

set(NETCDF_C_DIR "/usr/local/Cellar/netcdf/4.9.2_1")
set(NETCDF_FORTRAN_DIR "/usr/local/Cellar/netcdf-fortran/4.6.1")

FIND_PATH(NETCDFC_FOUND libnetcdf.a ${NETCDF_C_DIR}/lib)
FIND_PATH(NETCDFF_FOUND libnetcdff.a ${NETCDF_FORTRAN_DIR}/lib)
MESSAGE(" NETCDFC_FOUND = ${NETCDFC_FOUND}")
MESSAGE(" NETCDFF_FOUND = ${NETCDFF_FOUND}")

include_directories(${NETCDF_C_DIR}/include
                    ${NETCDF_FORTRAN_DIR}/include)
link_directories(${NETCDF_C_DIR}/lib
                ${NETCDF_FORTRAN_DIR}/lib
                ${PFUNIT_TOP_DIR}/lib)

add_executable(FATES_fire_exe ${fire_sources})

target_link_libraries(FATES_fire_exe 
                  netcdff 
                  netcdf
                  fates
                  csm_share
                  funit)

add_test(fire_test FATES_fire_exe)

# Tell CTest how to figure out that "STOP 1" fails for the current
define_Fortran_stop_failure(fire_test)

