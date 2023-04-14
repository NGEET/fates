# Local files
set(rad_sources FatesUnitTestRadiation.F90)

set(NETCDF_C_DIR "/usr/local/Cellar/netcdf/4.9.2")
set(NETCDF_FORTRAN_DIR "/usr/local/Cellar/netcdf-fortran/4.6.0")

FIND_PATH(NETCDFC_FOUND libnetcdf.a ${NETCDF_C_DIR}/lib)
FIND_PATH(NETCDFF_FOUND libnetcdff.a ${NETCDF_FORTRAN_DIR}/lib)
MESSAGE(" NETCDFC_FOUND = ${NETCDFC_FOUND}")
MESSAGE(" NETCDFF_FOUND = ${NETCDFF_FOUND}")

include_directories(${NETCDF_C_DIR}/include
                    ${NETCDF_FORTRAN_DIR}/include)
link_directories(${NETCDF_C_DIR}/lib
                ${NETCDF_FORTRAN_DIR}/lib)

add_executable(FATES_rad_exe ${rad_sources})

target_link_libraries(FATES_rad_exe PUBLIC
                  netcdff 
                  netcdf
                  fates
                  utils
                  csm_share)

target_include_directories(FATES_rad_exe PUBLIC
        "${PROJECT_BINARY_DIR}"
        "${PROJECT_SOURCE_DIR}/fates_main"
        )

add_test(rad_test FATES_rad_exe)

# Tell CTest how to figure out that "STOP 1" fails for the current
define_Fortran_stop_failure(rad_test)
