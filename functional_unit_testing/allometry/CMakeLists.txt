set(allom_sources FatesTestAllometry.F90)

add_executable(FATES_allom_exe ${allom_sources})

target_link_libraries(FATES_allom_exe
                  fates
                  csm_share)

add_test(allom_test FATES_allom_exe)

# Tell CTest how to figure out that "STOP 1" fails for the current
define_Fortran_stop_failure(allom_test)