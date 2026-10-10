add_pfunit_ctest(FatesPlantHydraulicsMod_drive_test 
  TEST_SOURCES FatesPlantHydraulicsMod_drive_test.pf
  LINK_LIBRARIES fates csm_share
)

add_pfunit_ctest(FatesPlantHydraulicsMod_sizedep_test 
  TEST_SOURCES FatesPlantHydraulicsMod_sizedep_test.pf
  LINK_LIBRARIES fates csm_share
)

add_pfunit_ctest(FatesPlantHydraulicsMod_geom_test 
  TEST_SOURCES FatesPlantHydraulicsMod_geom_test.pf
  LINK_LIBRARIES fates csm_share
)

add_pfunit_ctest(FatesPlantHydraulicsMod_cohort_test 
  TEST_SOURCES FatesPlantHydraulicsMod_cohort_test.pf
  LINK_LIBRARIES fates csm_share
)
