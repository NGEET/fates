set(pfunit_sources test_Ecotypes.pf)

add_pfunit_ctest(Ecotypes
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
