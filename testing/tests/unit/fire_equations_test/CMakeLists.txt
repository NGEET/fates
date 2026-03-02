set(pfunit_sources test_FireEquations.pf)

add_pfunit_ctest(FireEquations
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  