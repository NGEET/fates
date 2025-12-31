set(pfunit_sources test_FireFuel.pf)

add_pfunit_ctest(FireFuel
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  