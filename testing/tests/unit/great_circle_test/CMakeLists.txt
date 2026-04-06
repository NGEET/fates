set(pfunit_sources test_GreatCircle.pf)

add_pfunit_ctest(GreatCircle
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  