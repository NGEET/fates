set(pfunit_sources test_CountCohorts.pf)

add_pfunit_ctest(CountCohorts
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  