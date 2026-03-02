set(pfunit_sources test_SortCohorts.pf)

add_pfunit_ctest(SortCohorts
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  