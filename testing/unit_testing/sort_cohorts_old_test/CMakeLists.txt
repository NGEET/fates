set(pfunit_sources test_SortCohortsOld.pf)

add_pfunit_ctest(SortCohortsOld
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  