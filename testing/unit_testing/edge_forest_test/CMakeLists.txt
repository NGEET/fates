set(pfunit_sources test_EdgeForest.pf)

add_pfunit_ctest(EdgeForest
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
