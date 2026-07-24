set(pfunit_sources test_QuadraticRoots.pf)

add_pfunit_ctest(QuadraticRoots
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  