set(pfunit_sources test_InsertCohortOld.pf)

add_pfunit_ctest(InsertCohortsOld
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  