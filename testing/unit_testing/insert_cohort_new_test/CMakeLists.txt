set(pfunit_sources test_InsertCohortNew.pf)

add_pfunit_ctest(InsertCohortsNew
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  