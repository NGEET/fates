set(pfunit_sources test_InsertCohort.pf)

add_pfunit_ctest(InsertCohorts
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  