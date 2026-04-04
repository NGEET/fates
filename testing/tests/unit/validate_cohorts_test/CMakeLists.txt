set(pfunit_sources test_ValidateCohorts.pf)

add_pfunit_ctest(ValidateCohorts
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  