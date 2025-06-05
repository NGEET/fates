set(pfunit_sources test_Moss.pf)

add_pfunit_ctest(Moss
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  