set(pfunit_sources test_Hello.pf)

add_pfunit_ctest(Hello
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  