set(pfunit_sources test_HelloWorld.pf)

add_pfunit_ctest(HelloWorld
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  