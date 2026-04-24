set(pfunit_sources test_FireWeather.pf)

add_pfunit_ctest(FireWeather
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)