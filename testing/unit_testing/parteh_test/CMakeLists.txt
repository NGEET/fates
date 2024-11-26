set(pfunit_sources test_PRTAllometricCNP.pf)

add_pfunit_ctest(PRTAllometricCNP
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
