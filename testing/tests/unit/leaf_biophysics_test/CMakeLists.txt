set(pfunit_sources test_LeafBiophysics.pf)

add_pfunit_ctest(LeafBiophysics
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  