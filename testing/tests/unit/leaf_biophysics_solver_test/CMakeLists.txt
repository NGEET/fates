set(pfunit_sources test_LeafBiophysicsSolvers.pf)

add_pfunit_ctest(LeafBiophysicsSolvers
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
