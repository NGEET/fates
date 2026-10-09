set(pfunit_sources test_EDAccumulateFluxes.pf)

add_pfunit_ctest(EDAccumulateFluxes
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
