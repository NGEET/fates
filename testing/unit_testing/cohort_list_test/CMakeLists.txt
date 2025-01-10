set(pfunit_sources test_CohortLinkedList.pf)

add_pfunit_ctest(CohortLikedList
  TEST_SOURCES "${pfunit_sources}"
  LINK_LIBRARIES fates csm_share)
  