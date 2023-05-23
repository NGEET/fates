# Append sources from this directory to the list.
list(APPEND fates_sources 
      FatesUnitTestIOMod.F90)

# Export source list to parent
sourcelist_to_parent(fates_sources)