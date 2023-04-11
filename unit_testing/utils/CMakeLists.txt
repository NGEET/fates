# Append sources from this directory to the list.
list(APPEND utils_sources FatesUnitTestIOMod.F90 FatesUnitTestOrbitalMod.F90)

# Export source list to parent
sourcelist_to_parent(circle_area_sources)