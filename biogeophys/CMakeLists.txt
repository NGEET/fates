# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND fates_sources
  FatesHydroWTFMod.F90
  )

sourcelist_to_parent(fates_sources)
