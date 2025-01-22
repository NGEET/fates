# This file is required for unit testing, but is not used for production runs
list(APPEND fates_sources
  TwoStreamMLPEMod.F90
  FatesRadiationMemMod.F90
  FatesTwoStreamUtilsMod.F90
  )

sourcelist_to_parent(fates_sources)
