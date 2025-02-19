# This file is required for unit testing, but is not used for production runs
list(APPEND fates_sources
  FatesLitterMod.F90
  FatesCohortMod.F90
  FatesAllometryMod.F90
  DamageMainMod.F90
  EDCohortDynamicsMod.F90
  FatesPatchMod.F90)

sourcelist_to_parent(fates_sources)
