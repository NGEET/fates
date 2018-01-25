# Note that this is just used for unit testing; hence, we only need to add
# source files that are currently used in unit tests

list(APPEND clm_sources
  FatesGlobals.F90
  EDTypesMod.F90
  EDPftvarcon.F90
  FatesConstantsMod.F90
  FatesHydraulicsMemMod.F90
  FatesParametersInterface.F90
  FatesUtilsMod.F90
  )

sourcelist_to_parent(clm_sources)

