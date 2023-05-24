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

  list(APPEND fates_sources
  FatesConstantsMod.F90
  EDTypesMod.F90
  EDPftvarcon.F90
  EDParamsMod.F90
  FatesGlobals.F90 
  FatesInterfaceTypesMod.F90
  FatesParametersInterface.F90
  FatesParameterDerivedMod.F90
  FatesHydraulicsMemMod.F90
  FatesRunningMeanMod.F90
  FatesSizeAgeTypeIndicesMod.F90
  FatesIntegratorsMod.F90
  FatesUtilsMod.F90
  FatesSynchronizedParamsMod.F90
  ChecksBalancesMod.F90
  )

sourcelist_to_parent(clm_sources)
sourcelist_to_parent(fates_sources)
