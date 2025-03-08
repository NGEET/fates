# This file is required for unit testing, but is not used for production runs
list(APPEND fates_sources
  SFParamsMod.F90
  SFFireWeatherMod.F90
  SFNesterovMod.F90
  FatesFuelMod.F90
  FatesFuelClassesMod.F90
  SFEquationsMod.F90
  )

sourcelist_to_parent(fates_sources)
