program FatesTestAllometry

  !use FatesAllometryMod,           only : h2d_allom
  !use PRTParametersMod,            only : prt_params
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader

  call param_reader%param_read()

  print *, "Hello, allometry"

end program FatesTestAllometry