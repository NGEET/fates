program FatesTestAllometry

  !use FatesAllometryMod,           only : h2d_allom
  use PRTParametersMod,            only : prt_params
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader
  character(len=*), parameter        :: param_file = 'fates_params_default.nc'
  
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  

end program FatesTestAllometry