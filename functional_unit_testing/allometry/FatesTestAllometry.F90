program FatesTestAllometry

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesAllometryMod,           only : h_allom
  use PRTParametersMod,            only : prt_params
  use FatesUnitTestParamReaderMod, only : fates_unit_test_param_reader

  implicit none

  ! LOCALS:
  type(fates_unit_test_param_reader) :: param_reader
  character(len=*), parameter        :: param_file = 'fates_params_default.nc'
  integer                            :: numpft
  integer                            :: i
  real(r8)                           :: height

  ! CONSTANTS:
  real(r8) :: min_dbh = 0.5_r8   ! minimum DBH to calculate [cm]
  real(r8) :: max_dbh = 200.0_r8 ! maximum DBH to calculate [cm]
  
  ! read in parameter file
  call param_reader%Init(param_file)
  call param_reader%RetrieveParameters()
  numpft = size(prt_params%wood_density, dim=1) - 1

  do i = 1, numpft
    call h_allom(25.0_r8, i, height)
    print *, height
  end do
  
end program FatesTestAllometry