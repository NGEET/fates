module FatesUnitTestParamReaderMod

  use FatesParametersInterface, only : fates_param_reader_type
  use FatesParametersInterface, only : fates_parameters_type
  use PRTInitParamsFatesMod,    only : PRTRegisterParams

  implicit none
  private

  type, public, extends(fates_param_reader_type) :: fates_unit_test_param_reader

    character(:), allocatable :: filename ! local file name of parameter
    
    contains
      procedure, public :: Read => read_parameters
      procedure, public :: param_read  

  end type fates_unit_test_param_reader

  contains 

  subroutine read_parameters(this, fates_params)
    !
    ! DESCRIPTION:
    ! Read 'fates_params' parameters from storage
    !
    ! ARGUMENTS:
    class(fates_unit_test_param_reader)         :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    print *, "read in parameters"

    !call ParametersFromNetCDF(fates_paramfile, is_host_file, fates_params)

  end subroutine read_parameters

  ! --------------------------------------------------------------------------------------

  subroutine param_read(this)
    !
    ! DESCRIPTION:
    ! Read in fates parameters
    !
    ! ARGUMENTS:
    class(fates_unit_test_param_reader), intent(in) :: this 
    
    ! LOCALS:
    !class(fates_parameters_type), allocatable :: fates_params

    ! allocate and read in parameters
    !allocate(fates_params)
    !call fates_params%Init()
    !call PRTRegisterParams(fates_params)

    !call this%Read(fates_params)

  end subroutine param_read

end module FatesUnitTestParamReaderMod