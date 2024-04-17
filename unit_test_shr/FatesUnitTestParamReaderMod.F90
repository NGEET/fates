module FatesUnitTestParamReaderMod

  use FatesConstantsMod,          only : r8 => fates_r8
  use FatesParametersInterface,   only : fates_param_reader_type
  use FatesParametersInterface,   only : fates_parameters_type
  use FatesParametersInterface,   only : param_string_length
  use FatesParametersInterface,   only : max_dimensions, max_used_dimensions
  use FatesParametersInterface,   only : dimension_shape_scalar, dimension_shape_1d, dimension_shape_2d
  use EDParamsMod,                only : FatesRegisterParams, FatesReceiveParams
  use SFParamsMod,                only : SpitFireRegisterParams, SpitFireReceiveParams
  use PRTInitParamsFatesMod,      only : PRTRegisterParams, PRTReceiveParams
  use PRTParametersMod,           only : prt_params
  use FatesParameterDerivedMod,   only : param_derived
  use FatesSynchronizedParamsMod, only : FatesSynchronizedParamsInst
  use EDPftvarcon,                only : EDPftvarcon_inst
  use FatesUnitTestIOMod,         only : OpenNCFile, GetDimID, GetDimLen, GetVar, CloseNCFile

  implicit none
  private

  type, public, extends(fates_param_reader_type) :: fates_unit_test_param_reader

    character(:), allocatable :: filename ! local file name of parameter
    
    contains
      procedure, public :: Read => ReadParameters
      procedure, public :: Init 
      procedure, public :: RetrieveParameters  

  end type fates_unit_test_param_reader

  contains

  subroutine Init(this, param_file)
    !
    ! DESCRIPTION:
    ! Initialize the parameter reader class
    !
    ! ARGUMENTS:
    class(fates_unit_test_param_reader) :: this
    character(len=*)                    :: param_file 

    this%filename = trim(param_file)

  end subroutine Init

  ! --------------------------------------------------------------------------------------

  subroutine ReadParameters(this, fates_params)
    !
    ! DESCRIPTION:
    ! Read 'fates_params' parameters from storage
    !
    ! ARGUMENTS:
    class(fates_unit_test_param_reader)         :: this
    class(fates_parameters_type), intent(inout) :: fates_params

    ! LOCALS:
    real(r8), allocatable              :: data2d(:, :)
    real(r8), allocatable              :: data1d(:)
    real(r8)                           :: data_scalar
    integer                            :: ncid
    integer                            :: num_params
    integer                            :: dimension_shape
    integer                            :: i
    integer                            :: max_dim_size
    character(len=param_string_length) :: name
    integer                            :: dimension_sizes(max_dimensions)
    character(len=param_string_length) :: dimension_names(max_dimensions)
    logical                            :: is_host_param

    call OpenNCFile(this%filename, ncid, 'read')
    call SetParameterDimensions(ncid, fates_params)

    num_params = fates_params%num_params()
    do i = 1, num_params
      call fates_params%GetMetaData(i, name, dimension_shape, dimension_sizes,           &
        dimension_names, is_host_param)
        select case(dimension_shape)
          case(dimension_shape_scalar)
            call GetVar(ncid, name, data_scalar)
            call fates_params%SetData(i, data_scalar)
          case(dimension_shape_1d)
            call GetVar(ncid, name, data1d)
            call fates_params%SetData(i, data1d)
          case(dimension_shape_2d)
            call GetVar(ncid, name, data2d)
            call fates_params%SetData(i, data2d)
          case default
            write(*, '(a,a)') 'dimension shape:', dimension_shape
            write(*, '(a)') 'unsupported number of dimensions reading parameters.'
            stop
        end select
    end do

    if (allocated(data1d)) deallocate(data1d)
    if (allocated(data2d)) deallocate(data2d)

    call CloseNCFile(ncid)

  end subroutine ReadParameters

  ! --------------------------------------------------------------------------------------

  subroutine RetrieveParameters(this)
    !
    ! DESCRIPTION:
    ! Read in fates parameters
    !
    ! ARGUMENTS:
    class(fates_unit_test_param_reader), intent(in) :: this 
    
    ! LOCALS:
    class(fates_parameters_type), allocatable :: fates_params
    class(fates_parameters_type), allocatable :: fates_pft_params

    ! allocate and read in parameters
    allocate(fates_params)
    allocate(fates_pft_params)
    call fates_params%Init()
    call fates_pft_params%Init()

    call EDPftvarcon_inst%Init()
    
    call PRTRegisterParams(fates_params)
    call FatesRegisterParams(fates_params)  
    call SpitFireRegisterParams(fates_params) 
    call FatesSynchronizedParamsInst%RegisterParams(fates_params)
    call EDPftvarcon_inst%Register(fates_pft_params)

    call this%Read(fates_params)
    call this%Read(fates_pft_params)

    call FatesReceiveParams(fates_params)
    call SpitFireReceiveParams(fates_params)
    call PRTReceiveParams(fates_params)
    call FatesSynchronizedParamsInst%ReceiveParams(fates_params)
    call EDPftvarcon_inst%Receive(fates_pft_params)

    call fates_params%Destroy()
    call fates_pft_params%Destroy()
    deallocate(fates_params)
    deallocate(fates_pft_params)

    ! initialize derived parameters
    call param_derived%Init(size(prt_params%wood_density, dim=1))

  end subroutine RetrieveParameters

  ! --------------------------------------------------------------------------------------

  subroutine SetParameterDimensions(ncid, fates_params)
    !
    ! DESCRIPTION:
    ! Read in fates parameters
    !
    ! ARGUMENTS:
    integer,                      intent(in)    :: ncid         ! netcdf file ID
    class(fates_parameters_type), intent(inout) :: fates_params ! fates parameters class

    ! LOCALS:
    integer                            :: num_used_dimensions
    character(len=param_string_length) :: used_dimension_names(max_used_dimensions)
    integer                            :: used_dimension_sizes(max_used_dimensions)

    call fates_params%GetUsedDimensions(.false., num_used_dimensions, used_dimension_names)

    call GetUsedDimensionSizes(ncid, num_used_dimensions, used_dimension_names,          &
      used_dimension_sizes)

    call fates_params%SetDimensionSizes(.false., num_used_dimensions,                    &
      used_dimension_names, used_dimension_sizes)

  end subroutine SetParameterDimensions

  ! --------------------------------------------------------------------------------------

  subroutine GetUsedDimensionSizes(ncid, num_used_dimensions, dimension_names, dimension_sizes)
    !
    ! DESCRIPTION:
    ! Get dimension sizes for parameters
    !
    ! ARGUMENTS:
    integer,                            intent(in)  :: ncid
    integer,                            intent(in)  :: num_used_dimensions
    character(len=param_string_length), intent(in)  :: dimension_names(:)
    integer,                            intent(out) :: dimension_sizes(:)
 
    ! LOCALS
    integer :: d
    integer :: dim_id
 
    dimension_sizes(:) = 0
    do d = 1, num_used_dimensions
      call GetDimID(ncid, dimension_names(d), dim_id)
      call GetDimLen(ncid, dim_id, dimension_sizes(d))
    end do
 
  end subroutine GetUsedDimensionSizes

end module FatesUnitTestParamReaderMod