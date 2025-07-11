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
  use FatesInterfaceTypesMod,     only : nleafage
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
    real(r8), allocatable              :: data2d(:, :)                    ! data for 2D parameters
    real(r8), allocatable              :: data1d(:)                       ! data for 1D parameters
    real(r8)                           :: data_scalar                     ! data for scalar parameters
    integer                            :: ncid                            ! netcdf file ID
    integer                            :: num_params                      ! total number of parameters
    integer                            :: dimension_shape                 ! shape of parameter's dimension
    integer                            :: i                               ! looping index
    character(len=param_string_length) :: name                            ! parameter name
    integer                            :: dimension_sizes(max_dimensions) ! sizes of dimensions from parameter file
    character(len=param_string_length) :: dimension_names(max_dimensions) ! names of dimensions from parameter file 
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
    class(fates_unit_test_param_reader), intent(in) :: this        ! parameter reader class
     
    ! LOCALS:
    class(fates_parameters_type), allocatable :: fates_params      ! fates parameters (for non-pft parameters)
    class(fates_parameters_type), allocatable :: fates_pft_params  ! fates parameters (for pft parameters)

    ! allocate and read in parameters
    allocate(fates_params)
    allocate(fates_pft_params)
    call fates_params%Init()
    call fates_pft_params%Init()

    call EDPftvarcon_inst%Init()
    
    call FatesRegisterParams(fates_params)
    call SpitFireRegisterParams(fates_params) 
    call PRTRegisterParams(fates_params)
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
    
    nleafage = size(prt_params%leaf_long, dim=2)
  
    ! initialize derived parameters
    call param_derived%Init(size(prt_params%wood_density, dim=1))

  end subroutine RetrieveParameters

  ! --------------------------------------------------------------------------------------

  subroutine SetParameterDimensions(ncid, fates_params)
    !
    ! DESCRIPTION:
    ! Gets and sets the parameter dimensions for the fates parameters class
    !
    ! ARGUMENTS:
    integer,                      intent(in)    :: ncid         ! netcdf file ID
    class(fates_parameters_type), intent(inout) :: fates_params ! fates parameters class

    ! LOCALS:
    integer                            :: num_used_dimensions                       ! total number of dimensions
    character(len=param_string_length) :: used_dimension_names(max_used_dimensions) ! dimension names
    integer                            :: used_dimension_sizes(max_used_dimensions) ! dimension sizes

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
    ! Gets dimension sizes for parameters
    !
    ! ARGUMENTS:
    integer,                            intent(in)  :: ncid                ! netcdf file ID
    integer,                            intent(in)  :: num_used_dimensions ! number of dimensions
    character(len=param_string_length), intent(in)  :: dimension_names(:)  ! dimension names
    integer,                            intent(out) :: dimension_sizes(:)  ! dimension sizes
 
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