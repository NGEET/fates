module FatesParametersInterface
  ! NOTE(bja, 2017-01) this is part of the interface between fates and
  ! the host model. To avoid circular dependancies, it should not
  ! depend on any host modules.

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals, only : fates_log
  
  implicit none

  integer, parameter, public :: max_params = 250
  integer, parameter, public :: max_dimensions = 2
  integer, parameter, public :: max_used_dimensions = 25
  integer, parameter, public :: param_string_length = 40
  integer, parameter, public :: dimension_shape_scalar = 0
  integer, parameter, public :: dimension_shape_1d = 1
  integer, parameter, public :: dimension_shape_2d = 2

  ! FIXME(bja, 2017-01) these strings need to be changed to 'fates_'
  ! to namespace dimonsions and prevent name collisions if someone
  ! wants to write a single netcdf file containing host and fates
  ! parameters. Can't be done easily until this framework is being
  ! used to read variables.
  ! FIXME(bja, 2017-01) change 'param' to 'scalar'!
  character(len=*), parameter, public :: dimension_name_scalar = 'param'
  character(len=*), parameter, public :: dimension_name_pft = 'pft'
  character(len=*), parameter, public :: dimension_name_segment = 'segment'
  character(len=*), parameter, public :: dimension_name_cwd = 'NCWD'
  character(len=*), parameter, public :: dimension_name_lsc = 'litterclass'
  character(len=*), parameter, public :: dimension_name_fsc = 'litterclass'
  character(len=*), parameter, public :: dimension_name_allpfts = 'allpfts'
  
  type, private ::  parameter_type
     character(len=param_string_length) :: name
     logical :: sync_with_host
     integer :: dimension_shape
     integer :: dimension_sizes(max_dimensions)
     character(len=param_string_length) :: dimension_names(max_dimensions)
     real(r8), allocatable :: data(:, :)
  end type parameter_type

  type, public :: fates_parameters_type
     integer, private :: num_parameters
     type(parameter_type), private :: parameters(max_params)

   contains
     procedure, public :: Init
     procedure, public :: Destroy
     procedure, public :: RegisterParameter
     generic, public :: RetreiveParameter => RetreiveParameterScalar, RetreiveParameter1D, RetreiveParameter2D
     generic, public :: SetData => SetDataScalar, SetData1D, SetData2D
     procedure, public :: GetUsedDimensions
     procedure, public :: SetDimensionSizes
     procedure, public :: GetMaxDimensionSize
     procedure, public :: GetMetaData
     procedure, public :: num_params

     procedure, private :: RetreiveParameterScalar
     procedure, private :: RetreiveParameter1D
     procedure, private :: RetreiveParameter2D
     procedure, private :: SetDataScalar
     procedure, private :: SetData1D
     procedure, private :: SetData2D
     procedure, private :: FindIndex
     
  end type fates_parameters_type

contains

  !-----------------------------------------------------------------------
  subroutine Init(this)

    implicit none

    class(fates_parameters_type), intent(inout) :: this

    this%num_parameters = 0

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine Destroy(this)

    implicit none

    class(fates_parameters_type), intent(inout) :: this

    integer :: n
    do n = 1, this%num_parameters
       deallocate(this%parameters(n)%data)
    end do

  end subroutine Destroy

  !-----------------------------------------------------------------------
  subroutine RegisterParameter(this, name, dimension_shape, dimension_names, sync_with_host)
    
    implicit none

    class(fates_parameters_type), intent(inout) :: this
    character(len=param_string_length), intent(in) :: name
    integer, intent(in) :: dimension_shape
    character(len=param_string_length) :: dimension_names(1:)
    logical, intent(in), optional :: sync_with_host

    integer :: i, n, num_names
    
    this%num_parameters = this%num_parameters + 1
    i = this%num_parameters
    ! FIXME(bja, 2017-01) assert(i <= max_params)
    this%parameters(i)%name = name
    this%parameters(i)%dimension_shape = dimension_shape
    this%parameters(i)%dimension_sizes(:) = 0
    ! FIXME(bja, 2017-01) assert(size(dimension_names, 1) <= max_dimensions)
    num_names = min(max_dimensions, size(dimension_names, 1))
    this%parameters(i)%dimension_names(:) = ''
    do n = 1, num_names
       this%parameters(i)%dimension_names(n) = dimension_names(n)
    end do
    this%parameters(i)%sync_with_host = .false.
    if (present(sync_with_host)) then
       this%parameters(i)%sync_with_host = sync_with_host
    end if
    
  end subroutine RegisterParameter

  !-----------------------------------------------------------------------
  subroutine RetreiveParameterScalar(this, name, data)

    implicit none

    class(fates_parameters_type), intent(inout) :: this
    character(len=param_string_length), intent(in) :: name
    real(r8), intent(out) :: data

    integer :: i

    i = this%FindIndex(name)
    ! assert(size(data) == size(this%parameters(i)%data))
    data = this%parameters(i)%data(1, 1)

  end subroutine RetreiveParameterScalar

  !-----------------------------------------------------------------------
  subroutine RetreiveParameter1D(this, name, data)

    use abortutils, only : endrun

    implicit none

    class(fates_parameters_type), intent(inout) :: this
    character(len=param_string_length), intent(in) :: name
    real(r8), intent(out) :: data(:)

    integer :: i, d, size_dim_1

    i = this%FindIndex(name)
    if (size(data) /= size(this%parameters(i)%data(:, 1))) then
       write(fates_log(), *) 'ERROR : retreiveparameter1d : ', name, ' size inconsistent.'
       write(fates_log(), *) 'ERROR : size expected size = ', size(data)
       write(fates_log(), *) 'ERROR : data size received from file = ', size(this%parameters(i)%data(:, 1))
       write(fates_log(), *) 'ERROR : dimesions received from file'
       write(fates_log(), *) 'ERROR : names    size'
       do d = 1, max_dimensions
          write(fates_log(), *) this%parameters(i)%dimension_names(d), ', ', this%parameters(i)%dimension_sizes(d)
       end do
       call endrun(msg='size error retreiving 1d parameter.')
    end if
    data = this%parameters(i)%data(:, 1)

  end subroutine RetreiveParameter1D

  !-----------------------------------------------------------------------
  subroutine RetreiveParameter2D(this, name, data)

    implicit none

    class(fates_parameters_type), intent(inout) :: this
    character(len=param_string_length), intent(in) :: name
    real(r8), intent(out) :: data(:, :)

    integer :: i

    i = this%FindIndex(name)
    ! assert(size(data) == size(this%parameters(i)%data))
    data = this%parameters(i)%data

  end subroutine RetreiveParameter2D

  !-----------------------------------------------------------------------
  function FindIndex(this, name) result(i)

    implicit none

    class(fates_parameters_type), intent(in) :: this
    character(len=param_string_length), intent(in) :: name

    integer :: i
    
    do i = 1, this%num_parameters
       if (trim(this%parameters(i)%name) == trim(name)) then
          exit
       end if
    end do
    if (i > this%num_parameters) then
       ! error, parameter name not found.
    end if
    
  end function FindIndex

  !-----------------------------------------------------------------------
  integer function num_params(this)

    implicit none
    
    class(fates_parameters_type), intent(in) :: this

    num_params = this%num_parameters

  end function num_params
  
  !-----------------------------------------------------------------------
  subroutine GetUsedDimensions(this, is_host_file, num_used_dimensions, used_dimensions)
    ! Construct a list of the unique dimension names used by the
    ! parameters.
    
    implicit none
    
    class(fates_parameters_type), intent(inout) :: this
    logical, intent(in) :: is_host_file
    integer, intent(out) :: num_used_dimensions
    character(len=param_string_length), intent(out) :: used_dimensions(max_used_dimensions)

    integer :: p, d, i
    character(len=param_string_length) :: dim_name

    num_used_dimensions = 0
    do p = 1, this%num_parameters
       if (is_host_file .eqv. this%parameters(p)%sync_with_host) then
          do d = 1, max_dimensions
             dim_name = this%parameters(p)%dimension_names(d)
             if (len_trim(dim_name) /= 0) then
                ! non-empty dimension name, check if it needs to be added to the list.
                do i = 1, num_used_dimensions
                   if (used_dimensions(i) == dim_name) then
                      ! dimension is already in list. can stop searching
                      exit
                   end if
                end do

                if (i > num_used_dimensions) then
                   ! dimension name was not in the list, add it.
                   num_used_dimensions = num_used_dimensions + 1
                   used_dimensions(num_used_dimensions) = dim_name
                end if
             end if ! if dim_name
          end do ! do d
       end if ! if host_param
    end do ! do p

  end subroutine GetUsedDimensions
  
  !-----------------------------------------------------------------------
  subroutine SetDimensionSizes(this, is_host_file, num_used_dimensions, dimension_names, dimension_sizes)
    ! Construct a list of the unique dimension names used by the
    ! parameters.
    
    implicit none
    
    class(fates_parameters_type), intent(inout) :: this
    logical, intent(in) :: is_host_file
    integer, intent(in) :: num_used_dimensions
    character(len=param_string_length), intent(in) :: dimension_names(max_used_dimensions)
    integer, intent(in) :: dimension_sizes(max_used_dimensions)

    integer :: p, d, i
    character(len=param_string_length) :: dim_name

    do p = 1, this%num_parameters
       if (is_host_file .eqv. this%parameters(p)%sync_with_host) then
          do d = 1, max_dimensions
             dim_name = this%parameters(p)%dimension_names(d)
             if (len_trim(dim_name) /= 0) then
                ! non-empty dimension name, set the size
                do i = 1, num_used_dimensions
                   if (trim(dimension_names(i)) == trim(dim_name)) then
                      !write(*, *) '--> ', trim(this%parameters(p)%name), ' setting ', trim(dim_name), ' d = ', d, 'size = ', dimension_sizes(i)
                      this%parameters(p)%dimension_sizes(d) = dimension_sizes(i)
                      exit
                   end if
                end do
             end if ! if dim_name
          end do ! do dim
       end if ! if host_param
    end do ! do param

  end subroutine SetDimensionSizes
  
  !-----------------------------------------------------------------------
  subroutine GetMetaData(this, index, name, dimension_shape, dimension_sizes, is_host_param)

    implicit none
    
    class(fates_parameters_type), intent(in) :: this
    integer, intent(in) :: index
    character(len=param_string_length), intent(out) :: name
    integer, intent(out) :: dimension_shape
    integer, intent(out) :: dimension_sizes(max_dimensions)
    logical, intent(out) :: is_host_param

    name = this%parameters(index)%name
    dimension_shape = this%parameters(index)%dimension_shape
    dimension_sizes = this%parameters(index)%dimension_sizes
    is_host_param = this%parameters(index)%sync_with_host

  end subroutine GetMetaData
  
  !-----------------------------------------------------------------------
  function GetMaxDimensionSize(this) result(max_dim_size)

    implicit none

    class(fates_parameters_type), intent(in) :: this

    integer :: p, d, max_dim_size

    max_dim_size = 0

    do p = 1, this%num_params()
       do d = 1, max_dimensions
          max_dim_size = max(max_dim_size, this%parameters(p)%dimension_sizes(d))
       end do
    end do
    
  end function GetMaxDimensionSize

  !-----------------------------------------------------------------------
  subroutine SetDataScalar(this, index, data)

    implicit none
    
    class(fates_parameters_type), intent(inout) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: data

    allocate(this%parameters(index)%data(1, 1))
    this%parameters(index)%data(1, 1) = data

  end subroutine SetDataScalar
  
  !-----------------------------------------------------------------------
  subroutine SetData1D(this, index, data)

    use abortutils, only : endrun
    
    implicit none
    
    class(fates_parameters_type), intent(inout) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: data(:)

    integer :: size_dim_1, d

    size_dim_1 = this%parameters(index)%dimension_sizes(1)
    if (size(data) /= size_dim_1) then
       write(fates_log(), *) 'ERROR : setdata1d : ', this%parameters(index)%name, ' size inconsistent.'
       write(fates_log(), *) 'ERROR : expected size = ', size(data)
       write(fates_log(), *) 'ERROR : data size received from file = ', size_dim_1
       write(fates_log(), *) 'ERROR : dimesions received from file'
       write(fates_log(), *) 'ERROR : names    size'
       do d = 1, max_dimensions
          write(fates_log(), *) this%parameters(index)%dimension_names(d), ', ', this%parameters(index)%dimension_sizes(d)
       end do
       call endrun(msg='size error setting 1d parameter.')
    end if

    allocate(this%parameters(index)%data(size_dim_1, 1))
    this%parameters(index)%data(:, 1) = data(:)

  end subroutine SetData1D

  !-----------------------------------------------------------------------
  subroutine SetData2D(this, index, data)
    ! FIXME(bja, 2017-01) this is broken, needs data dimensions to work correctly!

    implicit none
    
    class(fates_parameters_type), intent(inout) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: data(:, :)

    ! NOTE(bja, 2017-01) This should work for fortran 2003? Or 2008?
    ! Either way, it works with intel and pgi being used in 2017-01,
    ! but is broken in gfortran 5.2 and earlier. That would copy the
    ! data as well....

    !X! allocate(this%parameters(index)%data, source=data)
    
    allocate(this%parameters(index)%data(size(data, 1), size(data, 2)))
    this%parameters(index)%data = data

  end subroutine SetData2D
end module FatesParametersInterface



