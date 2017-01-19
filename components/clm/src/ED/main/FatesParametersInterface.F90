module FatesParametersInterface
  ! NOTE(bja, 2017-01) this is part of the interface between fates and
  ! the host model. To avoid circular dependancies, it should not
  ! depend on any host modules.

  use FatesConstantsMod, only : r8 => fates_r8
  
  implicit none

  integer, parameter, public :: max_params = 250
  integer, parameter, public :: max_dimensions = 2
  integer, parameter, public :: param_string_length = 40
  integer, parameter, public :: dimension_shape_scalar = 0
  integer, parameter, public :: dimension_shape_1d = 1
  integer, parameter, public :: dimension_shape_2d = 2

  ! FIXME(bja, 2017-01) these strings need to be changed to 'fates_'
  ! to namespace dimonsions and prevent name collisions if someone
  ! wants to write a single netcdf file containing host and fates
  ! parameters. Can't be done easily until this framework is being
  ! used to read variables.
  character(len=*), parameter, public :: dimension_name_scalar = 'scalar'
  character(len=*), parameter, public :: dimension_name_pft = 'pft'
  character(len=*), parameter, public :: dimension_name_segment = 'segment'
  
  type, private ::  parameter_type
     character(len=param_string_length) :: name
     logical :: host_parameter
     integer :: dimension_shape
     character(len=param_string_length) :: dimension_names(max_dimensions)
     real(r8), allocatable :: data(:, :)
  end type parameter_type

  type, public :: fates_parameters_type
     integer, private :: num_parameters
     type(parameter_type), private :: parameters(max_params)

   contains
     procedure, public :: Init
     procedure, public :: RegisterParameter
     generic, public :: RetreiveParameter => RetreiveParameterScalar, RetreiveParameter1D, RetreiveParameter2D
     generic, public :: SetData => SetDataScalar, SetData1D, SetData2D
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
  subroutine RegisterParameter(this, name, dimension_shape, dimension_names, host_parameter)
    
    implicit none

    class(fates_parameters_type), intent(inout) :: this
    character(len=param_string_length), intent(in) :: name
    integer, intent(in) :: dimension_shape
    character(len=param_string_length) :: dimension_names(1:)
    logical, intent(in), optional :: host_parameter

    integer :: i, n, num_names
    
    this%num_parameters = this%num_parameters + 1
    i = this%num_parameters
    ! FIXME(bja, 2017-01) assert(i <= max_params)
    this%parameters(i)%name = name
    this%parameters(i)%dimension_shape = dimension_shape
    ! FIXME(bja, 2017-01) assert(size(dimension_names, 1) <= max_dimensions)
    num_names = min(max_dimensions, size(dimension_names, 1))
    do n = 1, num_names
       this%parameters(i)%dimension_names(n) = dimension_names(n)
    end do
    this%parameters(i)%host_parameter = .false.
    if (present(host_parameter)) then
       this%parameters(i)%host_parameter = .true.
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

    implicit none

    class(fates_parameters_type), intent(inout) :: this
    character(len=param_string_length), intent(in) :: name
    real(r8), intent(out) :: data(:)

    integer :: i

    i = this%FindIndex(name)
    ! assert(size(data) == size(this%parameters(i)%data))
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
  subroutine GetMetaData(this, index, name, dimension_shape)

    implicit none
    
    class(fates_parameters_type), intent(in) :: this
    integer, intent(in) :: index
    character(len=param_string_length), intent(out) :: name
    integer, intent(out) :: dimension_shape

    name = this%parameters(index)%name
    dimension_shape = this%parameters(index)%dimension_shape

  end subroutine GetMetaData
  
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
    ! FIXME(bja, 2017-01) this is broken, needs data dimensions to work correctly!

    implicit none
    
    class(fates_parameters_type), intent(inout) :: this
    integer, intent(in) :: index
    real(r8), intent(in) :: data(:)

    allocate(this%parameters(index)%data(size(data), 1))
    this%parameters(index)%data(:, 1) = data

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



