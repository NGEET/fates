module FatesUnitTestIOMod
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesGlobals,      only : fates_endrun
  use shr_kind_mod,      only : SHR_KIND_CL
  use netcdf

  implicit none
  private

  ! LOCALS
  integer, public, parameter :: type_double = 1 ! type
  integer, public, parameter :: type_int = 2    ! type
  integer, public, parameter :: type_char = 3   ! type
  
  interface GetVar
    module procedure GetVarScalarReal
    module procedure GetVar1DReal
    module procedure GetVar2DReal
    module procedure GetVar3DReal
    module procedure GetVar1DInt
    module procedure GetVar2DInt
    module procedure GetVar3DInt
  end interface

  interface WriteVar
    module procedure WriteVar1DReal
    module procedure WriteVar2DReal
    module procedure WriteVar1DInt
    module procedure WriteVar2DInt
    module procedure WriteVar1DChar
    module procedure WriteVar2DChar
  end interface

  public :: OpenNCFile
  public :: CloseNCFile
  public :: GetDimID
  public :: GetDimLen
  public :: GetVar
  public :: RegisterNCDims
  public :: RegisterVar
  public :: WriteVar
  public :: EndNCDef

  contains

  !=======================================================================================

  logical function CheckFile(filename, fmode)
    !
    ! DESCRIPTION:
    ! Checks to see if a file exists and checks against the mode
    !

    ! ARGUMENTS:
    character(len=*), intent(in) :: filename ! Name of file to open
    character(len=*), intent(in) :: fmode    ! File mode

    ! LOCALS:
    character(len=len(filename)) :: fname       ! Local filename (trimmed)
    integer                      :: ios         ! I/O status
    logical                      :: file_exists ! Does the file exist?

    ! trim filename of whitespace
    fname = trim(adjustl(filename))

    ! Does the file exist?
    inquire(file=fname, exist=file_exists)

    select case (fmode)
    case('read')

      if (.not. file_exists) then
        write(*,'(a,a,a)') "File ", fname(1:len_trim(fname)), " does not exist. Can't read."
        CheckFile = .false.
      else
        CheckFile = .true.
      end if

    case('readwrite')

      CheckFile = .true.

    case('write')
      if (file_exists) then
        write(*, '(a, a, a)') "File ", fname(1:len_trim(fname)), " exists. Cannot open write only."
      else
        CheckFile = .true.
        CheckFile = .false.
      end if
    case default
      write(*,'(a)') "Invalid file mode."
      CheckFile = .false.
    end select

  end function CheckFile

  !=======================================================================================

  subroutine Check(status)
    !
    ! DESCRIPTION:
    ! Checks status of netcdf operations

    ! ARGUMENTS:
    integer, intent(in) :: status ! return status code from a netcdf procedure

    if (status /= nf90_noerr) then
      write(*,*) trim(nf90_strerror(status))
      stop
    end if

  end subroutine Check

  ! =======================================================================================

  subroutine OpenNCFile(nc_file, ncid, fmode)
    !
    ! DESCRIPTION:
    ! Opens a netcdf file

    ! ARGUMENTS:
    character(len=*), intent(in)  :: nc_file ! file name
    integer,          intent(out) :: ncid    ! netcdf file unit number
    character(len=*)              :: fmode   ! file mode

    if (CheckFile(nc_file, fmode)) then
      ! depending on mode
      select case(fmode)
      case ('read')
        call Check(nf90_open(trim(nc_file), NF90_NOCLOBBER, ncid))
      case ('write')
        call Check(nf90_create(trim(nc_file), NF90_CLOBBER, ncid))
      case ('readwrite')
        call Check(nf90_create(trim(nc_file), NF90_CLOBBER, ncid))
      case DEFAULT
        write(*,*) 'Need to specify read, write, or readwrite'
        stop
      end select
    else
      write(*,*) 'Problem reading file'
      stop
    end if

  end subroutine OpenNCFile

  !=======================================================================================

  subroutine CloseNCFile(ncid)
    !
    ! DESCRIPTION:
    ! Closes a netcdf file

    ! ARGUMENTS:
    integer, intent(in) :: ncid ! netcdf file unit number

    call Check(nf90_close(ncid))

  end subroutine CloseNCFile

  !=======================================================================================

  subroutine GetDimID(ncid, var_name, dim_id)
    !
    ! DESCRIPTION:
    ! Gets dimension IDs for a variable ID
    !

    ! ARGUMENTS:
    integer,          intent(in)  :: ncid     ! netcdf file unit number
    character(len=*), intent(in)  :: var_name ! variable name
    integer,          intent(out) :: dim_id   ! dimension ID

    call Check(nf90_inq_dimid(ncid, var_name, dim_id))

  end subroutine GetDimID

  !=======================================================================================

  subroutine GetDimLen(ncid, dim_id, dim_len)
    !
    ! DESCRIPTION:
    ! Gets dimension lengths given a dimension ID
    !

    ! ARGUMENTS:
    integer, intent(in)  :: ncid    ! netcdf file unit number
    integer, intent(in)  :: dim_id  ! dimension ID
    integer, intent(out) :: dim_len ! dimension length

    call Check(nf90_inquire_dimension(ncid, dim_id, len=dim_len))

  end subroutine GetDimLen

  !=======================================================================================

  subroutine GetDims(ncid, varID, dim_lens)
    !
    ! DESCRIPTION:
    ! Get dimensions for a netcdf variable
    !

    ! ARGUMENTS
    integer,              intent(in)  :: ncid        ! netcdf file unit ID
    integer,              intent(in)  :: varID       ! variable ID
    integer, allocatable, intent(out) :: dim_lens(:) ! dimension lengths

    ! LOCALS:
    integer              :: numDims   ! number of dimensions
    integer, allocatable :: dimIDs(:) ! dimension IDs
    integer              :: i         ! looping index

    ! find dimensions of data
    call Check(nf90_inquire_variable(ncid, varID, ndims=numDims))

    ! allocate data to grab dimension information
    allocate(dim_lens(numDims))
    allocate(dimIDs(numDims))

    ! get dimIDs
    call Check(nf90_inquire_variable(ncid, varID, dimids=dimIDs))

    ! grab these dimensions
    do i = 1, numDims
      call Check(nf90_inquire_dimension(ncid, dimIDs(i), len=dim_lens(i)))
    end do

  end subroutine GetDims

  !=====================================================================================

  subroutine GetVarScalarReal(ncid, var_name, data)
    !
    ! DESCRIPTION:
    ! Read in variables for scalar real data
    !

    ! ARGUMENTS:
    integer,          intent(in)  :: ncid     ! netcdf file unit ID
    character(len=*), intent(in)  :: var_name ! variable name
    real(r8),         intent(out) :: data     ! data value

    ! LOCALS:
    integer              :: varID       ! variable ID
    integer, allocatable :: dim_lens(:) ! dimension lengths

    ! find variable ID first
    call Check(nf90_inq_varid(ncid, var_name, varID))

    ! read data
    call Check(nf90_get_var(ncid, varID, data))

  end subroutine GetVarScalarReal

  !=====================================================================================

  subroutine GetVar1DReal(ncid, var_name, data)
    !
    ! DESCRIPTION:
    ! Read in variables for 1D real data
    !

    ! ARGUMENTS:
    integer,                       intent(in)  :: ncid     ! netcdf file unit ID
    character(len=*),              intent(in)  :: var_name ! variable name
    real(r8),         allocatable, intent(out) :: data(:)  ! data values

    ! LOCALS:
    integer              :: varID       ! variable ID
    integer, allocatable :: dim_lens(:) ! dimension lengths

    ! find variable ID first
    call Check(nf90_inq_varid(ncid, var_name, varID))

    ! get dimensions of data
    call GetDims(ncid, varID, dim_lens)

    ! read data
    allocate(data(dim_lens(1)))
    call Check(nf90_get_var(ncid, varID, data))

  end subroutine GetVar1DReal

  !=====================================================================================

  subroutine GetVar1DInt(ncid, var_name, data)
    !
    ! DESCRIPTION:
    ! Read in variables for 1D integer data
    !

    ! ARGUMENTS:
    integer,                       intent(in)  :: ncid     ! netcdf file unit ID
    character(len=*),              intent(in)  :: var_name ! variable name
    integer,          allocatable, intent(out) :: data(:)  ! data values

    ! LOCALS:
    integer              :: varID       ! variable ID
    integer, allocatable :: dim_lens(:) ! dimension lengths

    ! find variable ID first
    call Check(nf90_inq_varid(ncid, var_name, varID))

    ! get dimensions of data
    call GetDims(ncid, varID, dim_lens)

    ! read data
    allocate(data(dim_lens(1)))
    call Check(nf90_get_var(ncid, varID, data))

  end subroutine GetVar1DInt

  !=====================================================================================

  subroutine GetVar2DReal(ncid, var_name, data)
    !
    ! DESCRIPTION:
    ! Read in variables for 2D real data
    !

    ! ARGUMENTS:
    integer,                       intent(in)  :: ncid      ! netcdf file unit ID
    character(len=*),              intent(in)  :: var_name  ! variable name
    real(r8),         allocatable, intent(out) :: data(:,:) ! data values

    ! LOCALS:
    integer              :: varID       ! variable ID
    integer, allocatable :: dim_lens(:) ! dimension lengths

    ! find variable ID first
    call Check(nf90_inq_varid(ncid, var_name, varID))

    ! get dimensions of data
    call GetDims(ncid, varID, dim_lens)

    ! read data
    allocate(data(dim_lens(1), dim_lens(2)))
    call Check(nf90_get_var(ncid, varID, data))

  end subroutine GetVar2DReal

  !=====================================================================================

  subroutine GetVar2DInt(ncid, var_name, data)
    !
    ! DESCRIPTION:
    ! Read in variables for 2D integer data
    !

    ! ARGUMENTS:
    integer,                       intent(in)  :: ncid      ! netcdf file unit ID
    character(len=*),              intent(in)  :: var_name  ! variable name
    integer,          allocatable, intent(out) :: data(:,:) ! data values

    ! LOCALS:
    integer              :: varID       ! variable ID
    integer, allocatable :: dim_lens(:) ! dimension lengths

    ! find variable ID first
    call Check(nf90_inq_varid(ncid, var_name, varID))

    ! get dimensions of data
    call GetDims(ncid, varID, dim_lens)

    ! read data
    allocate(data(dim_lens(1), dim_lens(2)))
    call Check(nf90_get_var(ncid, varID, data))

  end subroutine GetVar2DInt

  !=====================================================================================

  subroutine GetVar3DReal(ncid, var_name, data)
    !
    ! DESCRIPTION:
    ! Read in variables for 3D real data
    !

    ! ARGUMENTS:
    integer,                       intent(in)  :: ncid        ! netcdf file unit ID
    character(len=*),              intent(in)  :: var_name    ! variable name
    real(r8),         allocatable, intent(out) :: data(:,:,:) ! data values

    ! LOCALS:
    integer              :: varID       ! variable ID
    integer, allocatable :: dim_lens(:) ! dimension lengths

    ! find variable ID first
    call Check(nf90_inq_varid(ncid, var_name, varID))

    ! get dimensions of data
    call GetDims(ncid, varID, dim_lens)

    ! read data
    allocate(data(dim_lens(1), dim_lens(2), dim_lens(3)))
    call Check(nf90_get_var(ncid, varID, data))

  end subroutine GetVar3DReal

  !=====================================================================================

  subroutine GetVar3DInt(ncid, var_name, data)
    !
    ! DESCRIPTION:
    ! Read in variables for 3D integer data
    !

    ! ARGUMENTS:
    integer,                       intent(in)  :: ncid        ! netcdf file unit ID
    character(len=*),              intent(in)  :: var_name    ! variable name
    integer,          allocatable, intent(out) :: data(:,:,:) ! data values

    ! LOCALS:
    integer              :: varID       ! variable ID
    integer, allocatable :: dim_lens(:) ! dimension lengths

    ! find variable ID first
    call Check(nf90_inq_varid(ncid, var_name, varID))

    ! get dimensions of data
    call GetDims(ncid, varID, dim_lens)

    ! read data
    allocate(data(dim_lens(1), dim_lens(2), dim_lens(3)))
    call Check(nf90_get_var(ncid, varID, data))

  end subroutine GetVar3DInt

  !=====================================================================================

  subroutine RegisterNCDims(ncid, dim_names, dim_lens, num_dims, dim_IDs)
    !
    ! DESCRIPTION:
    ! Defines variables and dimensions
    !

    ! ARGUMENTS:
    integer,          intent(in)  :: ncid                ! netcdf file id
    character(len=*), intent(in)  :: dim_names(num_dims) ! dimension names
    integer,          intent(in)  :: dim_lens(num_dims)  ! dimension lengths
    integer,          intent(in)  :: num_dims            ! number of dimensions
    integer,          intent(out) :: dim_IDs(num_dims)   ! dimension IDs

    ! LOCALS:
    integer :: i ! looping index

    do i = 1, num_dims
      call Check(nf90_def_dim(ncid, dim_names(i), dim_lens(i), dim_IDs(i)))
    end do

  end subroutine RegisterNCDims

  !=====================================================================================

  subroutine RegisterVar(ncid, var_name, dimID, type, att_names, atts, num_atts, varID)
    !
    ! DESCRIPTION:
    ! Defines variables and dimensions
    !

    ! ARGUMENTS:
    integer,          intent(in)  :: ncid                ! netcdf file id
    character(len=*), intent(in)  :: var_name            ! variable name
    integer,          intent(in)  :: dimID(:)            ! dimension IDs
    integer,          intent(in)  :: type                ! type: int or double
    character(len=*), intent(in)  :: att_names(num_atts) ! attribute names
    character(len=*), intent(in)  :: atts(num_atts)      ! attribute values
    integer,          intent(in)  :: num_atts            ! number of attributes
    integer,          intent(out) :: varID               ! variable ID


    ! LOCALS:
    integer :: i       ! looping index
    integer :: nc_type ! netcdf type

    if (type == type_double) then
      nc_type = NF90_DOUBLE
    else if (type == type_int) then
      nc_type = NF90_INT
    else if (type == type_char) then
      nc_type = NF90_CHAR
    else
      write(*, *) "Must pick correct type"
      stop
    end if

    call Check(nf90_def_var(ncid, var_name, nc_type, dimID, varID))

    do i = 1, num_atts
      call Check(nf90_put_att(ncid, varID, att_names(i), atts(i)))
    end do

  end subroutine RegisterVar

  !  =====================================================================================

  subroutine EndNCDef(ncid)
    !
    ! DESCRIPTION:
    ! End defining of netcdf dimensions and variables
    !

    ! ARGUMENTS:
    integer, intent(in)  :: ncid ! netcdf file id

    call Check(nf90_enddef(ncid))

  end subroutine EndNCDef

  !  =====================================================================================

  subroutine WriteVar1DReal(ncid, varID, data)
    !
    ! DESCRIPTION:
    ! Write 1D real data
    !

    ! ARGUMENTS:
    integer,  intent(in) :: ncid    ! netcdf file id
    integer,  intent(in) :: varID   ! variable ID
    real(r8), intent(in) :: data(:) ! data to write

    call Check(nf90_put_var(ncid, varID, data(:)))

  end subroutine WriteVar1DReal

  !  =====================================================================================

  subroutine WriteVar2DReal(ncid, varID, data)
    !
    ! DESCRIPTION:
    ! Write 2D real data
    !

    ! ARGUMENTS:
    integer,  intent(in) :: ncid      ! netcdf file id
    integer,  intent(in) :: varID     ! variable ID
    real(r8), intent(in) :: data(:,:) ! data to write

    call Check(nf90_put_var(ncid, varID, data(:,:)))

  end subroutine WriteVar2DReal

  !  =====================================================================================

  subroutine WriteVar1DInt(ncid, varID, data)
    !
    ! DESCRIPTION:
    ! Write 1D integer data
    !

    ! ARGUMENTS:
    integer, intent(in) :: ncid    ! netcdf file id
    integer, intent(in) :: varID   ! variable ID
    integer, intent(in) :: data(:) ! data to write

    call Check(nf90_put_var(ncid, varID, data(:)))

  end subroutine WriteVar1DInt

  !  =====================================================================================

  subroutine WriteVar2DInt(ncid, varID, data)
    !
    ! DESCRIPTION:
    ! Write 2D integer data
    !

    ! ARGUMENTS:
    integer, intent(in) :: ncid      ! netcdf file id
    integer, intent(in) :: varID     ! variable ID
    integer, intent(in) :: data(:,:) ! data to write

    call Check(nf90_put_var(ncid, varID, data(:,:)))

  end subroutine WriteVar2DInt

  !  =====================================================================================
  
  subroutine WriteVar1DChar(ncid, varID, data)
    !
    ! DESCRIPTION:
    ! Write 1D character data
    !

    ! ARGUMENTS:
    integer,          intent(in) :: ncid    ! netcdf file id
    integer,          intent(in) :: varID   ! variable ID
    character(len=*), intent(in) :: data(:) ! data to write

    call Check(nf90_put_var(ncid, varID, data(:)))

  end subroutine WriteVar1DChar

  !  =====================================================================================

  subroutine WriteVar2DChar(ncid, varID, data)
    !
    ! DESCRIPTION:
    ! Write 2D character data
    !

    ! ARGUMENTS:
    integer,          intent(in) :: ncid      ! netcdf file id
    integer,          intent(in) :: varID     ! variable ID
    character(len=*), intent(in) :: data(:,:) ! data to write

    call Check(nf90_put_var(ncid, varID, data(:,:)))

  end subroutine WriteVar2DChar

  !  =====================================================================================

end module FatesUnitTestIOMod