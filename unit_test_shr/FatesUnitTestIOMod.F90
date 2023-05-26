module FatesUnitTestIOMod
  use FatesConstantsMod, only : r8 => fates_r8
  use shr_kind_mod,      only : SHR_KIND_CL
  use netcdf
  
  implicit none

  ! LOCALS
  integer, parameter :: BASE_UNIT = 10  ! Base unit for files the first time unit_number is called
  integer, parameter :: MAX_PATH = 256  ! Maximum path length
  integer, parameter :: MAX_CHAR = 80   ! Maximum length for messages
  integer            :: logf            ! Unit number for output log file

  interface GetVar
    module procedure GetVar1DReal
    module procedure GetVar2DReal
    module procedure GetVar3DReal
    module procedure GetVar1DInt
    module procedure GetVar2DInt
    module procedure GetVar3DInt
  end interface
  
  contains

    integer function UnitNumber()
      !
      ! DESCRIPTION:
      ! Generates a unit number to be used in opening files
      ! The first time the function is called, it returns BASE_UNIT
      !
      ! LOCALS:
      integer :: iunit          ! File unit (increments after first call)
      logical :: first = .true. ! First time this has been called?
      save

      if (first) then
        ! Set first to false and iunit to base unit on first call
        iunit = BASE_UNIT
        first = .false.
      else
        ! Otherwise, increment
        iunit = iunit + 1
      endif

      ! Set to output
      UnitNumber = iunit

    end function UnitNumber

    !=====================================================================================

    character(len=9) function FileMode(mode)
      !
      ! DESCRIPTION:
      ! Gets a file mode
      !

      ! ARGUMENTS:
      character(len=*), intent(in), optional :: mode     ! Optional mode ('r', 'w', 'rw')
      
      ! Get mode of open (read, write, or read/write)
      select case(mode)
      case ('r', 'R')
        FileMode = 'read'
      case ('w', 'W')
        FileMode = 'write'
      case ('rw', 'RW', 'wr', 'WR')
        FileMode = 'readwrite'
      case DEFAULT
        FileMode = 'readwrite'
      end select

    end function FileMode

    !=====================================================================================

    logical function CheckFile(filename, fmode)
      !
      ! DESCRIPTION:
      ! Checks to see if a file exists and checks against the mode
      !

      ! ARGUMENTS:
      character(len=*), intent(in) :: filename ! Name of file to open
      character(len=*), intent(in) :: fmode    ! File mode

      ! LOCALS:
      character(len=MAX_PATH)      :: fname       ! Local filename (trimmed)
      integer, dimension(MAX_PATH) :: farray      ! Array of characters of file name
      integer                      :: i           ! looping index
      integer                      :: ios         ! I/O status
      logical                      :: file_exists ! Does the file exist?
      
      ! trim filename of whitespace
      fname = trim(adjustl(filename))

      if (fmode == 'read' .or. fmode == 'readwrite') then
        ! Check for valid name of file
        farray = 0
        do i = 1, len_trim(fname)
          farray(i) = ichar(fname(i:i))
        enddo
        if (any(farray > MAX_PATH)) then
          write(logf,'(A)') "Invalid filename"
          CheckFile = .false.
          return 
        endif
      endif

      ! Does the file exist?
      inquire(file=fname, exist=file_exists)

      ! Open file if conditions are correct
      if (file_exists .and. fmode == 'write') then
        write(logf,'(A,A,A)') "File ", fname(1:len_trim(fname)),                         &
          " exists. Cannot open write only."
        CheckFile = .false.
      else if (.not. file_exists .and. fmode == 'read') then
        write(logf, '(A,A,A)') "File ", fname(1:len_trim(fname)),                        &
          " does not exist. Can't read."
        CheckFile = .false.
      else
        CheckFile = .true.
      endif

    end function CheckFile

    !=====================================================================================

    integer function OpenFile(filename, mode)
      !
      ! DESCRIPTION:
      ! Opens the file filename if it can, returns a unit number for it.
      ! The first time the function is called, it returns BASE_UNIT
      !

      ! ARGUMENTS:
      character(len = *), intent(in)           :: filename ! Name of file to open
      character(len = *), intent(in), optional :: mode     ! Optional mode ('r', 'w', 'rw')

      ! LOCALS:
      character(len=9)        :: fmode      ! file mode
      integer                 :: iunit      ! file unit number
      integer                 :: ios        ! I/O status
      character(len=MAX_PATH) :: fname      ! Local filename (trimmed)
      
      ! get the file mode, defaults to readwrite
      if (present(mode)) then 
        fmode = FileMode(mode)
      else 
        fmode = 'readwrite'
      end if 

      if (CheckFile(filename, fmode)) then 
        
        ! trim filename of whitespace
        fname = trim(adjustl(filename))

        iunit = UnitNumber()
        open(iunit, file=fname, action=fmode, iostat=ios)
        if (ios /= 0) then
          write(logf,'(A,A,A,I6)') "Problem opening",                                    &
            fname(1:len_trim(fname)), " ios: ", ios
          stop 
        endif
      else 
        stop
      end if

      OpenFile = iunit

    end function OpenFile 

    !=====================================================================================

    subroutine Check(status)
      !
      ! DESCRIPTION:
      ! Checks status of netcdf operations
    
      ! ARGUMENTS:
      integer, intent(in) :: status ! return status code from a netcdf procedure
      
      if (status /= nf90_noerr) then 
        write(logf,*) trim(nf90_strerror(status))
        stop 
      end if

    end subroutine Check

    !=====================================================================================

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
          call Check(nf90_create(trim(nc_file), 0, ncid))
        case ('readwrite')
          call Check(nf90_open(trim(nc_file), 0, ncid))
        case DEFAULT
          write(logf,*) 'Need to specify read, write, or readwrite'
          stop
        end select
      else 
        write(logf,*) 'Problem reading file'
        stop
      end if 

    end subroutine OpenNCFile

    !=====================================================================================

    subroutine CloseNCFile(ncid)
      !
      ! DESCRIPTION:
      ! Closes a netcdf file
    
      ! ARGUMENTS:
      integer, intent(in) :: ncid ! netcdf file unit number

      call Check(nf90_close(ncid))

    end subroutine CloseNCFile

    !=====================================================================================
  
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

    subroutine WriteVar2DReal(ncid, dim_names, data)
      !
      ! DESCRIPTION:
      ! Writes variables for 2D real data
      !

      ! ARGUMENTS: 
      integer,               intent(in) :: ncid         ! netcdf file id 
      character(len=*),      intent(in) :: dim_names(2) ! dimension names
      real(r8),              intent(in) :: data(:,:)    ! data to write out

      ! LOCALS:
      integer :: dimID(2) ! dimension IDs
      integer :: i        ! looping index

      do i = 1, 2
        call Check(nf90_def_dim(ncid, dim_names(i), size(data, i), dimID(i)))
        
      end do 

    end subroutine WriteVar2DReal

end module FatesUnitTestIOMod