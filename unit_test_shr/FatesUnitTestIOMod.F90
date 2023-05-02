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

  interface get_var
    module procedure get_var1D_real
    module procedure get_var2D_real
    module procedure get_var3D_real
    module procedure get_var1D_int
    module procedure get_var2D_int
    module procedure get_var3D_int
  end interface
  
  contains

    integer function unit_number()
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
      unit_number = iunit

    end function unit_number

    !:.........................................................................:

    integer function open_file(filename, mode)
      !
      ! DESCRIPTION:
      ! Opens the file filename if it can, returns a unit number for it.
      ! The first time the function is called, it returns BASE_UNIT
      !

      ! ARGUMENTS:
      character(len = *), intent(in)           :: filename ! Name of file to open
      character(len = *), intent(in), optional :: mode     ! Optional mode ('r', 'w', 'rw')

      ! LOCALS:
      character(len = 9)           :: fmode             ! File open mode
      logical                      :: file_exists       ! Does the file exist?
      character(len = MAX_PATH)    :: fname             ! Local filename (trimmed)
      integer                      :: i                 ! Looping index
      integer                      :: ios               ! I/O status
      integer                      :: iunit = BASE_UNIT ! File unit number
      integer, dimension(MAX_PATH) :: farray            ! Array of characters of file name

      ! Get mode of open (read, write, or read/write)
      ! Defaults to read/write
      if (present(mode)) then
        select case(mode)
        case ('r', 'R')
          fmode = 'read'
        case ('w', 'W')
          fmode = 'write'
        case ('rw', 'RW', 'wr', 'WR')
          fmode = 'readwrite'
        case DEFAULT
          fmode = 'readwrite'
        end select
      else
        fmode = 'readwrite'
      endif

      ! trim filename of whitespace
      fname = trim(adjustl(filename))

      if (fmode == 'read' .or. fmode == 'readwrite') then
        ! Check for valid name of file
        farray = 0
        do i = 1, len_trim(fname)
          farray(i) = ichar(fname(i:i))
        enddo
        if (any(farray > 126)) then
          write(logf,'(A)') "Invalid filename"
          stop 
        endif
      endif

      ! Does the file exist?
      inquire(file = fname, exist = file_exists)

      ! Open file if conditions are correct
      if (file_exists .and. fmode == 'write') then
        write(logf,'(A,A,A)') "File ", fname(1:len_trim(fname)),               &
          " exists. Cannot open write only."
        stop 
      else if (.not. file_exists .and. fmode == 'read') then
        write(logf, '(A,A,A)') "File ", fname(1:len_trim(fname)),              &
          " does not exist. Can't read."
        stop 
      else
        iunit = unit_number()
        open(iunit, file=fname, action=fmode, iostat=ios)
        if (ios /= 0) then
          write(logf,'(A,A,A,I6)') "Problem opening",                          &
            fname(1:len_trim(fname)), " ios: ", ios
          stop 
        endif
      endif

      open_file = iunit

    end function open_file 

    !:.........................................................................:

    subroutine check(status)
      !
      ! DESCRIPTION:
      ! Checks status of netcdf operations
    
      ! ARGUMENTS:
      integer, intent (in) :: status ! return status code from a netcdf procedure
      
      if (status /= nf90_noerr) then 
        write(logf,*) trim(nf90_strerror(status))
        stop 
      end if

    end subroutine check

    !:.........................................................................:
    
    subroutine read_patch_data(file, num_pft, canopy_area, elai, esai, nrad)
      !
      ! DESCRIPTION:
      ! Reads and return patch data
      !
    
      ! ARGUMENTS:
      character(len=MAX_PATH), intent(in)  :: file               ! patch file name
      integer,                 intent(in)  :: num_pft            ! number of PFTs
      real(r8), allocatable,   intent(out) :: canopy_area(:,:,:) ! canopy area profile
      real(r8), allocatable,   intent(out) :: elai(:,:,:)        ! exposed lai profile
      real(r8), allocatable,   intent(out) :: esai(:,:,:)        ! exposed sai profile
      integer,  allocatable,   intent(out) :: nrad(:,:)          ! number of exposed leaf layers

      ! LOCALS:
      integer :: ncid ! netcdf file unit number

      ! open file
      call check(nf90_open(trim(file), 0, ncid))
      
      ! read in data
      call get_var(ncid, 'can_area', canopy_area)
      call get_var(ncid, 'elai', elai)
      call get_var(ncid, 'esai', esai)
      call get_var(ncid, 'nrad', nrad)

      ! close file
      call check(nf90_close(ncid))

    end subroutine read_patch_data

    !:.........................................................................:

    subroutine get_dims(ncid, varID, dim_lens)
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
      call check(nf90_inquire_variable(ncid, varID, ndims=numDims))

      ! allocate data to grab dimension information
      allocate(dim_lens(numDims))
      allocate(dimIDs(numDims))

      ! get dimIDs
      call check(nf90_inquire_variable(ncid, varID, dimids=dimIDs))

      ! grab these dimensions
      do i = 1, numDims
        call check(nf90_inquire_dimension(ncid, dimIDs(i), len=dim_lens(i)))
      end do

    end subroutine get_dims

    !:.........................................................................:

    subroutine get_var1D_real(ncid, var_name, data)
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
      call check(nf90_inq_varid(ncid, var_name, varID))

      ! get dimensions of data
      call get_dims(ncid, varID, dim_lens)

      ! read data
      allocate(data(dim_lens(1)))
      call check(nf90_get_var(ncid, varID, data))
  
    end subroutine get_var1D_real

    !:.........................................................................:

    subroutine get_var1D_int(ncid, var_name, data)
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
      call check(nf90_inq_varid(ncid, var_name, varID))

      ! get dimensions of data
      call get_dims(ncid, varID, dim_lens)

      ! read data
      allocate(data(dim_lens(1)))
      call check(nf90_get_var(ncid, varID, data))
  
    end subroutine get_var1D_int

    !:.........................................................................:

    subroutine get_var2D_real(ncid, var_name, data)
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
      call check(nf90_inq_varid(ncid, var_name, varID))

      ! get dimensions of data
      call get_dims(ncid, varID, dim_lens)

      ! read data
      allocate(data(dim_lens(1), dim_lens(2)))
      call check(nf90_get_var(ncid, varID, data))
  
    end subroutine get_var2D_real

    !:.........................................................................:

    subroutine get_var2D_int(ncid, var_name, data)
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
      call check(nf90_inq_varid(ncid, var_name, varID))

      ! get dimensions of data
      call get_dims(ncid, varID, dim_lens)

      ! read data
      allocate(data(dim_lens(1), dim_lens(2)))
      call check(nf90_get_var(ncid, varID, data))
  
    end subroutine get_var2D_int

    !:.........................................................................:

    subroutine get_var3D_real(ncid, var_name, data)
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
      call check(nf90_inq_varid(ncid, var_name, varID))

      ! get dimensions of data
      call get_dims(ncid, varID, dim_lens)

      ! read data
      allocate(data(dim_lens(1), dim_lens(2), dim_lens(3)))
      call check(nf90_get_var(ncid, varID, data))
  
    end subroutine get_var3D_real

    !:.........................................................................:

    subroutine get_var3D_int(ncid, var_name, data)
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
      call check(nf90_inq_varid(ncid, var_name, varID))

      ! get dimensions of data
      call get_dims(ncid, varID, dim_lens)

      ! read data
      allocate(data(dim_lens(1), dim_lens(2), dim_lens(3)))
      call check(nf90_get_var(ncid, varID, data))
  
    end subroutine get_var3D_int

    !:.........................................................................:

  !   subroutine write_radiation_data(file, kdir, declin)
  !     !
  !     ! DESCRIPTION:
  !     ! Opens the file filename if it can, returns a unit number for it.
  !     ! The first time the function is called, it returns BASE_UNIT
  !     !
      
  !     ! ARGUMENTS:
  !     character(len=MAX_PATH), intent(in) :: file             ! output file name
  !     real(r8),                intent(in) :: kdir(num_pft,48) ! direct beam extinction coefficient
  !     real(r8),                intent(in) :: declin(48)

  !     ! LOCALS:
  !     integer :: ncid
  !     integer :: pft_dimid, time_dimid
  !     integer :: kdir_id, declin_id
  !     integer :: dimids(2)

  !     ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
  !     ! overwrite this file, if it already exists.
  !     call check(nf90_create(trim(file), NF90_CLOBBER, ncid))

  !     ! define dimensions
  !     call check(nf90_def_dim(ncid, "pft", num_pft, pft_dimid))
  !     call check(nf90_def_dim(ncid, "time", 48, time_dimid))

  !     ! define the variables
  !     dimids =  (/ pft_dimid, time_dimid/)
  !     call check(nf90_def_var(ncid, "kdir", NF90_REAL8, dimids, kdir_id))
  !     call check(nf90_def_var(ncid, "declin", NF90_REAL8, time_dimid, declin_id))

  !     ! end define mode
  !     call check(nf90_enddef(ncid))

  !     ! write to file. 
  !     call check(nf90_put_var(ncid, kdir_id, kdir))
  !     call check(nf90_put_var(ncid, declin_id, declin))

  !     ! close the file.
  !     call check(nf90_close(ncid))

  !   end subroutine write_radiation_data

end module FatesUnitTestIOMod