module FatesUnitTestIOMod
  use FatesConstantsMod,          only : r8 => fates_r8
  use FatesParametersInterface,   only : fates_parameters_type
  use FatesParametersInterface,   only : param_string_length, max_used_dimensions
  use EDParamsMod,                only : FatesRegisterParams, FatesReceiveParams
  use SFParamsMod,                only : SpitFireRegisterParams, SpitFireReceiveParams
  use PRTInitParamsFATESMod,      only : PRTRegisterParams, PRTReceiveParams
  use shr_kind_mod,               only : SHR_KIND_CL
  use FatesSynchronizedParamsMod, only : FatesSynchronizedParamsInst
  use netcdf
  
  implicit none

  ! LOCALS
  integer, parameter :: BASE_UNIT = 10  ! Base unit for files the first time unit_number is called
  integer, parameter :: MAX_PATH = 256  ! Maximum path length
  integer, parameter :: MAX_CHAR = 80   ! Maximum length for messages
  integer            :: logf            ! Unit number for output log file
  
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

    subroutine read_parameters(fates_paramfile)
      !
      ! DESCRIPTION:
      !		Reads in parameters from the FATES parameter file
      !
    
      ! ARGUMENTS:
      character(len=SHR_KIND_CL), intent(in) :: fates_paramfile ! parameter file name
   
      ! LOCALS:
      class(fates_parameters_type), allocatable :: fates_params

      allocate(fates_params)
      call fates_params%Init()
      call FatesRegisterParams(fates_params)
      call SpitFireRegisterParams(fates_params)
      call PRTRegisterParams(fates_params)
      call FatesSynchronizedParamsInst%RegisterParams(fates_params)

    end subroutine read_parameters

    !:.........................................................................:

    subroutine read_netcdf_params(filename, fates_params)
      !
      ! DESCRIPTION:
      !   Calls actual netcdf library methods for reading FATES parameters from netcdf
      !

      ! ARGUMENTS:
      character(len=*),             intent(in)    :: filename     ! full path of parameter file
      class(fates_parameters_type), intent(inout) :: fates_params ! fates parameters type
      
      ! LOCALS:
      logical :: file_exists  ! does the file exist?
      integer :: ncid         ! netcdf file unit number
      integer :: max_dim_size !

      ! check if file is on disk
      inquire(file=trim(adjustl(filename)), exist=file_exists)
      if (.not. file_exists) then 
        write(logf,'(a)') 'File ', filename, ' does not exist.'
        stop "Stopping"
      end if 
      
      ! open the file
      call check(nf90_open(trim(adjustl(filename)), 0, ncid))

      ! get and set the correct dimensions for the parameters
      call set_param_dimensions(ncid, fates_params)

      ! max_dim_size = fates_params%GetMaxDimensionSize()
      ! num_params = fates_params%num_params()
      ! do i = 1, num_params
      !    call fates_params%GetMetaData(i, name, dimension_shape, dimension_sizes, dimension_names, is_host_param)
      !    if (is_host_file .eqv. is_host_param) then
      !       select case(dimension_shape)
      !       case(dimension_shape_scalar)
      !          size_dim_1 = 1
      !          size_dim_2 = 1
      !       case(dimension_shape_1d)
      !          size_dim_1 = dimension_sizes(1)
      !          size_dim_2 = 1
      !       case(dimension_shape_2d)
      !          size_dim_1 = dimension_sizes(1)
      !          size_dim_2 = dimension_sizes(2)
      !       case default
      !          write(fates_log(),*) 'dimension shape:',dimension_shape
      !          call endrun(msg='unsupported number of dimensions reading parameters.')
   
      !       end select
      !       if (masterproc) then
      !          write(fates_log(), *) 'clmfates_interfaceMod.F90:: reading '//trim(name)
      !       end if
      !       call readNcdio(ncid, name, dimension_shape, dimension_names, subname, data(1:size_dim_1, 1:size_dim_2))
      !       call fates_params%SetData(i, data(1:size_dim_1, 1:size_dim_2))
      !    end if
      ! end do
      ! deallocate(data)
      ! call ncd_pio_closefile(ncid)

    end subroutine read_netcdf_params

    !:.........................................................................:

    subroutine set_param_dimensions(ncid, fates_params)
      !
      ! DESCRIPTION:
      !   get the list of dimensions used by the FATES parameters
      !

      ! ARGUMENTS:
      integer,                      intent(inout) :: ncid         ! netcdf unit number
      class(fates_parameters_type), intent(inout) :: fates_params ! fates parameters

      ! LOCALS:
      integer                            :: num_dimensions
      character(len=param_string_length) :: dimension_names(max_used_dimensions)
      integer                            :: dimension_sizes(max_used_dimensions)
      integer                            :: d
      integer                            :: dim_id
    
      dimension_sizes(:) = 0

      call fates_params%GetUsedDimensions(.false., num_dimensions,             &
        dimension_names)

      do d = 1, num_dimensions
        call check(nf90_inq_dimid(ncid, dimension_names(d), dim_id))
        call check(nf90_inquire_dimension(ncid, dim_id, len=dimension_sizes(d)))
      end do

      call fates_params%SetDimensionSizes(.false., num_dimensions,             &
        dimension_names, dimension_sizes)

    end subroutine set_param_dimensions

    !:.........................................................................:

    ! subroutine read_patch_data(file, canopy_area, elai, esai, nrad)
    !   !
    !   ! DESCRIPTION:
    !   ! Reads and return patch data
    
    !   ! ARGUMENTS:
    !   character(len=MAX_PATH), intent(in)  :: file                                  ! patch file name
    !   real(r8),                intent(out) :: canopy_area(num_can,num_pft,nlevleaf) ! canopy area profile
    !   real(r8),                intent(out) :: elai(num_can,num_pft,nlevleaf)        ! exposed lai profile
    !   real(r8),                intent(out) :: esai(num_can,num_pft,nlevleaf)        ! exposed sai profile
    !   real(r8),                intent(out) :: nrad(num_can,num_pft)                 ! number of exposed leaf layers

    !   ! LOCALS:
    !   real(r8) :: nrad_r(num_can,num_pft) ! number of exposed leaf layers
    !   integer  :: fidA, varID

    !   ! open file
    !   call check(nf90_open(trim(file), 0, fidA))

    !   ! read patch data values
    !   call check(nf90_inq_varid(fidA, "can_area", varID))
    !   call check(nf90_get_var(fidA, varID, canopy_area))

    !   call check(nf90_inq_varid(fidA, "elai", varID))
    !   call check(nf90_get_var(fidA, varID, elai))

    !   call check(nf90_inq_varid(fidA, "esai", varID))
    !   call check(nf90_get_var(fidA, varID, esai))

    !   call check(nf90_inq_varid(fidA, "nrad", varID))
    !   call check(nf90_get_var(fidA, varID, nrad))

    ! end subroutine read_patch_data

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