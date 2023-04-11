module FatesUnitTestIOMod
  use FatesConstantsMod, only : fates_r8 => r8
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

  !:...........................................................................:

  integer function open_file(filename, mode)
    !
    ! DESCRIPTION:
    ! Opens the file filename if it can, returns a unit number for it.
    ! The first time the function is called, it returns BASE_UNIT

    ! ARGUMENTS:
    character(len = *), intent(in)           :: filename ! Name of file to open
    character(len = *), intent(in), optional :: mode     ! Optional mode ('r', 'w', 'rw')

    ! LOCALS:
    character(len = 9)           :: fmode             ! File open mode
    logical                      :: file_exists       ! Does the file exist?
    character(len = MAX_PATH)    :: fname             ! Local filename (trimmed)
    character(len = MAX_CHAR)    :: message           ! Error message
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
            write(*,*) "Invalid filename"
            stop "Stopped"
        endif
    endif

    ! Does the file exist?
    inquire(file = fname, exist = file_exists)

    ! Open file if conditions are correct
    if (file_exists .and. fmode == 'write') then
        write(message, '(A,A,A)') "File ", fname(1:len_trim(fname)),           &
            " exists. Cannot open write only."
        write(*,*) message
        stop "Stopped"
    else if (.not. file_exists .and. fmode == 'read') then
        write(message, '(A,A,A)') "File ", fname(1:len_trim(fname)),           &
            " does not exist. Can't read."
        write(*,*) message
        stop "Stopped"
    else
        iunit = unit_number()
        open(iunit, file = fname, action = fmode, iostat = ios)
        if (ios /= 0) then
            write(message, '(A,A,A,I6)') "Problem opening",                    &
                fname(1:len_trim(fname)), " ios: ", ios
            write(*,*) message
            stop "Stopped"
        endif
    endif

    open_file = iunit

  end function open_file

  !:...........................................................................:

  subroutine write_radiation_data(file, kdir, declin)
    character(len=100), intent(in) :: file              ! output file name
    real(r8),           intent(in) :: kdir(num_pft, 48) ! direct beam extinction coefficient
    real(r8),           intent(in) :: declin(48)

    integer :: ncid
    integer :: pft_dimid, time_dimid
    integer :: kdir_id, declin_id
    integer :: dimids(2)

    ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
    ! overwrite this file, if it already exists.
    call check(nf90_create(trim(file), NF90_CLOBBER, ncid))

    ! define dimensions
    call check(nf90_def_dim(ncid, "pft", num_pft, pft_dimid))
    call check(nf90_def_dim(ncid, "time", 48, time_dimid))

    ! define the variables
    dimids =  (/ pft_dimid, time_dimid/)
    call check(nf90_def_var(ncid, "kdir", NF90_REAL8, dimids, kdir_id))
    call check(nf90_def_var(ncid, "declin", NF90_REAL8, time_dimid, declin_id))

    ! end define mode
    call check(nf90_enddef(ncid))

    ! write to file. 
    call check(nf90_put_var(ncid, kdir_id, kdir))
    call check(nf90_put_var(ncid, declin_id, declin))

    ! close the file.
    call check(nf90_close(ncid))

  end subroutine write_radiation_data

  !:...........................................................................:

  subroutine check(status)
    !
    ! DESCRIPTION:
    ! Checks status of netcdf operations
   
    ! ARGUMENTS:
    integer, intent (in) :: status ! return status code from a netcdf procedure
    
    if (status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if

  end subroutine check

    !:...........................................................................:

    subroutine read_patch_data(file, canopy_area, elai, esai, nrad)
    !
    ! DESCRIPTION:
    ! Reads and return patch data
   
    ! ARGUMENTS:
    character(len=100), intent(in)  :: file                                  ! patch file name
      real(r8),           intent(out) :: canopy_area(num_can,num_pft,nlevleaf) ! canopy area profile
      real(r8),           intent(out) :: elai(num_can,num_pft,nlevleaf)        ! exposed lai profile
      real(r8),           intent(out) :: esai(num_can,num_pft,nlevleaf)        ! exposed sai profile
      real(r8),           intent(out) :: nrad(num_can,num_pft)                 ! number of exposed leaf layers

      ! Data dictionary: declare variable types, definitions, and units
      real(r8) :: nrad_r(num_can,num_pft) ! number of exposed leaf layers
      integer  :: fidA, varID

      ! open file
      call check(nf90_open(trim(file), 0, fidA))

      ! read patch data values
      call check(nf90_inq_varid(fidA, "can_area", varID))
      call check(nf90_get_var(fidA, varID, canopy_area))

      call check(nf90_inq_varid(fidA, "elai", varID))
      call check(nf90_get_var(fidA, varID, elai))

      call check(nf90_inq_varid(fidA, "esai", varID))
      call check(nf90_get_var(fidA, varID, esai))

      call check(nf90_inq_varid(fidA, "nrad", varID))
      call check(nf90_get_var(fidA, varID, nrad))

  end subroutine read_patch_data

  !:...........................................................................:

  subroutine get_parameters(file, rhol, rhos, taul, taus, xl, ci)
      !
      ! Purpose:
      !		Returns parameter values 
      !
      ! Record of revisions:
      !		Date 		   Programmer			  Description of change
      !		========	 =============		=====================
      !	  03/20/23   A. C. Foster     Original code
      !

      ! Data dictionary: declare calling parameter types & definitions
      character(len=100), intent(in)  :: file                  ! parameter file name
      real(r8),           intent(out) :: rhol(num_pft,num_swb) ! leaf reflectance
      real(r8),           intent(out) :: rhos(num_pft,num_swb) ! stem reflectance
      real(r8),           intent(out) :: taul(num_pft,num_swb) ! leaf transmittance
      real(r8),           intent(out) :: taus(num_pft,num_swb) ! stem transmittance
      real(r8),           intent(out) :: xl(num_pft)           ! leaf orientation index
      real(r8),           intent(out) :: ci(num_pft)           ! clumping index

      ! read in all parameters
      call read_param_netcdf(file, rhol(:,1), rhol(:,2), rhos(:,1), rhos(:,2), &
        taul(:,1), taul(:,2), taus(:,1), taus(:,2), xl, ci)

  end subroutine get_parameters

  !:...........................................................................:

  subroutine read_param_netcdf(file, rholvis, rholnir, rhosvis, rhosnir,       &
      taulvis, taulnir, tausvis, tausnir, xl, ci)
      !
      ! Purpose:
      !		Reads a FATES netcdf parameter file
      !
      ! Record of revisions:
      !		Date 		   Programmer			  Description of change
      !		========	 =============		=====================
      !	  03/20/23   A. C. Foster     Original code
      !

    ! Data dictionary: declare calling parameter types & definitions
    character(len=100), intent(in)  :: file             ! parameter file name 
    real(r8),           intent(out) :: rholvis(num_pft) ! vis leaf reflectance
    real(r8),           intent(out) :: rholnir(num_pft) ! nir leaf reflectance
    real(r8),           intent(out) :: rhosvis(num_pft) ! vis stem reflectance
    real(r8),           intent(out) :: rhosnir(num_pft) ! nir stem reflectance
    real(r8),           intent(out) :: taulvis(num_pft) ! vis leaf transmittance
    real(r8),           intent(out) :: taulnir(num_pft) ! nir leaf transmittance
    real(r8),           intent(out) :: tausvis(num_pft) ! vis stem transmittance
    real(r8),           intent(out) :: tausnir(num_pft) ! nir stem transmittance
    real(r8),           intent(out) :: xl(num_pft)      ! leaf orientation index
    real(r8),           intent(out) :: ci(num_pft)      ! clumping index

    ! Data dictionary: declare variable types, definitions, and units
    integer :: fidA, paramID
      
    ! open file
    call check(nf90_open(trim(file), 0, fidA))
  
    ! read parameter values
    call check(nf90_inq_varid(fidA, "fates_rad_leaf_rhovis", paramID))
    call check(nf90_get_var(fidA, paramID, rholvis))

    call check(nf90_inq_varid(fidA, "fates_rad_leaf_rhonir", paramID))
    call check(nf90_get_var(fidA, paramID, rholnir))

    call check(nf90_inq_varid(fidA, "fates_rad_stem_rhovis", paramID))
    call check(nf90_get_var(fidA, paramID, rhosvis))

    call check(nf90_inq_varid(fidA, "fates_rad_stem_rhonir", paramID))
    call check(nf90_get_var(fidA, paramID, rhosnir))

    call check(nf90_inq_varid(fidA, "fates_rad_leaf_tauvis", paramID))
    call check(nf90_get_var(fidA, paramID, taulvis))

    call check(nf90_inq_varid(fidA, "fates_rad_leaf_taunir", paramID))
    call check(nf90_get_var(fidA, paramID, taulnir))

    call check(nf90_inq_varid(fidA, "fates_rad_stem_tauvis", paramID))
    call check(nf90_get_var(fidA, paramID, tausvis))

    call check(nf90_inq_varid(fidA, "fates_rad_stem_taunir", paramID))
    call check(nf90_get_var(fidA, paramID, tausnir))

    call check(nf90_inq_varid(fidA, "fates_rad_leaf_xl", paramID))
    call check(nf90_get_var(fidA, paramID, xl))

    call check(nf90_inq_varid(fidA, "fates_rad_leaf_clumping_index", paramID))
    call check(nf90_get_var(fidA, paramID, ci))

    ! close file
    call check(nf90_close(fidA))

  end subroutine read_param_netcdf

  !:...........................................................................:

end module FatesUnitTestIOMod