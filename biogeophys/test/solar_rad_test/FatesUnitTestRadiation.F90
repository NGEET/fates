program FatesUnitTestRadiation
  !
  ! DESCRIPTION:
  !		Test the FATES radiation schemes
  !

  use FatesUnitTestIOMod,      only : logf, MAX_PATH, open_file
  use FatesUnitTestOrbitalMod, only : SHR_KIND_R8
  use FatesUnitTestOrbitalMod, only : get_orbital_vals, shr_orb_cosz
  use FatesUnitTestOrbitalMod, only : shr_orb_decl
  use FatesConstantsMod,       only : r8 => fates_r8
  use EDPftvarcon,             only : EDPftvarcon_type
  !use EDParamsMod,             only : nlevleaf, nclmax

  implicit none
    
  ! LOCALS:
  character(len=MAX_PATH) :: param_file ! parameter file name
  character(len=MAX_PATH) :: patch_file ! patch data file name
  character(len=MAX_PATH) :: out_file   ! output file name
  integer                 :: year, jday ! year and day of year to simulate
  real(r8)                :: lat, lon   ! latitude/longitude to simulate [degrees]
  real(r8)                :: fcansno    ! fraction of canopy covered by snow [0-1]
  real(r8), allocatable   :: rhol(:,:)  ! leaf reflectance [0-1]
  real(r8), allocatable   :: rhos(:,:)  ! stem reflectance [0-1]
  real(r8), allocatable   :: taul(:,:)  ! leaf transmittance [0-1]
  real(r8), allocatable   :: taus(:,:)  ! stem transmittance [0-1]
  real(r8), allocatable   :: xl(:)      ! leaf orientation index
  real(r8), allocatable   :: ci(:)      ! clumping index

  ! PARAMETERS
  integer, parameter :: numSWb = 2 ! number of shortwave bands to simulate

  interface

    subroutine read_radiation_namelist(year, jday, lat, lon, fcansno,          &
        param_file, patch_file, out_file)

      use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file
      use FatesConstantsMod,  only : r8 => fates_r8

      implicit none

      character(len=MAX_PATH), intent(out) :: param_file, patch_file, out_file 
      integer,                 intent(out) :: year, jday 
      real(r8),                intent(out) :: lat, lon
      real(r8),                intent(out) :: fcansno

    end subroutine read_radiation_namelist

    subroutine read_radiation_params(file, numSWb, rhol, rhos, taul, taus, xl, &
      ci)

      use FatesUnitTestIOMod, only : MAX_PATH, check, read_in_parameter
      use FatesConstantsMod,  only : r8 => fates_r8
      use netcdf

      implicit none
    
      character(len=MAX_PATH), intent(in)  :: file          
      integer,                 intent(in)  :: numSWb    
      real(r8), allocatable,   intent(out) :: rhol(:,:) 
      real(r8), allocatable,   intent(out) :: rhos(:,:) 
      real(r8), allocatable,   intent(out) :: taul(:,:) 
      real(r8), allocatable,   intent(out) :: taus(:,:) 
      real(r8), allocatable,   intent(out) :: xl(:)     
      real(r8), allocatable,   intent(out) :: ci(:)     

    end subroutine read_radiation_params

  end interface

  call read_radiation_namelist(year, jday, lat, lon, fcansno,                  &
    param_file, patch_file, out_file)
  
  ! open log file
  logf = open_file("log.txt")

  ! read in FATES parameter file
  call read_radiation_params(param_file, numSWb, rhol, rhos, taul, taus, xl, ci)

  ! read in patch data
  call read_patch_data(patch_file, canopy_area_profile, elai_profile,          &
    esai_profile, nrad_r)

  ! call get_orbital_vals(year, logf, eccen, mvelp, obliqr, lambm0, mvelpp)

  ! cosz(:) = 0.0
  ! k_dir(:,:) = 0.0
  
  ! ! for each half-hourly time step in the day
  ! do i = 0, 47
  !     calday = jday + i*0.02083333_SHR_KIND_R8
  !     call shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, declin, eccf)
  !     cosz(i) = shr_orb_cosz(calday, lat, lon, declin)

  !     if (cosz(i) > 0.0_r8) then 
  !         ! call norman radiation scheme
  !         call PatchNormanRadiation(rhol, rhos, taul, taus, xl, clumping_index, &
  !             canopy_area_profile, elai_profile, esai_profile, fcansno,         &
  !             cosz(i), nrad, 2, k_dir(:,i))
  !     end if
  ! end do 

  ! call write_radiation_data(out_file, k_dir, cosz)

  close(logf)

end program FatesUnitTestRadiation

!:.............................................................................:

subroutine read_radiation_namelist(year, jday, lat, lon, fcansno, param_file,  &
    patch_file, out_file)
  !
  ! DESCRIPTION:
  !		read in the namelist associated with the radiation unit tests and 
  !   initialize values
  !

  use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file
  use FatesConstantsMod,  only : r8 => fates_r8
  
  implicit none

  ! ARGUMENTS:
  character(len=MAX_PATH), intent(out) :: param_file ! parameter file name
  character(len=MAX_PATH), intent(out) :: patch_file ! patch data file name
  character(len=MAX_PATH), intent(out) :: out_file   ! output file name
  integer,                 intent(out) :: year, jday ! year and day of year
  real(r8),                intent(out) :: lat, lon   ! latitude and longitude [degrees]
  real(r8),                intent(out) :: fcansno    ! fraction of canopy covered by snow [0-1]
  
  ! LOCALS:
  character(len=MAX_PATH) :: rad_nl = 'radiation_nl' ! radiation namelist name
  character(len=MAX_CHAR) :: message                 ! Error message
  character(len=MAX_CHAR) :: msg                     ! I/O Error message
  integer                 :: rad_nl_file             ! unit number for namelist
  integer                 :: ios                     ! I/O status
  

  ! Namelist of radiation parameters
  namelist /radiation/ year, jday, lat, lon, fcansno, param_file, patch_file,  &
    out_file

  ! Now read parameters namelist
  rad_nl_file = open_file(trim(rad_nl), 'r')
  read(rad_nl_file, radiation, iostat=ios, iomsg=msg)
  
  if (ios /= 0) then
    ! Problem reading file - tell user.
    write(message, '(A, I6, A, A)') "Error reading radiation namelist file",   &
      ios, "IOMSG: ", msg
    write(*,*) message
    stop "Stopped"
  end if

  close(rad_nl_file)

end subroutine read_radiation_namelist

!:.............................................................................:

subroutine read_radiation_params(file, numSWb, rhol, rhos, taul, taus, xl, ci)
  !
  ! DESCRIPTION: 
  ! read in the parameters we need for this test
  !

  use FatesUnitTestIOMod, only : MAX_PATH, check, read_in_parameter
  use FatesConstantsMod,  only : r8 => fates_r8
  use netcdf

  implicit none

  ! ARGUMENTS:
  character(len=MAX_PATH), intent(in)  :: file           ! parameter file name
  integer,                 intent(in)  :: numSWb         ! number of shortwave bands to simulate
  real(r8), allocatable,   intent(out) :: rhol(:,:)      ! leaf reflectance [0-1]
  real(r8), allocatable,   intent(out) :: rhos(:,:)      ! stem reflectance [0-1]
  real(r8), allocatable,   intent(out) :: taul(:,:)      ! leaf transmittance [0-1]
  real(r8), allocatable,   intent(out) :: taus(:,:)      ! stem transmittance [0-1]
  real(r8), allocatable,   intent(out) :: xl(:)          ! leaf orientation index
  real(r8), allocatable,   intent(out) :: ci(:)          ! clumping index

  ! LOCALS:
  integer               :: funit      ! file unit number
  real(r8), allocatable :: rholvis(:) ! leaf visible reflectance [0-1]
  real(r8), allocatable :: rholnir(:) ! leaf NIR reflectance [0-1]
  real(r8), allocatable :: rhosvis(:) ! stem visible reflectance [0-1]
  real(r8), allocatable :: rhosnir(:) ! stem NIR reflectance [0-1]
  real(r8), allocatable :: taulvis(:) ! leaf visible transmittance [0-1]
  real(r8), allocatable :: taulnir(:) ! leaf NIR transmittance [0-1]
  real(r8), allocatable :: tausvis(:) ! stem visible transmittance [0-1]
  real(r8), allocatable :: tausnir(:) ! stem NIR transmittance [0-1]

  ! open file
  call check(nf90_open(trim(file), 0, funit))

  ! read in parameters
  call read_in_parameter(funit, 'fates_rad_leaf_rhovis', 'fates_pft', rholvis)
  call read_in_parameter(funit, 'fates_rad_leaf_rhonir', 'fates_pft', rholnir)
  call read_in_parameter(funit, 'fates_rad_stem_rhovis', 'fates_pft', rhosvis)
  call read_in_parameter(funit, 'fates_rad_stem_rhonir', 'fates_pft', rhosnir)
  call read_in_parameter(funit, 'fates_rad_leaf_tauvis', 'fates_pft', taulvis)
  call read_in_parameter(funit, 'fates_rad_leaf_taunir', 'fates_pft', taulnir)
  call read_in_parameter(funit, 'fates_rad_stem_tauvis', 'fates_pft', tausvis)
  call read_in_parameter(funit, 'fates_rad_stem_taunir', 'fates_pft', tausnir)
  call read_in_parameter(funit, 'fates_rad_leaf_xl', 'fates_pft', xl)
  call read_in_parameter(funit, 'fates_rad_leaf_clumping_index', 'fates_pft',  &
    ci)

  ! allocate the arrays correctly
  allocate(rhol(size(rholvis, 1), numSWb))
  allocate(rhos(size(rhosvis, 1), numSWb))
  allocate(taul(size(taulvis, 1), numSWb))
  allocate(taus(size(tausvis, 1), numSWb))

  ! put arrays together
  rhol(:,1) = rholvis
  rhol(:,2) = rholnir 
  rhos(:,1) = rhosvis
  rhos(:,2) = rhosnir 
  taul(:,1) = taulvis 
  taul(:,2) = taulnir 
  taus(:,1) = tausvis
  taus(:,2) = tausnir

  ! close file
  call check(nf90_close(funit))

end subroutine read_radiation_params



