program FatesUnitTestRadiation
  !
  ! DESCRIPTION:
  !		Test the FATES radiation schemes
  !

  use FatesUnitTestIOMod,      only : logf, MAX_PATH
  use FatesUnitTestOrbitalMod, only : SHR_KIND_R8
  use FatesUnitTestOrbitalMod, only : get_orbital_vals, shr_orb_cosz
  use FatesUnitTestOrbitalMod, only : shr_orb_decl
  !use EDParamsMod,             only : nlevleaf, nclmax

  implicit none
    
  ! LOCALS:
  character(len=256) :: param_file, patch_file, out_file  ! file names
  integer                 :: numpft, numSWb
  ! real(r8) :: rhol(num_pft,num_swb)   ! leaf reflectance (0-1)
  ! real(r8) :: rhos(num_pft,num_swb)   ! stem reflectance (0-1)
  ! real(r8) :: taul(num_pft,num_swb)   ! leaf transmittance (0-1)
  ! real(r8) :: taus(num_pft,num_swb)   ! stem transmittance (0-1)
  ! real(r8) :: xl(num_pft)             ! leaf/stem orientation index (-1-1)
  ! real(r8) :: clumping_index(num_pft) ! clumping index (0-1)
  ! real(r8) :: fcansno                 ! fraction of canopy covered by snow
  ! real(r8) :: canopy_area_profile(num_can,num_pft,nlevleaf) ! fraction of crown area per canopy area in each layer
  ! real(r8) :: elai_profile(num_can,num_pft,nlevleaf)       ! exposed leaf area in each canopy layer, pft, and leaf layer
  ! real(r8) :: esai_profile(num_can,num_pft,nlevleaf)       ! exposed stem area in each canopy layer, pft, and leaf layer
  ! real(r8) :: nrad_r(num_can,num_pft)   ! number of exposed leaf layers for each canopy layer and pft
  ! integer  :: nrad(num_can,num_pft)   ! number of exposed leaf layers for each canopy layer and pft
  ! real(r8) :: k_dir(num_pft,48)          ! direct beam extinction coefficient
  ! real(SHR_KIND_R8) :: eccen    ! orbital eccentricity
  ! real(SHR_KIND_R8) :: mvelp    ! moving vernal equinox long
  ! real(SHR_KIND_R8) :: obliqr   ! Earths obliquity in rad
  ! real(SHR_KIND_R8) :: lambm0   ! Mean long of perihelion at vernal equinox (radians)
  ! real(SHR_KIND_R8) :: mvelpp   ! moving vernal equinox long of perihelion plus pi (rad)
  ! real(SHR_KIND_R8) :: eccf     ! Earth-sun distance factor (ie. (1/r)**2)
  ! real(SHR_KIND_R8) :: calday   ! calendar day (including fraction)
  ! real(SHR_KIND_R8) :: declin   ! solar declination (radians)
  ! real(SHR_KIND_R8) :: cosz(48) ! cosine of solar zenith angle (radians)
  ! real(r8)          :: lat, lon
  ! integer           :: jday     ! julian day
  ! integer           :: year     ! year
  ! integer           :: i        ! looping index

  interface

    subroutine read_radiation_namelist(numpft, numSWb, param_file, patch_file, &
      out_file)

      use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file

      implicit none

      character(len=MAX_PATH), intent(out) :: param_file, patch_file, out_file  ! file names
      integer,                 intent(out) :: numpft, numSWb

    end subroutine read_radiation_namelist

  end interface

  print *, 'Hello, unit test'
    
  ! file_names
  ! logf = open_file("log.txt")

  ! ! set julian day, year, and lat/lon here
  ! jday = 165
  ! year = 2000
  ! lat = 45.76_r8*pi_const/180.0_r8
  ! lon = 237.67_r8*pi_const/180.0_r8
  ! fcansno = 0.0_r8

  ! ! get patch and parameter values, as well as orbital parameters
  ! call get_parameters(file_in, rhol, rhos, taul, taus, xl, clumping_index)
  ! call read_patch_data(patch_file, canopy_area_profile, elai_profile,        &
  !     esai_profile, nrad_r)
  ! nrad = int(nrad_r)
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

end program FatesUnitTestRadiation

!:.............................................................................:

subroutine read_radiation_namelist(numpft, numSWb, param_file, patch_file,     &
  out_file)

  use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file
  
  implicit none

  !
  ! DESCRIPTION:
  !		read in the namelist associated with the radiation unit tests and 
  !   initialize values
  !

  ! ARGUMENTS:
  character(len=MAX_PATH), intent(out) :: param_file, patch_file, out_file  ! file names
  integer,                 intent(out) :: numpft, numSWb

  ! LOCALS:
  character(len=MAX_PATH) :: rad_nl = 'radiation_nl' ! radiation namelist name
  character(len=MAX_PATH) :: message                 ! Error message
  character(len=MAX_PATH) :: msg                     ! I/O Error message
  integer                 :: rad_nl_file             ! unit number for namelist
  integer                 :: ios                     ! I/O status
  

  ! Namelist of radiation parameters
  namelist /radiation/ numpft, numSWb, param_file, patch_file, out_file


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