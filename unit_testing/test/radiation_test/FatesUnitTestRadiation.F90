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
  character(len=MAX_PATH) :: param_file, patch_file, out_file  ! file names
  integer                 :: numSWb                            ! number of shortwave bands to simulate
  integer                 :: year, jday                        ! year and day of year to simulate
  real(r8)                :: lat, lon                          ! latitude/longitude to simulate [degrees]
  real(r8)                :: fcansno                           ! fraction of canopy covered by snow [0-1]

  interface

    subroutine read_radiation_namelist(numSWb, year, jday, lat, lon, fcansno,  &
        param_file, patch_file, out_file)

      use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file

      implicit none

      character(len=MAX_PATH), intent(out) :: param_file, patch_file, out_file 
      integer,                 intent(out) :: numSWb, year, jday 
      real(r8),                intent(out) :: lat, lon
      real(r8),                intent(out) :: fcansno

    end subroutine read_radiation_namelist

  end interface

  call read_radiation_namelist(numSWB, year, jday, lat, lon, fcansno,          &
    param_file, patch_file, out_file)
  
  ! open log file
  logf = open_file("log.txt")

  ! get patch and parameter values, as well as orbital parameters
  !call get_parameters(param_file, rhol, rhos, taul, taus, xl, clumping_index)
  
  ! call read_patch_data(patch_file, canopy_area_profile, elai_profile,        &
  !     esai_profile, nrad_r)
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

subroutine read_radiation_namelist(numSWb, year, jday, lat, lon, fcansno,      &
    param_file, patch_file, out_file)

  use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file
  
  implicit none

  !
  ! DESCRIPTION:
  !		read in the namelist associated with the radiation unit tests and 
  !   initialize values
  !

  ! ARGUMENTS:
  character(len=MAX_PATH), intent(out) :: param_file, patch_file, out_file  ! file names
  integer,                 intent(out) :: numSWb                            ! number of shortwave bands
  integer,                 intent(out) :: year, jday                        ! year and day of year
  real(r8),                intent(out) :: lat, lon                          ! latitude and longitude [degrees]
  real(r8),                intent(out) :: fcansno                           ! fraction of canopy covered by snow [0-1]

  ! LOCALS:
  character(len=MAX_PATH) :: rad_nl = 'radiation_nl' ! radiation namelist name
  character(len=MAX_PATH) :: message                 ! Error message
  character(len=MAX_PATH) :: msg                     ! I/O Error message
  integer                 :: rad_nl_file             ! unit number for namelist
  integer                 :: ios                     ! I/O status
  

  ! Namelist of radiation parameters
  namelist /radiation/ numSWb, year, jday, lat, lon, fcansno, param_file,      &
    patch_file, out_file

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