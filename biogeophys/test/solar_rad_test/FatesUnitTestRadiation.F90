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

  end interface

  !:...........................................................................:
  
  ! open log file
  logf = open_file("log.txt")

  ! read in namelist to get some runtime parameters
  call read_radiation_namelist(year, jday, lat, lon, fcansno, param_file,      &
    patch_file, out_file)
  
  ! read in parameter file and initialize EDPFTvarcon_inst 
  


  ! read in patch data
 !call read_patch_data(patch_file, canopy_area_profile, elai_profile,          &
 !   esai_profile, nrad_r)

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

  use FatesUnitTestIOMod, only : MAX_PATH, MAX_CHAR, open_file, logf
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
    write(logf, '(A, I6, A, A)') "Error reading radiation namelist file",      &
      ios, "IOMSG: ", msg
    stop "Stopped"
  end if

  close(rad_nl_file)

end subroutine read_radiation_namelist

!:.............................................................................:



