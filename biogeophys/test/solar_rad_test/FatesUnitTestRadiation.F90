program FatesUnitTestRadiation
  !
  ! DESCRIPTION:
  !		Test the FATES radiation schemes
  !

  use FatesUnitTestIOMod,        only : logf, MAX_PATH, OpenFile
  use FatesConstantsMod,         only : r8 => fates_r8
  use EDPftvarcon,               only : EDPftvarcon_inst
  use FatesPatchMod,             only : fates_patch_type
  use FatesUnitTestRadiationMod, only : fates_rad_test

  implicit none
    
  ! LOCALS:
  type(fates_patch_type), pointer :: fatesPatch                   ! patch object
  type(fates_rad_test)            :: fatesRad                     ! rad type 

  ! PARAMETERS
  integer, parameter              :: current_tod = 0              ! current time [seconds past 0Z]
  character(len=MAX_PATH)         :: rad_nl_file = 'radiation_nl' ! radiation runtime file
  character(len=MAX_PATH)         :: pft_nl_file = 'pft_nl'       ! pft runtime file

  !---------------------------------------------------------------------------------------
  
  ! open log file
  logf = OpenFile("log.txt")

  ! initialize radiation unit test parameters
  call fatesRad%Init(EDPftvarcon_inst, rad_nl_file, pft_nl_file)

  ! initialize a patch and set values
  allocate(fatesPatch)
  call fatesRad%InitPatch(fatesPatch)

  ! ! calculate orbital values based on input year
  ! call get_orbital_vals(year, logf, eccen, mvelp, obliqr, lambm0, mvelpp)

  ! ! initialize values
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





