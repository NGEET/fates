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
  use EDSurfaceRadiationMod,     only : PatchNormanRadiation
  use FatesInterfaceTypesMod,    only : hlm_numSWb

  use shr_kind_mod,              only : SHR_KIND_R8

  implicit none
    
  ! LOCALS:
  type(fates_patch_type), pointer :: fatesPatch     ! patch object
  type(fates_rad_test)            :: fatesRad       ! rad type
  integer                         :: i              ! looping index
  real(SHR_KIND_R8)               :: cosz           ! cosine of solar zenith angle
  real(r8), allocatable           :: alb_parb_out(:) 
  real(r8), allocatable           :: albi_parb_out(:)
  real(r8), allocatable           :: fabd_parb_out(:)
  real(r8), allocatable           :: fabi_parb_out(:)
  real(r8), allocatable           :: ftdd_parb_out(:)
  real(r8), allocatable           :: ftid_parb_out(:)
  real(r8), allocatable           :: ftii_parb_out(:)

  ! PARAMETERS
  integer,                 parameter :: current_tod = 0              ! current time [seconds past 0Z]
  character(len=MAX_PATH), parameter :: rad_nl_file = 'radiation_nl' ! radiation runtime file
  character(len=MAX_PATH), parameter :: pft_nl_file = 'pft_nl'       ! pft runtime file

  !---------------------------------------------------------------------------------------
  
  ! open log file
  logf = OpenFile("log.txt")

  ! initialize radiation unit test parameters
  call fatesRad%Init(EDPftvarcon_inst, rad_nl_file, pft_nl_file)

  hlm_numSWb = fatesRad%num_swb
  allocate(alb_parb_out(hlm_numSWb))
  allocate(albi_parb_out(hlm_numSWb))
  allocate(fabd_parb_out(hlm_numSWb))
  allocate(fabi_parb_out(hlm_numSWb))
  allocate(ftdd_parb_out(hlm_numSWb))
  allocate(ftid_parb_out(hlm_numSWb))
  allocate(ftii_parb_out(hlm_numSWb))

  ! initialize a patch and set values
  allocate(fatesPatch)
  call fatesRad%InitPatch(fatesPatch)

  ! calculate orbital parameters based on input year
  call fatesRad%GetOrbitalParams()
  
  ! for each half-hourly time step in the day
  do i = 0, 47
    call fatesRad%CalcCosZ(i, cosz)
    
    if (cosz > 0.0_r8) then 
      ! call norman radiation scheme
      call PatchNormanRadiation(fatesPatch, alb_parb_out, albi_parb_out, fabd_parb_out,  &
        fabi_parb_out, ftdd_parb_out, ftid_parb_out, ftii_parb_out)
    end if

  end do 

  ! call write_radiation_data(out_file, k_dir, cosz)

  close(logf)

end program FatesUnitTestRadiation





