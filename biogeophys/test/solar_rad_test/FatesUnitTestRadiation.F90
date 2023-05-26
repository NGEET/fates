program FatesUnitTestRadiation
  !
  ! DESCRIPTION:
  !		Test the FATES radiation schemes
  !

  use FatesUnitTestIOMod,        only : logf, MAX_PATH, OpenFile, Check
  use FatesConstantsMod,         only : r8 => fates_r8
  use EDPftvarcon,               only : EDPftvarcon_inst
  use FatesPatchMod,             only : fates_patch_type
  use FatesUnitTestRadiationMod, only : fates_rad_test, ZeroPatchRadVars
  use EDSurfaceRadiationMod,     only : PatchNormanRadiation
  use FatesInterfaceTypesMod,    only : hlm_numSWb, numpft

  use netcdf

  use shr_kind_mod,              only : SHR_KIND_R8

  implicit none
    
  ! LOCALS:
  type(fates_patch_type), pointer :: fatesPatch   ! patch object
  type(fates_rad_test)            :: fatesRad     ! rad type
  integer                         :: i            ! looping index
  real(r8), allocatable           :: time_counter(:)
  real(r8), allocatable           :: alb_dir(:,:) ! surface albedo (direct) [0-1]
  real(r8), allocatable           :: alb_dif(:,:) ! surface albedo (diffuse) [0-1]
  real(r8), allocatable           :: abs_dir(:,:) ! flux absorbed by canopy per unit direct flux
  real(r8), allocatable           :: abs_dif(:,:) ! flux absorbed by canopy per unit diffuse flux
  real(r8), allocatable           :: ftdd(:,:)    ! down direct flux below canopy per unit direct flux
  real(r8), allocatable           :: ftid(:,:)    ! down diffuse flux below canopy per unit direct flx
  real(r8), allocatable           :: ftii(:,:)    ! down diffuse flux below canopy per unit diffuse flx

  ! PARAMETERS
  character(len=MAX_PATH), parameter :: rad_nl_file = 'radiation_nl' ! radiation runtime file
  character(len=MAX_PATH), parameter :: pft_nl_file = 'pft_nl'       ! pft runtime file
  integer :: ncid
  integer :: dimID_swb, dimID_time_counter
  integer :: time_counter_ID
  integer :: alb_dir_ID

  !---------------------------------------------------------------------------------------
  
  ! open log file
  logf = OpenFile("log.txt", mode='rw')

  ! initialize radiation unit test parameters
  call fatesRad%Init(EDPftvarcon_inst, rad_nl_file, pft_nl_file)

  ! allocate arrays based on input number of SWBs
  hlm_numSWb = fatesRad%num_swb
  allocate(alb_dir(fatesRad%num_steps,hlm_numSWb))
  allocate(alb_dif(fatesRad%num_steps,hlm_numSWb))
  allocate(abs_dir(fatesRad%num_steps,hlm_numSWb))
  allocate(abs_dif(fatesRad%num_steps,hlm_numSWb))
  allocate(ftdd(fatesRad%num_steps,hlm_numSWb))
  allocate(ftid(fatesRad%num_steps,hlm_numSWb))
  allocate(ftii(fatesRad%num_steps,hlm_numSWb))
  allocate(time_counter(fatesRad%num_steps))

  ! set this now
  numpft = fatesRad%num_pft

  ! initialize a patch and set values
  allocate(fatesPatch)
  call fatesRad%InitPatch(fatesPatch)

  ! calculate orbital parameters based on input year
  call fatesRad%GetOrbitalParams()
  
  ! for each time step 
  do i = 1, fatesRad%num_steps
    call fatesRad%CalcCosZ(i, fatesPatch%solar_zenith_angle)

    call ZeroPatchRadVars(fatesPatch)
    time_counter(i) = i

    fatesPatch%gnd_alb_dif(1:hlm_numSWb) = fatesRad%ground_albedo
    fatesPatch%gnd_alb_dir(1:hlm_numSWb) = fatesRad%ground_albedo
    
    if (fatesPatch%solar_zenith_angle > 0.0_r8) then 
      
      fatesPatch%solar_zenith_flag = .true.

      alb_dir(i,:) = 0._r8  
      alb_dif(i,:) = 0._r8  
      abs_dir(i,:) = 0._r8  
      abs_dif(i,:) = 0._r8  
      ftdd(i,:)    = 1._r8 
      ftid(i,:)    = 1._r8 
      ftii(i,:)    = 1._r8 

      if (maxval(fatesPatch%nrad(1,:)) == 0) then
        
        ! no leaf layers
        fatesPatch%radiation_error   = 0.0_r8
        alb_dir(i,1:hlm_numSWb)      = fatesRad%ground_albedo
        alb_dif(i,1:hlm_numSWb)      = fatesRad%ground_albedo
        ftid(i,1:hlm_numSWb)         = 0.0_r8

      else 

        ! call norman radiation scheme
        call PatchNormanRadiation(fatesPatch, alb_dir(i,:), alb_dif(i,:), abs_dir(i,:),  &
          abs_dif(i,:), ftdd(i,:), ftid(i,:), ftii(i,:))
      end if 
    end if
  end do 

  ! write output
  !call fatesRad%WriteRadiationData(alb_dir)

!print *, alb_dir(:,2)
  
call Check(NF90_CREATE(TRIM(fatesRad%out_file),NF90_CLOBBER,ncid))
call Check(nf90_def_dim(ncid,"time_counter", fatesRad%num_steps ,dimID_time_counter))

! Variable definitions
call Check(NF90_DEF_VAR(ncid,"time_counter", NF90_INT, (/dimID_time_counter/),time_counter_ID))
call Check(NF90_DEF_VAR(ncid,"alb_dir", NF90_DOUBLE, (/dimID_time_counter/), alb_dir_ID))

! Attributes :
call Check(NF90_PUT_ATT(ncid, time_counter_ID, "time_origin", "2000-06-14 00:00:00"))
call Check(NF90_PUT_ATT(ncid, time_counter_ID, "units", "seconds since 2000-06-14 00:00:00"))
call Check(NF90_PUT_ATT(ncid, time_counter_ID, "calendar","gregorian"))
call Check(NF90_PUT_ATT(ncid, time_counter_ID, "long_name","Time axis"))

call Check(NF90_PUT_ATT(ncid, alb_dir_ID, "coordinates", "time_counter"))
call Check(NF90_PUT_ATT(ncid, alb_dir_ID, "units", "0-1"))
call Check(NF90_PUT_ATT(ncid, alb_dir_ID, "long_name", "direct albedo"))
call Check(NF90_PUT_ATT(ncid, alb_dir_ID, "standard_name", "alb_dir"))

call Check(NF90_ENDDEF(ncid))

call Check(NF90_PUT_VAR(ncid, time_counter_ID, time_counter))
call Check(NF90_PUT_VAR(ncid, alb_dir_ID, alb_dir(:,2)))

call Check(nf90_close(ncid))

close(logf)

end program FatesUnitTestRadiation





