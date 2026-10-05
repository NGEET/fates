program FatesLeafLevelPhoto
  !
  ! DESCRIPTION:
  ! Leaf-level photosynthesis sensitivity sweep: five independent sweeps (PAR,
  ! CO2, VPD, leaf temperature, and soil water content), each varying one
  ! driver variable while holding everything else at a fixed default,
  ! evaluated for a single PFT
  !
  ! The Kumarathunge et al. (2019) temperature-acclimation model uses two running-mean
  ! reference temperatures (t_growth, t_home) in addition to the instantaneous leaf
  ! temperature being swept. Both are held fixed at the same default leaf temperature
  ! for every sweep and every point.
  !
  ! The soil water content sweep does not vary btran directly. It sweeps a soil
  ! water content fraction in [0, 1], maps it onto a soil matric potential, since 
  ! these drivers have no soil texture from which a real retention curve could be built), 
  ! and derives btran from that via BtranFromSMP, which is production's own 
  ! EDBtranMod.F90::btran_ed formula specialized to the single unfrozen layer at 
  ! root fraction 1 these drivers assume
  !

  use FatesConstantsMod,           only : r8 => fates_r8
  use FatesConstantsMod,           only : t_water_freeze_k_1atm
  use FatesArgumentUtils,          only : command_line_arg
  use FatesUnitTestParamReaderMod, only : ReadParameters, CheckLeafRespParams
  use EDPftvarcon,                 only : EDPftvarcon_inst
  use FatesParameterDerivedMod,    only : param_derived
  use FatesFactoryMod,             only : InitializeGlobals
  use LeafBiophysicsMod,           only : lb_params
  use LeafBiophysicsMod,           only : QSat
  use FatesInterfaceTypesMod,      only : hlm_maintresp_leaf_model
  use FatesConstantsMod,           only : lmrmodel_ryan_1991, lmrmodel_atkin_etal_2017
  use LeafBiophysicsMod,           only : FvCB1980, medlyn_model, net_assim_model
  use LeafBiophysicsMod,           only : photosynth_acclim_model_kumarathunge_etal_2019
  use FatesTestLeafPhotoMod,       only : EvaluateLeafPhotosynthesis, LeafNitrogenContent
  use FatesTestEnvironmentMod,     only : environment_type
  use FatesTestEnvironmentMod,     only : default_nscaler, default_rdark_scaler
  use FatesTestEnvironmentMod,     only : CanopyVaporPressure
  use FatesTestEnvironmentMod,     only : BtranFromSMP, SoilMatricPotential
  use FatesUnitTestIOMod,          only : OpenNCFile, RegisterNCDims, CloseNCFile
  use FatesUnitTestIOMod,          only : WriteVar, RegisterVarAtts, EndNCDef
  use FatesUnitTestIOMod,          only : type_double, type_int

  implicit none

  ! LOCALS:
  type(environment_type)        :: env                  ! prescribed atmospheric boundary conditions
  character(len=:), allocatable :: param_file           ! input parameter file
  character(len=:), allocatable :: out_file             ! output file name
  character(len=:), allocatable :: pft_arg              ! command-line pft argument
  real(r8),         allocatable :: par_vals(:)          ! swept PAR values [umol/m2/s]
  real(r8),         allocatable :: co2_vals(:)          ! swept CO2 partial pressure values [Pa]
  real(r8),         allocatable :: vpd_vals(:)          ! swept leaf-to-air VPD values [Pa]
  real(r8),         allocatable :: temp_vals(:)         ! swept leaf temperature values [K]
  real(r8),         allocatable :: soilfrac_vals(:)     ! swept soil water content, fraction of saturation [0-1]
  real(r8),         allocatable :: veg_esat_bytemp(:)   ! saturation vapor pressure at each swept leaf temperature [Pa]
  real(r8),         allocatable :: can_vpress_bytemp(:) ! canopy air vapor pressure at each swept leaf temperature, fixed default VPD [Pa]
  real(r8),         allocatable :: can_vpress_byvpd(:)  ! canopy air vapor pressure at each swept VPD (= veg_esat - vpd), fixed default leaf temperature [Pa]
  real(r8),         allocatable :: btran_bysoilfrac(:)  ! btran derived from each swept soil water content fraction [0-1]
  real(r8)                      :: lnc_top              ! leaf N content at the canopy top [gN/m2 leaf]
  real(r8)                      :: vcmax25top           ! top-of-canopy carboxylation rate at 25degC [umol/m2/s]
  real(r8)                      :: jmax25top            ! top-of-canopy electron transport rate at 25degC [umol/m2/s]
  real(r8)                      :: kp25top              ! top-of-canopy initial slope of C4 CO2 response at 25degC [umol/m2/s]
  real(r8)                      :: smpsc                ! soil matric potential at full stomatal closure [mm, negative]
  real(r8)                      :: smpso                ! soil matric potential at full stomatal opening [mm, negative]
  real(r8)                      :: qs_dummy             ! saturation specific humidity output from QSat (unused here)
  real(r8)                      :: smp                  ! soil matric potential at a swept soil water content fraction [mm]
  integer                       :: target_pft           ! PFT to use
  integer                       :: n_par                ! PAR sweep array size 
  integer                       :: n_co2                ! CO2 sweep array size 
  integer                       :: n_vpd                ! VPD sweep array size 
  integer                       :: n_temp               ! temperature sweep array size 
  integer                       :: n_soilfrac           ! soilfrac sweep array size
  integer                       :: i                    ! looping index
  
  ! sweep output, each dimensioned (n_sweep_points)
  real(r8), allocatable :: anet_bypar(:), agross_bypar(:), gs_bypar(:), ci_bypar(:)
  real(r8), allocatable :: anet_byco2(:), agross_byco2(:), gs_byco2(:), ci_byco2(:)
  real(r8), allocatable :: anet_byvpd(:), agross_byvpd(:), gs_byvpd(:), ci_byvpd(:)
  real(r8), allocatable :: anet_bytemp(:), agross_bytemp(:), gs_bytemp(:), ci_bytemp(:)
  real(r8), allocatable :: anet_bysoilfrac(:), agross_bysoilfrac(:), gs_bysoilfrac(:), ci_bysoilfrac(:)
  
  ! CONSTANTS:
  integer, parameter :: default_pft = 1 ! default PFT index to evaluate (1-based)
  
  ! sweep ranges
  real(r8), parameter :: min_temp = 8.0_r8,    max_temp = 40.0_r8,    temp_inc = 0.5_r8   ! [degC]
  real(r8), parameter :: min_par  = 0.0_r8,    max_par  = 2500.0_r8,  par_inc  = 5.0_r8   ! [umol/m2/s]
  real(r8), parameter :: min_vpd  = 500.0_r8,  max_vpd  = 2500.0_r8,  vpd_inc  = 20.0_r8  ! [Pa] (0.5-2.5 kPa)
  real(r8), parameter :: min_co2  = 50.0_r8,  max_co2  = 1500.0_r8,  co2_inc  = 5.0_r8   ! [umol/mol]
  real(r8), parameter :: soilfrac_inc = 0.02_r8 ! [0-1]
  
  ! --------------------------------------------------------------------------------------
  ! Initial setup - read parameter file, initialize globals
  ! --------------------------------------------------------------------------------------
  
  ! read in parameter file name from command line
  param_file = command_line_arg(1)
  
  ! output file name, depends on either arg2 or is just default
  if (command_argument_count() >= 2) then
    out_file = trim(command_line_arg(2))
  else
    out_file = 'leaf_level_photo_out.nc'
  end if

  ! pft, depends on either arg3 or is just default
  if (command_argument_count() >= 3) then 
    pft_arg = command_line_arg(3)
    read(pft_arg,*) target_pft
  else 
    target_pft = default_pft
  end if
  
  ! read in parameter file
  call ReadParameters(param_file)
  
  ! initialize global values
  call InitializeGlobals()
  
  ! host-model-namelist-controlled leaf biophysics switches
  hlm_maintresp_leaf_model = lmrmodel_ryan_1991
  lb_params%electron_transport_model = FvCB1980
  lb_params%stomatal_model = medlyn_model
  lb_params%stomatal_assim_model = net_assim_model
  lb_params%photo_tempsens_model = photosynth_acclim_model_kumarathunge_etal_2019
  
  ! do a parameter check if for Atkin parameters
  call CheckLeafRespParams()
  
  ! leaf N content for target_pft - constant for the whole run
  lnc_top = LeafNitrogenContent(target_pft)
  
  ! photosynthetic capacity parameters for target_pft
  vcmax25top = EDPftvarcon_inst%vcmax25top(target_pft,1)
  jmax25top = param_derived%jmax25top(target_pft,1)
  kp25top = param_derived%kp25top(target_pft,1)
  
  ! soil stress parameters for target_pft
  smpsc = EDPftvarcon_inst%smpsc(target_pft)
  smpso = EDPftvarcon_inst%smpso(target_pft)
  
  ! set atmospheric defaults
  call env%Init()
  
  ! --------------------------------------------------------------------------------------
  ! build the swept-value arrays
  ! --------------------------------------------------------------------------------------
  
  n_par = int((max_par - min_par)/par_inc) + 1
  n_co2 = int((max_co2 - min_co2)/co2_inc) + 1
  n_vpd  = int((max_vpd - min_vpd)/vpd_inc) + 1
  n_temp = int((max_temp - min_temp)/temp_inc) + 1
  n_soilfrac = int((1.0_r8 - 0.0_r8)/soilfrac_inc) + 1
  
  allocate(par_vals(n_par))
  allocate(co2_vals(n_co2))
  allocate(vpd_vals(n_vpd))
  allocate(temp_vals(n_temp))
  allocate(soilfrac_vals(n_soilfrac))
  
  do i = 1, n_par
    par_vals(i) = min_par + par_inc*real(i-1, r8)
  end do
  do i = 1, n_co2
    ! ppm -> Pa, at the default (Init()-prescribed) canopy pressure
    co2_vals(i) = ((min_co2 + co2_inc*real(i-1, r8))/1.0e6_r8) * env%can_press
  end do
  do i = 1, n_vpd
    vpd_vals(i) = min_vpd + vpd_inc*real(i-1, r8)
  end do
  do i = 1, n_temp
    temp_vals(i) = (min_temp + temp_inc*real(i-1, r8)) + t_water_freeze_k_1atm
  end do
  do i = 1, n_soilfrac
    soilfrac_vals(i) = 0.0_r8 + soilfrac_inc*real(i-1, r8)
  end do
  
  allocate(veg_esat_bytemp(n_temp), can_vpress_bytemp(n_temp))
  allocate(can_vpress_byvpd(n_vpd))
  allocate(btran_bysoilfrac(n_soilfrac))
  
  allocate(anet_bypar(n_par), agross_bypar(n_par), gs_bypar(n_par), ci_bypar(n_par))
  allocate(anet_byco2(n_co2), agross_byco2(n_co2), gs_byco2(n_co2), ci_byco2(n_co2))
  allocate(anet_byvpd(n_vpd), agross_byvpd(n_vpd), gs_byvpd(n_vpd), ci_byvpd(n_vpd))
  allocate(anet_bytemp(n_temp), agross_bytemp(n_temp), gs_bytemp(n_temp), ci_bytemp(n_temp))
  allocate(anet_bysoilfrac(n_soilfrac), agross_bysoilfrac(n_soilfrac), gs_bysoilfrac(n_soilfrac), ci_bysoilfrac(n_soilfrac))
  
  ! ---------------------------------------------------------------------
  ! PAR sweep
  ! ---------------------------------------------------------------------
  do i = 1, n_par
    call EvaluateLeafPhotosynthesis(target_pft, par_vals(i), env%tempk,   &
      env%tempk, env%tempk, env%can_press, env%can_co2_ppress,            &
      env%can_o2_ppress, env%veg_esat, env%can_vpress, env%gb,            &
      default_nscaler, default_rdark_scaler, env%dayl_factor, env%btran,  &
      vcmax25top, jmax25top, kp25top, lnc_top, agross_bypar(i),           &
      anet_bypar(i), gs_bypar(i), ci_bypar(i))
  end do
  
  ! ---------------------------------------------------------------------
  ! CO2 sweep
  ! ---------------------------------------------------------------------
  do i = 1, n_co2
    call EvaluateLeafPhotosynthesis(target_pft, env%par, env%tempk,      &
      env%tempk, env%tempk, env%can_press, co2_vals(i),                  &
      env%can_o2_ppress, env%veg_esat, env%can_vpress, env%gb,           &
      default_nscaler, default_rdark_scaler, env%dayl_factor, env%btran, &
      vcmax25top, jmax25top, kp25top, lnc_top, agross_byco2(i),          &
      anet_byco2(i), gs_byco2(i), ci_byco2(i))
  end do
  
  ! ---------------------------------------------------------------------
  ! VPD sweep - leaf temperature fixed at the default, so veg_esat is
  ! constant and can_vpress is derived directly from the swept VPD
  ! ---------------------------------------------------------------------
  do i = 1, n_vpd
    can_vpress_byvpd(i) = CanopyVaporPressure(env%veg_esat, vpd=vpd_vals(i))
    call EvaluateLeafPhotosynthesis(target_pft, env%par, env%tempk,      &
      env%tempk, env%tempk, env%can_press, env%can_co2_ppress,           &
      env%can_o2_ppress, env%veg_esat, can_vpress_byvpd(i), env%gb,      &
      default_nscaler, default_rdark_scaler, env%dayl_factor, env%btran, &
      vcmax25top, jmax25top, kp25top, lnc_top, agross_byvpd(i),          &
      anet_byvpd(i), gs_byvpd(i), ci_byvpd(i))
  end do
  
  ! ---------------------------------------------------------------------
  ! leaf temperature sweep - t_growth/t_home held at the default leaf
  ! temperature throughout; VPD held fixed at the default as
  ! leaf temperature varies
  ! ---------------------------------------------------------------------
  do i = 1, n_temp
    call QSat(temp_vals(i), env%can_press, qs_dummy, veg_esat_bytemp(i))
    can_vpress_bytemp(i) = CanopyVaporPressure(veg_esat_bytemp(i))
    call EvaluateLeafPhotosynthesis(target_pft, env%par, temp_vals(i),     &
      env%tempk, env%tempk, env%can_press, env%can_co2_ppress,             &
      env%can_o2_ppress, veg_esat_bytemp(i), can_vpress_bytemp(i), env%gb, &
      default_nscaler, default_rdark_scaler, env%dayl_factor, env%btran,   &
      vcmax25top, jmax25top, kp25top, lnc_top, agross_bytemp(i),           &
      anet_bytemp(i), gs_bytemp(i), ci_bytemp(i))
  end do
  
  ! ---------------------------------------------------------------------
  ! soil water content sweep - btran derived from the real smpsc/smpso
  ! ramp at each swept fraction 
  ! ---------------------------------------------------------------------
  do i = 1, n_soilfrac
    smp = SoilMatricPotential(soilfrac_vals(i), smpsc)
    btran_bysoilfrac(i) = BtranFromSMP(smp, smpsc, smpso)
    call EvaluateLeafPhotosynthesis(target_pft, env%par, env%tempk,  &
      env%tempk, env%tempk, env%can_press, env%can_co2_ppress,       &
      env%can_o2_ppress, env%veg_esat, env%can_vpress, env%gb,       &
      default_nscaler, default_rdark_scaler, env%dayl_factor,        &
      btran_bysoilfrac(i), vcmax25top, jmax25top, kp25top, lnc_top,  &
      agross_bysoilfrac(i), anet_bysoilfrac(i), gs_bysoilfrac(i), ci_bysoilfrac(i))
  end do
  
  ! --------------------------------------------------------------------------------------
  ! write output
  ! --------------------------------------------------------------------------------------
  
  call WriteOutput()
  
contains 

  ! ==========================================================================

  subroutine WriteOutput()
    !
    ! DESCRIPTION:
    ! Writes every sweep's swept values, derived vapor-pressure/btran
    ! diagnostics (temperature/VPD/soil-water-content sweeps only), and
    ! agross/anet/gs/ci to netCDF.

    ! LOCALS:
    integer              :: ncid           ! netcdf file id
    character(len=20)    :: dim_names(5)   ! dimension names
    integer              :: dimIDs(5)      ! dimension IDs
    integer              :: parID, co2ID, vpdID, tempID, soilfracID
    integer              :: vegesatbytempID, canvpressbytempID
    integer              :: canvpressbyvpdID
    integer              :: btranbysoilfracID
    integer              :: agrossbyparID, anetbyparID, gsbyparID, cibyparID
    integer              :: agrossbyco2ID, anetbyco2ID, gsbyco2ID, cibyco2ID
    integer              :: agrossbyvpdID, anetbyvpdID, gsbyvpdID, cibyvpdID
    integer              :: agrossbytempID, anetbytempID, gsbytempID, cibytempID
    integer              :: agrossbysoilfracID, anetbysoilfracID, gsbysoilfracID, cibysoilfracID

    dim_names = [character(len=20) :: 'par', 'co2', 'vpd', 'temp', 'soilfrac']

    call OpenNCFile(trim(out_file), ncid, 'readwrite')

    call RegisterNCDims(ncid, dim_names, (/n_par, n_co2, n_vpd, n_temp, n_soilfrac/), 5, dimIDs)

    call RegisterVarAtts(ncid, dim_names(1), dimIDs(1:1), type_double, 'umol m-2 s-1',  &
      'swept incident PAR', parID)
    call RegisterVarAtts(ncid, dim_names(2), dimIDs(2:2), type_double, 'Pa',            &
      'swept CO2 partial pressure', co2ID)
    call RegisterVarAtts(ncid, dim_names(3), dimIDs(3:3), type_double, 'Pa',            &
      'swept leaf-to-air vapor pressure deficit', vpdID)
    call RegisterVarAtts(ncid, dim_names(4), dimIDs(4:4), type_double, 'K',             &
      'swept leaf temperature', tempID)
    call RegisterVarAtts(ncid, dim_names(5), dimIDs(5:5), type_double, '-',             &
      'swept soil water content, fraction of saturation', soilfracID)

    call RegisterVarAtts(ncid, 'veg_esat_bytemp', dimIDs(4:4), type_double, 'Pa',       &
      'saturation vapor pressure at each swept leaf temperature', vegesatbytempID)
    call RegisterVarAtts(ncid, 'can_vpress_bytemp', dimIDs(4:4), type_double, 'Pa',     &
      'canopy air vapor pressure at each swept leaf temperature, fixed default VPD',    &
      canvpressbytempID)
    call RegisterVarAtts(ncid, 'can_vpress_byvpd', dimIDs(3:3), type_double, 'Pa',      &
      'canopy air vapor pressure at each swept VPD, fixed default leaf temperature',    &
      canvpressbyvpdID)
    call RegisterVarAtts(ncid, 'btran_bysoilfrac', dimIDs(5:5), type_double, '-',       &
      'btran derived from each swept soil water content fraction', btranbysoilfracID)

    call RegisterVarAtts(ncid, 'agross_bypar', (/dimIDs(1)/), type_double,   &
      'umolC m-2 s-1', 'gross photosynthesis vs. PAR', agrossbyparID,                   &
      coordinates='par')
    call RegisterVarAtts(ncid, 'anet_bypar', (/dimIDs(1)/), type_double,     &
      'umolC m-2 s-1', 'net photosynthesis vs. PAR', anetbyparID, coordinates='par')
    call RegisterVarAtts(ncid, 'gs_bypar', (/dimIDs(1)/), type_double,       &
      'umol H2O m-2 s-1', 'stomatal conductance vs. PAR', gsbyparID,                    &
      coordinates='par')
    call RegisterVarAtts(ncid, 'ci_bypar', (/dimIDs(1)/), type_double, 'Pa', &
      'intracellular CO2 vs. PAR', cibyparID, coordinates='par')

    call RegisterVarAtts(ncid, 'agross_byco2', (/dimIDs(2)/), type_double,   &
      'umolC m-2 s-1', 'gross photosynthesis vs. CO2', agrossbyco2ID,                   &
      coordinates='co2')
    call RegisterVarAtts(ncid, 'anet_byco2', (/dimIDs(2)/), type_double,     &
      'umolC m-2 s-1', 'net photosynthesis vs. CO2', anetbyco2ID, coordinates='co2')
    call RegisterVarAtts(ncid, 'gs_byco2', (/dimIDs(2)/), type_double,       &
      'umol H2O m-2 s-1', 'stomatal conductance vs. CO2', gsbyco2ID,                    &
      coordinates='co2')
    call RegisterVarAtts(ncid, 'ci_byco2', (/dimIDs(2)/), type_double, 'Pa', &
      'intracellular CO2 vs. CO2', cibyco2ID, coordinates='co2')

    call RegisterVarAtts(ncid, 'agross_byvpd', (/dimIDs(3)/), type_double,   &
      'umolC m-2 s-1', 'gross photosynthesis vs. VPD', agrossbyvpdID,                   &
      coordinates='vpd')
    call RegisterVarAtts(ncid, 'anet_byvpd', (/dimIDs(3)/), type_double,     &
      'umolC m-2 s-1', 'net photosynthesis vs. VPD', anetbyvpdID, coordinates='vpd')
    call RegisterVarAtts(ncid, 'gs_byvpd', (/dimIDs(3)/), type_double,       &
      'umol H2O m-2 s-1', 'stomatal conductance vs. VPD', gsbyvpdID,                    &
      coordinates='vpd')
    call RegisterVarAtts(ncid, 'ci_byvpd', (/dimIDs(3)/), type_double, 'Pa', &
      'intracellular CO2 vs. VPD', cibyvpdID, coordinates='vpd')

    call RegisterVarAtts(ncid, 'agross_bytemp', (/dimIDs(4)/), type_double,  &
      'umolC m-2 s-1', 'gross photosynthesis vs. leaf temperature', agrossbytempID,     &
      coordinates='temperature')
    call RegisterVarAtts(ncid, 'anet_bytemp', (/dimIDs(4)/), type_double,    &
      'umolC m-2 s-1', 'net photosynthesis vs. leaf temperature', anetbytempID,         &
      coordinates='temperature')
    call RegisterVarAtts(ncid, 'gs_bytemp', (/dimIDs(4)/), type_double,      &
      'umol H2O m-2 s-1', 'stomatal conductance vs. leaf temperature', gsbytempID,      &
      coordinates='temperature')
    call RegisterVarAtts(ncid, 'ci_bytemp', (/dimIDs(4)/), type_double,      &
      'Pa', 'intracellular CO2 vs. leaf temperature', cibytempID,                       &
      coordinates='temperature')

    call RegisterVarAtts(ncid, 'agross_bysoilfrac', (/dimIDs(5)/),           &
      type_double, 'umolC m-2 s-1', 'gross photosynthesis vs. soil water content',      &
      agrossbysoilfracID, coordinates='soilfrac')
    call RegisterVarAtts(ncid, 'anet_bysoilfrac', (/dimIDs(5)/),             &
      type_double, 'umolC m-2 s-1', 'net photosynthesis vs. soil water content',        &
      anetbysoilfracID, coordinates='soilfrac')
    call RegisterVarAtts(ncid, 'gs_bysoilfrac', (/dimIDs(5)/), type_double,  &
      'umol H2O m-2 s-1', 'stomatal conductance vs. soil water content',                &
      gsbysoilfracID, coordinates='soilfrac')
    call RegisterVarAtts(ncid, 'ci_bysoilfrac', (/dimIDs(5)/), type_double,  &
      'Pa', 'intracellular CO2 vs. soil water content', cibysoilfracID,                 &
      coordinates='soilfrac')

    call EndNCDef(ncid)

    call WriteVar(ncid, parID, par_vals(:))
    call WriteVar(ncid, co2ID, co2_vals(:))
    call WriteVar(ncid, vpdID, vpd_vals(:))
    call WriteVar(ncid, tempID, temp_vals(:))
    call WriteVar(ncid, soilfracID, soilfrac_vals(:))

    call WriteVar(ncid, vegesatbytempID, veg_esat_bytemp(:))
    call WriteVar(ncid, canvpressbytempID, can_vpress_bytemp(:))
    call WriteVar(ncid, canvpressbyvpdID, can_vpress_byvpd(:))
    call WriteVar(ncid, btranbysoilfracID, btran_bysoilfrac(:))

    call WriteVar(ncid, agrossbyparID, agross_bypar(:))
    call WriteVar(ncid, anetbyparID, anet_bypar(:))
    call WriteVar(ncid, gsbyparID, gs_bypar(:))
    call WriteVar(ncid, cibyparID, ci_bypar(:))

    call WriteVar(ncid, agrossbyco2ID, agross_byco2(:))
    call WriteVar(ncid, anetbyco2ID, anet_byco2(:))
    call WriteVar(ncid, gsbyco2ID, gs_byco2(:))
    call WriteVar(ncid, cibyco2ID, ci_byco2(:))

    call WriteVar(ncid, agrossbyvpdID, agross_byvpd(:))
    call WriteVar(ncid, anetbyvpdID, anet_byvpd(:))
    call WriteVar(ncid, gsbyvpdID, gs_byvpd(:))
    call WriteVar(ncid, cibyvpdID, ci_byvpd(:))

    call WriteVar(ncid, agrossbytempID, agross_bytemp(:))
    call WriteVar(ncid, anetbytempID, anet_bytemp(:))
    call WriteVar(ncid, gsbytempID, gs_bytemp(:))
    call WriteVar(ncid, cibytempID, ci_bytemp(:))

    call WriteVar(ncid, agrossbysoilfracID, agross_bysoilfrac(:))
    call WriteVar(ncid, anetbysoilfracID, anet_bysoilfrac(:))
    call WriteVar(ncid, gsbysoilfracID, gs_bysoilfrac(:))
    call WriteVar(ncid, cibysoilfracID, ci_bysoilfrac(:))

    call CloseNCFile(ncid)

  end subroutine WriteOutput
  
end program FatesLeafLevelPhoto

! ----------------------------------------------------------------------------------------
