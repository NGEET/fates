module CrownFireEquationsMod

  ! ============================================================================
  ! Helper methods for FATES crown fire model
  ! Most are from Scott & Reinhardt 2001 &
  ! Giuseppe 2024
  ! ============================================================================

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nearzero


  implicit none
  private

  public :: PassiveCrownFireIntensity
  public :: HeatReleasePerArea
  public :: CrowningIndex
  public :: CrownFireIntensity
  public :: LiveFuelMoistureContent
  public :: MaxHeight
  public :: CrownFireBehaveFM10
  public :: BiomassBin
  public :: CrownFireCFB

  contains

  real(r8) function PassiveCrownFireIntensity(canopy_base_height, canopy_water_content)
    ! DESCRIPTION:
    ! Calculate the energy threshold for igniting crown fuels [kW/m or kJ/m/s]
    ! EQ. 11 in Scott & Reinhardt 2001

    ! ARGUMENTS:
    real(r8), intent(in) :: canopy_base_height   ! canopy base height at which biomass density > minimum density 0.011 kg/m3
    real(r8), intent(in) :: canopy_water_content ! canopy water content [%]
  

    ! Locals:
    real(r8)            :: crown_ignition_energy  ! surface fire intensity required to ignite crown fuels [kJ/kg]

    ! Note: crown_ignition_energy to be calculated based on PFT foliar moisture content from FATES-Hydro
    ! or create foliar moisture % based on BTRAN
    ! Use foliar_moisture(currentCohort%pft) and compute weighted PFT average with Eq 3 Van Wagner 1977
    ! in place of canopy_water_content parameter

    ! Eq 3 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
    crown_ignition_energy = 460.0_r8 + 25.9_r8 * canopy_water_content

    ! Crown fuel ignition potential (kW/m), Eq 4 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
    ! FI = (Czh)**3/2 where z=canopy base height,h=heat of crown ignite energy, FI=fire intensity
    ! 0.01 = C, empirical constant Van Wagner 1977 Eq 4 for 6m canopy base height, 100% FMC, FI 2500kW/m
    ! passive_crown_FI = min fire intensity to ignite canopy fuel (kW/m or kJ/m/s)
    PassiveCrownFireIntensity = (0.01_r8 * canopy_base_height * crown_ignition_energy)**1.5_r8

  end function PassiveCrownFireIntensity

  !---------------------------------------------------------------------------------------
  real(r8) function HeatReleasePerArea(SAV, i_r)
    !
    ! DESCRIPTION:
    ! Calculate heat release per unit area of surface fire (HPA, kJ/m2)
    ! EQ. 2 in Scott & Reinhardt 2001

    ! ARGUMENTS:
    real(r8), intent(in)  :: SAV      ! fuel surface area to volume ratio [/cm]
    real(r8), intent(in)  :: i_r      ! reaction intensity [kJ/m2/min]


    ! Locals:
    real(r8)   :: time_r       ! fire residence time [min]
    if (SAV < nearzero) then
      time_r = 0.0_r8
    else
      time_r = 12.595_r8 / SAV  ! EQ 3 in Scott & Reinhardt 2001
    end if

    HeatReleasePerArea = i_r * time_r

  end function HeatReleasePerArea


  !---------------------------------------------------------------------------------------
  real(r8) function CrowningIndex(eps, q_ig, i_r, canopy_bulk_density)
  !
  ! DESCRIPTION:
  ! Calculate open wind speed [km/hr] at which a fully active crown fire is sustained
  ! EQ. 20 in Scott & Reinhardt 2001
  ! XLG: currently we are ignoring slope effect 
  !
    ! ARGUMENTS:
    real(r8), intent(in) :: eps                  ! effective heating number [unitless]
    real(r8), intent(in) :: q_ig                 ! heat of preignition [kJ/kg] 
    real(r8), intent(in) :: i_r                  ! reaction intensity [kJ/m2/min]
    real(r8), intent(in) :: canopy_bulk_density  ! canopy fuel bulk density [kg biomass / m3]
    ! Locals:
    real(r8)   :: CI_temp      ! temporary variables

    if(i_r <= nearzero .or. canopy_bulk_density <= nearzero) then
      CI_temp = 0.0_r8
    else
      CI_temp = (164.8_r8 * eps * q_ig) / (i_r * canopy_bulk_density) - 1.0_r8
    end if
    CrowningIndex = 0.0457_r8 * ((CI_temp / 0.001612_r8)**0.7_r8)

  end function CrowningIndex


  !---------------------------------------------------------------------------------------

  real(r8) function CrownFireIntensity(HPA, canopy_fuel,patch_area, CFB, ROS_final)
  !
  ! Description
  ! Calculate fire intentisy for crown fire using
  ! EQ. 22 in Scott & Reinhardt 2001
  !
  ! ARGUMENTS:
    real(r8), intent(in) :: HPA             ! heat release per unit area [kJ/m2]
    real(r8), intent(in) :: canopy_fuel     ! canopy fuel load [kg biomass] 
    real(r8), intent(in) :: patch_area      ! current patch area, to calculate fuel load density in kg/m2
    real(r8), intent(in) :: CFB             ! crown fraction burnt [fraction]
    real(r8), intent(in) :: ROS_final       ! final rate of spread after a crown fire happens [m/min]
    ! Locals:
    real(r8), parameter  :: H_canopy = 18000.0_r8 ! heat yield for canopy fuels [kJ/kg biomass]

    CrownFireIntensity = (HPA + (canopy_fuel / patch_area * H_canopy * CFB)) * &
                         ROS_final / 60.0_r8 

  end function CrownFireIntensity

  !---------------------------------------------------------------------------------------

  real(r8) function LiveFuelMoistureContent(lai, smp, min_lfmc, coeff_lfmc, &
                    smp_alpha, lai_beta, gamma_int)
    !
    ! DESCRIPTION
    ! Calculates live fuel moisture content  
    ! using EQ. 4 in McNorton and Giuseppe 2024 'A global fuel characteristic model
    ! and dataset for wildfire prediction'
    ! LFMC = max_lfmc - min_lfmc * exp(-swc_alpha*swc + lai_beta*lai + gamma_int*swc*lai)
    ! swc is mass-based soil water content in the original Eq., we now switch to soil matric potential (smp).
    ! This switch require refitting this equation to min_lfmc + coeff_lfmc*exp(smp_alpha*smp + lai_beta*lai +
    ! gamma_int*smp*lai)
    ! we can make all parameters to be PFT specific 
    !
    ! ARGUMENTS:
    real(r8),         intent(in)    :: lai                   ! cohort level leaf area index [m2 m-2]
    real(r8),         intent(in)    :: smp                   ! soil matric potential [MPa]
    real(r8),         intent(in)    :: min_lfmc              ! min. live fuel moisure content [%]
    real(r8),         intent(in)    :: coeff_lfmc            ! this will adds to the min to increase LMFC based on soil water and plant phenology [unitelss] 
    real(r8),         intent(in)    :: smp_alpha             ! model coef. associated with swc [unitless]
    real(r8),         intent(in)    :: lai_beta              ! model coef. associated with lai [unitless]
    real(r8),         intent(in)    :: gamma_int             ! model coef. for interaction effect of swc and lai on LFMC [unitless]

    ! Locals:
    real(r8)          :: effect_temp        ! the temporary variable for calculating effect of SWC and LAI

    effect_temp = smp_alpha*smp + lai_beta*lai + gamma_int*smp*lai
    LiveFuelMoistureContent = min_lfmc + coeff_lfmc*exp(effect_temp)
      

  end function LiveFuelMoistureContent

  !---------------------------------------------------------------------------------------
  subroutine MaxHeight(height, max_height)
  !
  ! DESCRIPTION
  ! Search for max height across cohorts on patch
  !
  ! ARGUMENTS:
  real(r8),           intent(in)     :: height             ! cohort height [m]
  real(r8),           intent(inout)  :: max_height         ! max height of all cohorts [m]

  if(height > max_height)then
    max_height = height
  end if

  end subroutine MaxHeight


  !---------------------------------------------------------------------------------------


  subroutine CrownFireBehaveFM10(drying_ratio, fire_weather_index, miner_total, part_dens, wind, &
                                 canopy_bulk_density, ROS_active, CI)
    !
    ! DESCRIPTION
    ! Calculate theoretical rate of spread for a active crown fire using
    ! fuel model 10 and Rothermel's ROS model
    !
    use SFEquationsMod,    only : OptimumPackingRatio, ReactionIntensity
    use SFEquationsMod,    only : HeatofPreignition, EffectiveHeatingNumber
    use SFEquationsMod,    only : WindFactor, PropagatingFlux
    use SFEquationsMod,    only : ForwardRateOfSpread

    ! ARGUMENTS:
    real(r8), intent(in)  :: drying_ratio        ! SPITFIRE fuel parameters controlling fuel dying rate [unitless]
    real(r8), intent(in)  :: fire_weather_index  ! Nesterov fire weather index 
    real(r8), intent(in)  :: miner_total         ! SPITFIRE parameter to set fractional mineral content per unit biomass [fraction]
    real(r8), intent(in)  :: part_dens           ! SPITFIRE parameter to set particle density for fuels [kg/m2]
    real(r8), intent(in)  :: wind                ! Site wind speed [m/s]
    real(r8), intent(in)  :: canopy_bulk_density ! Canopy fuel bulk density [kg biomass / m3]               ! 
    real(r8), intent(out) :: ROS_active          ! Theoretical rate of spread a fully active crown fire using fuel model 10  [m/min]
    real(r8), intent(out) :: CI                  ! Open wind speed to sustain an active crown fire using fuel model 10 [m/min]

    ! Local variables:

    real(r8)                        :: fuel_1h           ! 1 hour fuel load using FM 10 [kg]
    real(r8)                        :: fuel_10h          ! 10 hour fuel load using FM 10 [kg]
    real(r8)                        :: fuel_100h         ! 100 hour fuel load using FM 10 [kg]
    real(r8)                        :: fuel_live         ! live fuel load using FM 10 [kg]
    real(r8)                        :: total_fuel        ! 1h + 10 h + 100h fuel using FM 10 [kg]
    real(r8)                        :: net_fuel          ! total_fuel excluding mineral content [kg] 
    real(r8)                        :: fuel_depth        ! fuel bed depth using FM 10 [m]
    real(r8)                        :: fuel_bd           ! fuel bulk density using FM 10 [kg biomass/m3]
    real(r8)                        :: fuel_sav1h        ! SAV for 1 hour fuel for FM 10 [/cm]
    real(r8)                        :: fuel_sav10h       ! SAV for 10 hour fuel for FM 10 [/cm]
    real(r8)                        :: fuel_sav100h      ! SAV for 100 hour fuel for FM 10 [/cm]
    real(r8)                        :: fuel_savlive      ! SAV for live fuel for FM 10 [/cm]
    real(r8)                        :: fuel_sav          ! mean fuel surface area to volume ratio using FM 10 [/cm]
    real(r8)                        :: fuel_eff_moist    ! fuel effective moisture content using FM 10
    real(r8)                        :: fuel_moist1h      ! FMC of 1 hour fuel for FM 10 [fraction]
    real(r8)                        :: fuel_moist10h     ! FMC of 10 hour fuel for FM 10 [fraction]
    real(r8)                        :: fuel_moist100h    ! FMC of 100 hour fuel for FM 10 [fraction]
    real(r8)                        :: fuel_moist_live   ! FMC of live fuels for FM 10 [fraction]
    real(r8)                        :: midflame_wind     ! 40% of open wind speed 
    real(r8)                        :: beta_fm10         ! packing ratio derived for fuel model 10 [unitless]
    real(r8)                        :: beta_op_fm10      ! optimum packing ratio for FM 10 [unitless]
    real(r8)                        :: beta_ratio_fm10   ! relative packing ratio for FM 10 [unitless]
    real(r8)                        :: i_r_fm10          ! reaction intensity for FM 10 [kJ/m2/min]
    real(r8)                        :: xi_fm10           ! propagating flux ratio for FM 10 [unitless]
    real(r8)                        :: eps_fm10          ! effective heating number for FM 10 [unitless]
    real(r8)                        :: phi_wind_fm10     ! wind factor for FM 10 [unitless]
    real(r8)                        :: q_ig_fm10         ! heat of pre-ignition for FM 10 [kJ/kg]

    ! Parameters for fuel model 10 to describe fuel characteristics; and some constants 
    ! fuel loading, MEF, and depth from Anderson 1982 Aids to determining fuel models for fire behavior
    ! SAV values from BEHAVE model Burgan & Rothermel (1984) 
    ! XLG: Can consider change all FM 10 parameters to SF parameters so users can use costumized fuel model for their study

    real(r8),parameter  :: fuel_1h_ton     = 3.01_r8                   ! FM 10 1-hr fuel loading (US tons/acre)
    real(r8),parameter  :: fuel_10h_ton    = 2.0_r8                    ! FM 10 10-hr fuel loading (US tons/acre)             
    real(r8),parameter  :: fuel_100h_ton   = 5.01_r8                   ! FM 10 100-hr fuel loading (US tons/acre)
    real(r8),parameter  :: fuel_live_ton = 2.0_r8                      ! FM 10 live fuel loading (US tons/acre)
    real(r8),parameter  :: fuel_mef     = 0.25_r8                      ! FM 10 moisture of extinction (volumetric), XLG: should we use this from FM 10??
    real(r8),parameter  :: fuel_depth_ft= 1.0_r8                       ! FM 10 fuel depth (ft)
    real(r8),parameter  :: sav_1h_ft   = 2000.0_r8                     ! BEHAVE model 1-hr SAV (ft2/ft3)
    real(r8),parameter  :: sav_10h_ft  = 109.0_r8                      ! BEHAVE model 10-hr SAV (ft2/ft3)             
    real(r8),parameter  :: sav_100h_ft = 30.0_r8                       ! BEHAVE model 100-hr SAV (ft2/ft3)
    real(r8),parameter  :: sav_live_ft  = 1650.0_r8                    ! BEHAVE model live SAV (ft2/ft3)
    real(r8),parameter  :: tonnes_acre_to_kg_m2 = 0.2241701_r8         ! convert tons/acre to kg/m2
    real(r8),parameter  :: sqft_cubicft_to_sqm_cubicm = 0.03280844_r8  ! convert ft2/ft3 to m2/m3
    real(r8),parameter  :: ft_to_meter = 0.3048_r8                     ! convert ft to meter
    real(r8),parameter  :: km_per_hr_to_m_per_min = 16.6667_r8         ! convert km/hour to m/min for wind speed

    fuel_1h     = fuel_1h_ton * tonnes_acre_to_kg_m2
    fuel_10h    = fuel_10h_ton * tonnes_acre_to_kg_m2
    fuel_100h   = fuel_100h_ton * tonnes_acre_to_kg_m2
    fuel_live   = fuel_live_ton * tonnes_acre_to_kg_m2
    total_fuel  = (fuel_1h + fuel_10h + fuel_100h + fuel_live) ! kg biomass /m2
    fuel_sav1h  = sav_1h_ft * sqft_cubicft_to_sqm_cubicm
    fuel_sav10h = sav_10h_ft * sqft_cubicft_to_sqm_cubicm
    fuel_sav100h = sav_100h_ft * sqft_cubicft_to_sqm_cubicm
    fuel_savlive  = sav_live_ft * sqft_cubicft_to_sqm_cubicm
    fuel_moist1h     = exp(-1.0_r8 * ((fuel_sav1h/drying_ratio) * fire_weather_index))
    fuel_moist10h    = exp(-1.0_r8 * ((fuel_sav10h/drying_ratio) * fire_weather_index))
    fuel_moist100h   = exp(-1.0_r8 * ((fuel_sav100h/drying_ratio) * fire_weather_index))
    fuel_moist_live  = exp(-1.0_r8 * ((fuel_savlive/drying_ratio) * fire_weather_index))
    fuel_depth       = fuel_depth_ft * ft_to_meter           !convert to meters
    fuel_bd          = total_fuel/fuel_depth                 !fuel bulk density (kg biomass/m3)

    fuel_sav         = fuel_sav1h *(fuel_1h/total_fuel) + fuel_sav10h*(fuel_10h/total_fuel) + & 
                       fuel_sav100h*(fuel_100h/total_fuel) + fuel_savlive*(fuel_live/total_fuel)

    fuel_eff_moist   = fuel_moist1h *(fuel_1h/total_fuel) + fuel_moist10h*(fuel_10h/total_fuel) + & 
                       fuel_moist100h*(fuel_100h/total_fuel) + fuel_moist_live*(fuel_live/total_fuel)

    net_fuel  = total_fuel * (1.0_r8 - miner_total)

    beta_fm10 = fuel_bd / part_dens
    beta_op_fm10 = OptimumPackingRatio(fuel_sav)
    if(beta_op_fm10 < nearzero) then
      beta_ratio_fm10 = 0.0_r8
    else
      beta_ratio_fm10 = beta_fm10 / beta_op_fm10
    end if

    i_r_fm10 = ReactionIntensity(net_fuel, fuel_sav, beta_ratio_fm10, &
                                 fuel_eff_moist, fuel_mef)

    q_ig_fm10 = HeatofPreignition(fuel_eff_moist) ! XLG: I'm not sure if we should use a constant eff_moist from FM 10 

    eps_fm10 = EffectiveHeatingNumber(fuel_sav)

    midflame_wind = wind * 0.40_r8 ! Scott & Reinhardt 2001 use 40% of open wind speed as effective wind speed
                                   
    phi_wind_fm10 = WindFactor(midflame_wind, beta_ratio_fm10, fuel_sav)

    xi_fm10 = PropagatingFlux(beta_fm10, fuel_sav)

    ! Calculate ROS_active, used for determining whether there is active or passtive crown fire
    ROS_active = ForwardRateOfSpread(fuel_bd, eps_fm10, q_ig_fm10, i_r_fm10, &
                                        xi_fm10, phi_wind_fm10)
   ! apply the multiplier of 3.34 for the correlation between durface ROS and active crown fire ROS 
   ! EQ. 8 
    ROS_active = ROS_active * 3.34_r8

    ! Calculate crowning index, which is used for calculating ROS_SA
    CI = CrowningIndex(eps_fm10, q_ig_fm10, i_r_fm10, &
                       canopy_bulk_density )
    CI = CI * km_per_hr_to_m_per_min  ! convert to m/min

  end subroutine CrownFireBehaveFM10

  !---------------------------------------------------------------------------------------


  subroutine BiomassBin(cbh, height, crown_depth, canopy_fuel_1h, biom_matrix)
    !
    ! DESCRIPTION:
    ! accumulate biomass at 1 m interval across cohorts and return 
    ! the biomass array
    !
    ! ARGUMENTS:
    real(r8), intent(in)                 :: cbh                 ! canopy base height [m]
    real(r8), intent(in)                 :: height              ! cohort height [m]
    real(r8), intent(in)                 :: crown_depth         ! crown length of a cohort [m]
    real(r8), intent(in)                 :: canopy_fuel_1h      ! leaf + 1-hour woody biomass of a cohort [kg biomass]
    real(r8), intent(inout)              :: biom_matrix(0:)     ! array that holds biomass by 1m interval across all cohorts on a patch [kg biomass]

    ! LOCALS:
    integer           :: h_idx                ! looping index
    real(r8)          :: crown_fuel_per_m     ! biomass by 1m interval [kg biomass]

    ! calculate biomass by 1m interval 
    crown_fuel_per_m = canopy_fuel_1h / crown_depth

    do h_idx = int(cbh), int(height)
      biom_matrix(h_idx) = biom_matrix(h_idx) + crown_fuel_per_m
    end do


  end subroutine BiomassBin

  !---------------------------------------------------------------------------------------

  subroutine CrownFireCFB(ROS_active, ROS_critical, ROS_front, &
    ROS_init, ROS_SA, active_crownfire, passive_crownfire, crown_frac_burnt)
    !
    ! DESCRIPTION:
    ! Determine if it is passive or active crown fire and
    ! return crown fraction burned EQ. 28 in Scott & Reinhardt 2001
    !
    ! ARGUMENTS:
    real(r8),  intent(in)          :: ROS_active             ! active crown fire ROS using FM10 [m/min]
    real(r8),  intent(in)          :: ROS_critical           ! critical ROS for sustaining active crown fire [m/min]
    real(r8),  intent(in)          :: ROS_front              ! surface ROS [m/min]
    real(r8),  intent(in)          :: ROS_init               ! ROS to initiate a crown fire [m/min]
    real(r8),  intent(in)          :: ROS_SA                 ! ROS at crowning index wind speed when ROS_active = ROS_SA [m/min]
    integer,   intent(out)         :: active_crownfire       ! active crown fire = 1
    integer,   intent(out)         :: passive_crownfire      ! passive crown fire = 1
    real(r8),  intent(out)         :: crown_frac_burnt       ! crown fraction burnt [0-1]
    !
    if(ROS_active >= ROS_critical)then
      active_crownfire = 1
      passive_crownfire = 0
      ! for active crown fire we set CFB to 1
      crown_frac_burnt = 1.0_r8
    else if(ROS_active < ROS_critical .and. &
    ROS_front > ROS_init .and. ROS_front < ROS_SA)then ! calculation of ROSs are disconnected so check to make sure these are true
      active_crownfire = 0
      passive_crownfire = 1
      ! calculate crown fraction burnt EQ. 28 in Scott & Reinhardt 2001
      crown_frac_burnt = min(1.0_r8, (ROS_front - ROS_init) / (ROS_SA - ROS_init))
    else
      active_crownfire = 0
      passive_crownfire = 0
      crown_frac_burnt = 0.0_r8
    end if


  end subroutine CrownFireCFB


end module CrownFireEquationsMod