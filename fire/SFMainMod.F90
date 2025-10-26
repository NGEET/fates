module SFMainMod

   ! ============================================================================
   ! All subroutines related to the SPITFIRE fire routine.
   ! Code originally developed by Allan Spessa & Rosie Fisher as part of the NERC-QUEST project.
   ! ============================================================================

   use FatesConstantsMod,      only : r8 => fates_r8
   use FatesConstantsMod,      only : itrue, ifalse
   use FatesConstantsMod,      only : pi_const
   use FatesConstantsMod,      only : nocomp_bareground, nearzero
   use FatesGlobals,           only : fates_log
   use FatesGlobals          , only : endrun => fates_endrun
   use shr_log_mod           , only : errMsg => shr_log_errMsg
   use FatesInterfaceTypesMod, only : hlm_masterproc
   use FatesInterfaceTypesMod, only : hlm_spitfire_mode
   use FatesInterfaceTypesMod, only : hlm_sf_nofire_def
   use FatesInterfaceTypesMod, only : hlm_sf_scalar_lightning_def
   use FatesInterfaceTypesMod, only : hlm_sf_successful_ignitions_def
   use FatesInterfaceTypesMod, only : hlm_sf_anthro_ignitions_def
   use FatesInterfaceTypesMod, only : bc_in_type
   use EDPftvarcon,            only : EDPftvarcon_inst
   use PRTParametersMod,       only : prt_params
   use PRTGenericMod,          only : element_pos
   use EDtypesMod,             only : ed_site_type
   use FatesPatchMod,          only : fates_patch_type
   use FatesCohortMod,         only : fates_cohort_type
   use EDtypesMod,             only : AREA
   use FatesLitterMod,         only : litter_type
   use FatesFuelClassesMod,    only : num_fuel_classes
   use PRTGenericMod,          only : carbon12_element
   use FatesInterfaceTypesMod, only : numpft
   use FatesAllometryMod,      only : CrownDepth
   use FatesFuelClassesMod,    only : fuel_classes

   implicit none
   private

   public :: DailyFireModel
   public :: UpdateFuelCharacteristics

   ! ======================================================================================

contains

   subroutine DailyFireModel(currentSite, bc_in)
      !
      !  DESCRIPTION:
      !  Runs the daily fire model

      ! ARGUMENTS:
      type(ed_site_type), intent(inout), target :: currentSite ! site object
      type(bc_in_type),   intent(in)            :: bc_in       ! BC in object

      ! LOCALS:
      type (fates_patch_type), pointer :: currentPatch ! patch object

      if (hlm_spitfire_mode > hlm_sf_nofire_def) then
         call UpdateFireWeather(currentSite, bc_in)
         call UpdateFuelCharacteristics(currentSite)
         call CalculateIgnitionsandFDI(currentSite, bc_in)
         call CalculateSurfaceRateOfSpread(currentSite)
         call CalculateSurfaceFireIntensity(currentSite)
         call CalculateAreaBurnt(currentSite)
         call CalculateRxFireAreaBurnt(currentSite)
         call CalculatePostFireMortality(currentSite)
      end if

   end subroutine DailyFireModel

   !---------------------------------------------------------------------------------------

   subroutine UpdateFireWeather(currentSite, bc_in)
      !
      !  DESCRIPTION:
      !  Updates the site's fire weather index, burn window for prescribed fire, and calculates effective windspeed based on
      !   vegetation characteristics
      !
      !  Currently we use tree and grass fraction averaged over whole grid (site) to
      !  prevent extreme divergence

      use FatesConstantsMod,      only : tfrz => t_water_freeze_k_1atm
      use FatesConstantsMod,      only : sec_per_day, sec_per_min
      use EDTypesMod,             only : CalculateTreeGrassAreaSite
      use FatesInterfaceTypesMod, only : hlm_use_managed_fire
      use SFParamsMod,            only : SF_val_rxfire_tpup, SF_val_rxfire_tplw, SF_val_rxfire_rhup, &
         SF_val_rxfire_rhlw, SF_val_rxfire_wdup, SF_val_rxfire_wdlw

      ! ARGUMENTS:
      type(ed_site_type), intent(inout), target :: currentSite
      type(bc_in_type),   intent(in)            :: bc_in

      ! LOCALS:
      type(fates_patch_type), pointer :: currentPatch   ! patch object
      real(r8)                        :: temp_C         ! daily averaged temperature [deg C]
      real(r8)                        :: precip         ! daily precip [mm/day]
      real(r8)                        :: rh             ! daily relative humidity [%]
      real(r8)                        :: wind           ! wind speed [m/s]
      real(r8)                        :: tree_fraction  ! site-level tree fraction [0-1]
      real(r8)                        :: grass_fraction ! site-level grass fraction [0-1]
      real(r8)                        :: bare_fraction  ! site-level bare ground fraction [0-1]
      integer                         :: iofp           ! index of oldest the fates patch

      ! NOTE that the boundary conditions of temperature, precipitation and relative humidity
      ! are available at the patch level. We are currently using a simplification where the whole site
      ! is simply using the values associated with the first patch.
      ! which probably won't have much impact, unless we decide to ever calculated fire weather for each patch.

      currentPatch => currentSite%oldest_patch

      ! If the oldest patch is a bareground patch (i.e. nocomp mode is on) use the first vegetated patch
      ! for the iofp index (i.e. the next younger patch)
      if (currentPatch%nocomp_pft_label == nocomp_bareground) then
         currentPatch => currentPatch%younger
      endif

      iofp = currentPatch%patchno
      temp_C = currentPatch%tveg24%GetMean() - tfrz
      precip = bc_in%precip24_pa(iofp)*sec_per_day
      rh = bc_in%relhumid24_pa(iofp)
      wind = bc_in%wind24_pa(iofp)

      ! convert to m/min
      currentSite%wind = wind*sec_per_min

      ! update fire weather index
      call currentSite%fireWeather%UpdateIndex(temp_C, precip, rh, wind)

      ! update prescribed fire burn window
      call currentSite%fireWeather%UpdateRxfireBurnWindow(hlm_use_managed_fire, temp_C, rh, wind, &
         SF_val_rxfire_tpup, SF_val_rxfire_tplw, SF_val_rxfire_rhup, SF_val_rxfire_rhlw,    &
         SF_val_rxfire_wdup, SF_val_rxfire_wdlw)

      ! calculate site-level tree, grass, and bare fraction
      call CalculateTreeGrassAreaSite(currentSite, tree_fraction, grass_fraction, bare_fraction)

      ! update effective wind speed
      call currentSite%fireWeather%UpdateEffectiveWindSpeed(wind*sec_per_min, tree_fraction, &
         grass_fraction, bare_fraction)

   end subroutine UpdateFireWeather

   !---------------------------------------------------------------------------------------

   subroutine UpdateFuelCharacteristics(currentSite)
      !
      !  DESCRIPTION:
      !  Updates fuel characteristics on each patch of the site
      !

      use SFParamsMod, only : SF_val_drying_ratio, SF_val_SAV, SF_val_FBD

      ! ARGUMENTS:
      type(ed_site_type), intent(in), target :: currentSite  ! site object

      ! LOCALS:
      type(fates_patch_type), pointer :: currentPatch ! FATES patch
      type(litter_type),      pointer :: litter       ! pointer to patch litter class
      real(r8) :: MEF_trunks, fuel_moisture_trunks

      currentPatch => currentSite%oldest_patch
      do while(associated(currentPatch))

         if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

            ! calculate live grass [kgC/m2]
            call currentPatch%UpdateLiveGrass()

            ! update fuel loading [kgC/m2]
            litter => currentPatch%litter(element_pos(carbon12_element))
            call currentPatch%fuel%UpdateLoading(sum(litter%leaf_fines(:)),                  &
               litter%ag_cwd(1), litter%ag_cwd(2), litter%ag_cwd(3), litter%ag_cwd(4),        &
               currentPatch%livegrass)
            ! update weighting factors
            call currentPatch%fuel%CalculateWeightingFactor(SF_val_SAV, SF_val_part_dens)

            ! sum up fuel classes and calculate fractional loading for each
            call currentPatch%fuel%SumLoading(SF_val_SAV, SF_val_part_dens)
            call currentPatch%fuel%CalculateFractionalLoading(SF_val_SAV, SF_val_part_dens)

            ! calculate fuel moisture [m3/m3]
            call currentPatch%fuel%UpdateFuelMoisture(SF_val_SAV, SF_val_drying_ratio,       &
               currentSite%fireWeather)

            ! calculate geometric properties
            call currentPatch%fuel%AverageBulkDensity(SF_val_FBD)
            call currentPatch%fuel%AverageSAV(SF_val_SAV)

         end if
         currentPatch => currentPatch%younger
      end do

   end subroutine UpdateFuelCharacteristics

   !---------------------------------------------------------------------------------------

   subroutine CalculateIgnitionsandFDI(currentSite, bc_in)
      !
      !  DESCRIPTION:
      !  Calculates ignitions and fire danger index (FDI) for a site
      !

      use FatesInterfaceTypesMod, only : hlm_spitfire_mode
      use EDParamsMod,            only : cg_strikes
      use EDParamsMod,            only : ED_val_nignitions
      use SFParamsMod,            only : SF_val_fdi_alpha
      use FatesConstantsMod,      only : years_per_day

      ! ARGUMENTS:
      type(ed_site_type), intent(inout), target :: currentSite ! site object
      type(bc_in_type),   intent(in)            :: bc_in       ! BC in object

      ! LOCALS:
      type(fates_patch_type), pointer :: currentPatch            ! patch object
      real(r8)                        :: cloud_to_ground_strikes ! fraction of cloud-to-ground strikes [0-1]
      real(r8)                        :: anthro_ignitions        ! anthropogenic ignitions [count/km2/day]
      integer                         :: iofp                    ! patch index

      ! CONSTANTS:
      real(r8), parameter :: igns_per_person_month = 0.0035_r8  ! potential human ignition counts (alpha in Li et al. 2012) (#/person/month)
      real(r8), parameter :: approx_days_per_month = 30.0_r8    ! approximate days per month [days]

      ! initialize site parameters to zero
      currentSite%NF_successful = 0.0_r8

      ! Equation 7 from Venevsky et al GCB 2002 (modification of equation 8 in Thonicke et al. 2010)
      ! FDI 0.1 = low, 0.3 moderate, 0.75 high, and 1 = extreme ignition potential for alpha 0.000337
      if (hlm_spitfire_mode == hlm_sf_successful_ignitions_def) then
         ! READING "SUCCESSFUL IGNITION" DATA
         ! force ignition potential to be extreme
         ! cloud_to_ground_strikes = 1 means using 100% of incoming observed ignitions
         currentSite%FDI = 1.0_r8
         cloud_to_ground_strikes = 1.0_r8
      else
         ! USING LIGHTNING STRIKE DATA
         currentSite%FDI  = 1.0_r8 - exp(-SF_val_fdi_alpha*currentSite%fireWeather%fire_weather_index)
         cloud_to_ground_strikes = cg_strikes
      end if

      ! if the oldest patch is a bareground patch (i.e. nocomp mode is on) use the first vegetated patch
      ! for the iofp index (i.e. the next younger patch)
      currentPatch => currentSite%oldest_patch
      if (currentPatch%nocomp_pft_label == nocomp_bareground)then
         currentPatch => currentPatch%younger
      endif
      iofp = currentPatch%patchno

      ! NF = number of lighting strikes per day per km2 scaled by cloud to ground strikes
      if (hlm_spitfire_mode == hlm_sf_scalar_lightning_def) then
         currentSite%NF = ED_val_nignitions*years_per_day*cloud_to_ground_strikes
      else
         ! use external daily lightning ignition data
         currentSite%NF = bc_in%lightning24(iofp)*cloud_to_ground_strikes
      end if

      ! calculate anthropogenic ignitions according to Li et al. (2012)
      ! add to ignitions by lightning
      if (hlm_spitfire_mode == hlm_sf_anthro_ignitions_def) then
         ! anthropogenic ignitions (count/km2/day)
         !           =  (ignitions/person/month)*6.8*population_density**0.43/approximate days per month
         anthro_ignitions = igns_per_person_month*6.8_r8*bc_in%pop_density(iofp)**0.43_r8/approx_days_per_month
         currentSite%NF = currentSite%NF + anthro_ignitions
      end if

   end subroutine CalculateIgnitionsandFDI

   !---------------------------------------------------------------------------------------

   subroutine CalculateSurfaceRateOfSpread(currentSite)
      !
      !  DESCRIPTION:
      !  Calculates potential rate of spread based on fuel characteristics for
      !  each patch of a site
      !

      use SFParamsMod,    only : SF_val_miner_total, SF_val_part_dens, SF_val_SAV
      use SFEquationsMod, only : OptimumPackingRatio, ReactionIntensity
      use SFEquationsMod, only : HeatofPreignition, EffectiveHeatingNumber
      use SFEquationsMod, only : WindFactor, HeatSink, PropagatingFlux
      use SFEquationsMod, only : ForwardRateOfSpread, BackwardRateOfSpread

      ! ARGUMENTS:
      type(ed_site_type), intent(in), target :: currentSite ! site object

      ! LOCALS:
      type(fates_patch_type), pointer :: currentPatch ! patch object
      real(r8)                        :: beta         ! packing ratio [unitless]
      real(r8)                        :: beta_op      ! optimum packing ratio [unitless]
      real(r8)                        :: beta_ratio   ! relative packing ratio [unitless]
      real(r8)                        :: IR_dead      ! reaction intensity of dead fuels [kJ m-2 min-1]
      real(r8)                        :: IR_live      ! reaction intensity of live fuels [kJ m-2 min-1]
      real(r8)                        :: i_r          ! summed reaction intensity of dead and live fuels [kJ/m2/min]
      real(r8)                        :: xi           ! propagating flux ratio [unitless]
      real(r8)                        :: heat_sink    ! total heat required to ignite per unit volume fuel bed [kJ m-3]
      real(r8)                        :: eps(num_fuel_classes)  ! effective heating number by fuel class [unitless]
      real(r8)                        :: phi_wind               ! wind factor [unitless]
      real(r8)                        :: q_ig(num_fuel_classes) ! heat of pre-ignition by fuel class [kJ/kg]
      integer                         :: i            ! looping index

      currentPatch => currentSite%oldest_patch
      do while(associated(currentPatch))
         if (currentPatch%nocomp_pft_label /= nocomp_bareground .and.                       &
            currentPatch%fuel%non_trunk_loading > nearzero) then

            ! fraction of fuel array volume occupied by fuel, i.e. compactness of fuel bed [unitless]
            ! Rothermel 1972 Eq. 31
            beta = currentPatch%fuel%bulk_density_weighted/SF_val_part_dens

            ! optimum packing ratio [unitless]
            beta_op = OptimumPackingRatio(currentPatch%fuel%SAV_weighted)

            ! relative packing ratio [unitless]
            if (beta_op < nearzero) then
               beta_ratio = 0.0_r8
            else
               beta_ratio = beta/beta_op
            end if

            ! remove mineral content from fuel load per Thonicke 2010
            currentPatch%fuel%non_trunk_loading = currentPatch%fuel%non_trunk_loading*(1.0_r8 - SF_val_miner_total)
            currentPatch%fuel%weighted_loading_dead = currentPatch%fuel%weighted_loading_dead*(1.0_r8 - SF_val_miner_total)
            currentPatch%fuel%weighted_loading_live = currentPatch%fuel%weighted_loading_live*(1.0_r8 - SF_val_miner_total)

            ! reaction intensity [kJ/m2/min]
            IR_dead = ReactionIntensity(currentPatch%fuel%weighted_loading_dead/0.45_r8,             &
               currentPatch%fuel%SAV_weighted, beta_ratio,                                    &
               currentPatch%fuel%average_moisture_dead, currentPatch%fuel%MEF_dead)

            IR_live = ReactionIntensity(currentPatch%fuel%weighted_loading_live/0.45_r8,             &
               currentPatch%fuel%SAV_weighted, beta_ratio,                                    &
               currentPatch%fuel%average_moisture_live, currentPatch%fuel%MEF_live)
            i_r = IR_dead + IR_live

            do i = 1, num_fuel_classes
               ! heat of preignition per fuel class [kJ/kg]
               q_ig(i) = HeatofPreignition(currentPatch%fuel%moisture(i))

               ! effective heating number per fuel class [unitless]
               eps(i) = EffectiveHeatingNumber(SF_val_SAV(i))
            end do
            ! total heat required to ignite per unit fuel bed
            heat_sink = HeatSink(q_ig, eps, currentPatch%fuel%weighting_factor, &
               currentPatch%fuel%bulk_density_weighted, currentPatch%fuel%wf_dead, &
               currentPatch%fuel%wf_live)

            ! wind factor [unitless]
            phi_wind = WindFactor(currentSite%fireWeather%effective_windspeed, beta_ratio,      &
               currentPatch%fuel%SAV_weighted)

            ! propagating flux [unitless]
            xi = PropagatingFlux(beta, currentPatch%fuel%SAV_weighted)

            ! forward rate of spread [m/min]
            currentPatch%ROS_front = ForwardRateOfSpread(heat_sink, i_r, xi, phi_wind)

            ! backwards rate of spread [m/min]
            !  backward ROS wind not changed by vegetation - so use wind, not effective_windspeed
            currentPatch%ROS_back = BackwardRateOfSpread(currentPatch%ROS_front,             &
               currentSite%wind)

         end if
         currentPatch => currentPatch%younger
      end do

   end subroutine CalculateSurfaceRateOfSpread

   !---------------------------------------------------------------------------------------

   subroutine CalculateSurfaceFireIntensity(currentSite)
      !
      !  DESCRIPTION:
      !  Calculates surface fireline intensity for each patch of a site
      !  Use calculated fire intensity to determine if prescribed fire or
      !  wildfire happens
      !

      use SFEquationsMod,    only : FireIntensity
      use SFParamsMod,       only : SF_val_fire_threshold
      use SFParamsMod,       only : SF_val_rxfire_max_threshold, SF_val_rxfire_min_threshold
      use SFParamsMod,       only : SF_val_rxfire_fuel_max, SF_val_rxfire_fuel_min
      use FatesRxFireMod,    only : is_prescribed_burn

      ! ARGUMENTS:
      type(ed_site_type), intent(inout), target :: currentSite

      ! LOCALS:
      type(fates_patch_type), pointer :: currentPatch                    ! patch object
      real(r8)                        :: fuel_consumed(num_fuel_classes) ! fuel consumed [kgC/m2]
      logical                         :: is_rxfire                       ! is it a prescribed fire?
      logical                         :: rxfire_fuel_check               ! is fuel within thresholds for prescribed burn
      logical                         :: fi_check  ! is (potential) fire intensity high enough for fire to actually happen?
      logical                         :: has_ignition                    ! is ignition greater than zero?

      currentPatch => currentSite%oldest_patch
      do while (associated(currentPatch))

         currentPatch%fuel%frac_burnt(:) = 0.0_r8

         if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

            call currentPatch%fuel%CalculateFuelBurnt(fuel_consumed)
            call currentPatch%fuel%CalculateResidenceTime(currentPatch%tau_l)

            ! calculate overall fuel consumed by spreading fire
            ! ignore 1000-hr fuels (i.e. trunks)
            currentPatch%TFC_ROS = sum(fuel_consumed) - fuel_consumed(fuel_classes%trunks())

            ! initialize patch parameters to zero
            currentPatch%FI = 0.0_r8        ! either nonrx or rx FI
            currentPatch%nonrx_fire = 0     ! only wildfire
            currentPatch%rx_fire = 0        ! only rx fire
            currentPatch%rx_FI = 0.0_r8
            currentPatch%nonrx_FI = 0.0_r8

            has_ignition = currentSite%NF > 0.0_r8

            if (has_ignition .or. currentSite%fireWeather%rx_flag == itrue) then

               ! fire intensity [kW/m]
               currentPatch%FI = FireIntensity(currentPatch%TFC_ROS/0.45_r8, currentPatch%ROS_front/60.0_r8)
               fi_check = currentPatch%FI > SF_val_fire_threshold

               ! check if prescribed fire can occur based on fuel load
               rxfire_fuel_check = currentPatch%fuel%non_trunk_loading > SF_val_rxfire_fuel_min .and. &
                  currentPatch%fuel%non_trunk_loading < SF_val_rxfire_fuel_max

               if (currentSite%fireWeather%rx_flag == itrue .and. rxfire_fuel_check) then

                  ! record burnable area after fuel load check
                  currentSite%rxfire_area_fuel = currentSite%rxfire_area_fuel + currentPatch%area

                  ! determine fire type
                  ! prescribed fire and wildfire cannot happen on the same patch
                  is_rxfire = is_prescribed_burn(currentPatch%FI, currentSite%NF, &
                     SF_val_rxfire_min_threshold, SF_val_rxfire_max_threshold, SF_val_fire_threshold)

                  if (is_rxfire) then
                     currentSite%rxfire_area_fi = currentSite%rxfire_area_fi + currentPatch%area ! record burnable area after FI check
                     currentPatch%rx_fire = 1

                  else if (has_ignition .and. fi_check) then  ! (potential) intensity is greater than kW/m energy threshold
                     currentPatch%nonrx_fire = 1
                  end if

               else if (has_ignition .and. fi_check)  then ! not a patch suitable for conducting prescribed fire or rxfire is not even turned on, but (potential) intensity is greater than kW/m energy threshold
                  currentPatch%nonrx_fire = 1
               end if

               ! assign fire intensities and ignitions based on fire type
               if (currentPatch%nonrx_fire == itrue) then
                  currentSite%NF_successful = currentSite%NF_successful + &
                     currentSite%NF*currentSite%FDI*currentPatch%area/area
                  currentPatch%nonrx_FI = currentPatch%FI
               else if (currentPatch%rx_fire == itrue) then
                  currentPatch%rx_FI = currentPatch%FI
               end if
            end if
         end if
         currentPatch => currentPatch%younger
      end do

   end subroutine CalculateSurfaceFireIntensity

   !---------------------------------------------------------------------------------------

   subroutine CalculateAreaBurnt(currentSite)
      !
      !  DESCRIPTION:
      !  Calculates area burnt for each patch of a site
      !
      use FatesConstantsMod, only : m2_per_km2
      use SFEquationsMod,    only : FireDuration, LengthToBreadth
      use SFEquationsMod,    only : AreaBurnt, FireSize

      ! ARGUMENTS:
      type(ed_site_type), intent(inout), target :: currentSite

      ! LOCALS:
      type(fates_patch_type), pointer :: currentPatch                    ! patch object
      real(r8)                        :: tree_fraction_patch             ! treed fraction on patch [0-1]
      real(r8)                        :: length_to_breadth               ! length to breadth ratio of fire ellipse (unitless)
      real(r8)                        :: fire_size                       ! size of fire [m2]
      real(r8)                        :: area_burnt                      ! area burnt [m2/km2]

      ! CONSTANTS:
      real(r8), parameter :: max_frac_burnt = 0.99_r8 ! maximum fraction burnt on patch

      currentPatch => currentSite%oldest_patch
      do while (associated(currentPatch))

         if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

            ! initialize patch parameters to zero
            currentPatch%FD = 0.0_r8
            currentPatch%nonrx_frac_burnt = 0.0_r8

            if (currentPatch%nonrx_fire == 1) then

               ! fire duration [min]
               currentPatch%FD = FireDuration(currentSite%FDI)

               ! length-to-breadth ratio of fire ellipse [unitless]
               tree_fraction_patch  = currentPatch%total_tree_area/currentPatch%area
               length_to_breadth = LengthToBreadth(currentSite%fireWeather%effective_windspeed, tree_fraction_patch)

               ! fire size [m2]
               fire_size = FireSize(length_to_breadth, currentPatch%ROS_back, &
                  currentPatch%ROS_front, currentPatch%FD)

               ! area burnt [m2/km2]
               area_burnt = AreaBurnt(fire_size, currentSite%NF, currentSite%FDI)

               ! convert to area burned per area patch per day
               ! i.e., fraction of the patch burned on that day
               currentPatch%nonrx_frac_burnt = min(max_frac_burnt, area_burnt/m2_per_km2)

            end if
         end if
         currentPatch => currentPatch%younger
      end do

   end subroutine CalculateAreaBurnt

   !---------------------------------------------------------------------------------------

   subroutine CalculateRxFireAreaBurnt (currentSite)
      !
      !  DESCRIPTION:
      !  Returns burned fraction for prescribed fire per patch by first checking
      !  if total burnable fraction at site level is greater than user defined fraction of site area
      !  if yes, calculate burned fraction as (user defined frac / total burnable frac)
      !
      use SFParamsMod, only : SF_val_rxfire_AB       ! user defined prescribed fire area in fraction per day to reflect burning capacity
      use SFParamsMod, only : SF_val_rxfire_min_frac ! minimum fraction of land needs to be burnable for conducting prescribed fire

      ! ARGUMENTS
      type(ed_site_type), intent(inout), target :: currentSite

      ! LOCALS
      type(fates_patch_type), pointer :: currentPatch
      real(r8)                        :: total_burnable_frac ! total fractional land area that can apply prescribed fire after condition checks at site level

      ! initialize site variables
      currentSite%rxfire_area_final = 0.0_r8
      total_burnable_frac = 0.0_r8

      ! update total burnable fraction
      total_burnable_frac = currentSite%rxfire_area_fi/AREA

      currentPatch => currentSite%oldest_patch

      do while (associated(currentPatch))
         if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
            currentPatch%fire = 0 ! fire, either rx or non-rx
            currentPatch%frac_burnt = 0.0_r8 ! rx_frac_burnt + nonrx_frac_burnt
            currentPatch%rx_frac_burnt = 0.0_r8
            if (currentPatch%rx_fire == itrue .and. &
               total_burnable_frac >= SF_val_rxfire_min_frac ) then
               currentSite%rxfire_area_final = currentSite%rxfire_area_final + currentPatch%area ! the final burned total land area
               currentPatch%rx_frac_burnt = min(0.99_r8, SF_val_rxfire_AB/total_burnable_frac)
            else
               currentPatch%rx_fire = 0 ! update rxfire occurence at patch
               currentPatch%rx_FI = 0.0_r8
            end if

            ! update patch level fire occurence and total frac burnt
            currentPatch%fire = currentPatch%nonrx_fire + currentPatch%rx_fire
            currentPatch%frac_burnt = currentPatch%nonrx_frac_burnt + currentPatch%rx_frac_burnt

            ! currentPatch%fire cannot be >1, which indicates both rx and wildfire are happening
            ! we currently do not allow this to happen on the same patch yet
            if (currentPatch%fire > 1) then
               write(fates_log(),*) 'Both wildfire and management fire are happening at same patch'
               write(fates_log(),*) 'rxfire =', currentPatch%rx_fire
               write(fates_log(),*) 'wildfire =', currentPatch%nonrx_fire
               call endrun(msg=errMsg(__FILE__, __LINE__))
            end if
         end if
         currentPatch => currentPatch%younger
      end do

   end subroutine CalculateRxFireAreaBurnt

   !---------------------------------------------------------------------------------------

   subroutine CalculatePostFireMortality(currentSite)
      !
      !  DESCRIPTION:
      !  Calculates mortality (for woody PFTs) due to fire from crown scorching and cambial damage
      !
      use SFEquationsMod, only : ScorchHeight, CrownFireMortality, CrownFractionBurnt
      use SFEquationsMod, only : CambialMortality, TotalFireMortality

      ! ARGUMENTS:
      type(ed_site_type), intent(in), target :: currentSite ! site object

      ! LOCALS:
      type(fates_patch_type),  pointer :: currentPatch    ! patch object
      type(fates_cohort_type), pointer :: currentCohort   ! cohort object
      real(r8)                         :: crown_depth     ! crown depth [m]
      integer                          :: i_pft           ! looping index

      currentPatch => currentSite%oldest_patch
      do while (associated(currentPatch))
         if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
            if (currentPatch%fire == 1) then

               ! calculate scorch height [m]
               do i_pft = 1, numpft
                  if (prt_params%woody(i_pft) == itrue) then
                     currentPatch%Scorch_ht(i_pft) = ScorchHeight(EDPftvarcon_inst%fire_alpha_SH(i_pft), &
                        currentPatch%FI)
                  else
                     currentPatch%Scorch_ht(i_pft) = 0.0_r8
                  end if
               end do

               ! calculate fire-related mortality
               currentCohort => currentPatch%tallest
               do while (associated(currentCohort))

                  currentCohort%fraction_crown_burned = 0.0_r8
                  currentCohort%fire_mort = 0.0_r8
                  currentCohort%crownfire_mort = 0.0_r8
                  currentCohort%cambial_mort = 0.0_r8

                  if (prt_params%woody(currentCohort%pft) == itrue) then

                     ! calculate crown fraction burned [0-1]
                     call CrownDepth(currentCohort%height, currentCohort%pft, crown_depth)
                     currentCohort%fraction_crown_burned = CrownFractionBurnt(currentPatch%Scorch_ht(currentCohort%pft),  &
                        currentCohort%height, crown_depth)

                     ! shrink canopy to account for burnt section
                     ! currentCohort%canopy_trim = min(currentCohort%canopy_trim, 1.0_r8 - currentCohort%fraction_crown_burned)

                     ! calculate cambial mortality rate [0-1]
                     currentCohort%cambial_mort = CambialMortality(EDPftvarcon_inst%bark_scaler(currentCohort%pft), &
                        currentCohort%dbh, currentPatch%tau_l)

                     ! calculate crown fire mortality [0-1]
                     currentCohort%crownfire_mort = CrownFireMortality(EDPftvarcon_inst%crown_kill(currentCohort%pft), &
                        currentCohort%fraction_crown_burned)

                     ! total fire mortality [0-1]
                     currentCohort%fire_mort = TotalFireMortality(currentCohort%crownfire_mort, &
                        currentCohort%cambial_mort)

                  end if
                  currentCohort => currentCohort%shorter
               end do
            end if
         end if
         currentPatch => currentPatch%younger
      end do

   end subroutine CalculatePostFireMortality

   !---------------------------------------------------------------------------------------

end module SFMainMod
