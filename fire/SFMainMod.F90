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
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : struct_organ
  use FatesInterfaceTypesMod, only : numpft
  use FatesAllometryMod,      only : CrownDepth
  use FatesFuelClassesMod,    only : fuel_classes
  
  implicit none
  private


  public :: DailyFireModel
  public :: UpdateFuelCharacteristics

  integer :: write_SF = ifalse   ! for debugging
  logical :: debug = .false.     ! for debugging

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
      call UpdateCanopyFuelCharacteristics(currentSite)
      call CalculateIgnitionsandFDI(currentSite, bc_in)
      call CalculateSurfaceRateOfSpread(currentSite)
      call CalculateSurfaceFireIntensity(currentSite)
      call PassiveActiveCrownFireCheck(currentSite)
      call CalculateAreaBurnt(currentSite)
      call CalculateRxfireAreaBurnt(currentSite)
      call crown_scorching(currentSite)
      call crown_damage(currentSite)
      call cambial_damage_kill(currentSite)
      call post_fire_mortality(currentSite)
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

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : sec_per_day, sec_per_min
    use EDTypesMod,        only : CalculateTreeGrassAreaSite
    use EDParamsMod,       only : rxfire_switch
    use SFParamsMod,       only : SF_val_rxfire_tpup, SF_val_rxfire_tplw, SF_val_rxfire_rhup, &
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
    call currentSite%fireWeather%UpdateRxfireBurnWindow(rxfire_switch, temp_C, rh, wind, &
    SF_val_rxfire_tpup, SF_val_rxfire_tplw, SF_val_rxfire_rhup, SF_val_rxfire_rhlw, &
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

    use SFParamsMod,     only : SF_val_drying_ratio, SF_val_SAV, SF_val_FBD

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
            
        ! sum up fuel classes and calculate fractional loading for each
        call currentPatch%fuel%SumLoading()
        call currentPatch%fuel%CalculateFractionalLoading()
          
        ! calculate fuel moisture [m3/m3]
        call currentPatch%fuel%UpdateFuelMoisture(SF_val_SAV, SF_val_drying_ratio,       &
          currentSite%fireWeather)
        
        ! calculate geometric properties
        call currentPatch%fuel%AverageBulkDensity_NoTrunks(SF_val_FBD)
        call currentPatch%fuel%AverageSAV_NoTrunks(SF_val_SAV)
            
      end if 
      currentPatch => currentPatch%younger
    end do 

  end subroutine UpdateFuelCharacteristics

  !---------------------------------------------------------------------------------------

  subroutine UpdateCanopyFuelCharacteristics(currentSite)
    !
    ! DESCRIPTION:
    ! Calculate canopy fuel load (sum of all 1 hour live fuels, in kg biomass),
    ! canopy base height (minimum canopy height at which canopy fuel bed density > 0.011 kg/m3 for 
    ! vertical fire spread through canopy), and canopy bulk desity (kg biomass/m3)
    !
    use FatesLitterMod,  only : ncwd
    use SFParamsMod,     only : SF_val_CWD_frac
    use FatesLitterMod,  only : adjust_SF_CWD_frac
  

    ! ARGUMENTS:
    type(ed_site_type), intent(in), target :: currentSite  ! site object

   
    type(fates_patch_type), pointer  :: currentPatch   ! FATES patch
    type(fates_cohort_type), pointer :: currentCohort  ! FATES cohort

    ! Locals:
    real(r8) ::  crown_depth          ! depth of crown [m]
    real(r8) ::  cbh_co               ! canopy base height of cohort [m]
    real(r8) ::  max_height           ! max cohort height on patch (m)
    real(r8) ::  woody_c              ! above-ground tree struct and sapwood biomass in cohort (kgC)
    real(r8) ::  leaf_c               ! leaf carbon (kgC)
    real(r8) ::  sapw_c               ! sapwood carbon (kgC)
    real(r8) ::  struct_c             ! structure carbon (kgC)
    real(r8) ::  crown_fuel_per_m     ! crown fuel per 1m section in cohort
    real(r8) ::  SF_val_CWD_frac_adj(ncwd)  ! adjusted fractional allocation of woody biomass to coarse wood debris pool

    real(r8), dimension(:), allocatable :: biom_matrix   ! matrix to track biomass from bottom to top
    integer :: h_idx                                     ! index 

    real(r8), parameter :: carbon_2_biomass = 0.45_r8

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch))
      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then
        !zero Patch level variables
        max_height                          = 0.0_r8
      
        ! find the max cohort height to set the upper bounds of biom_matrix
        currentCohort=>currentPatch%tallest
        do while(associated(currentCohort))
          if ( int(prt_params%woody(currentCohort%pft)) == itrue) then !trees
            if (currentCohort%height > max_height) then
              max_height = currentCohort%height
            end if
          end if ! trees only
          currentCohort => currentCohort%shorter;
        end do ! end cohort loop
  
        !allocate and initialize biom_matrix
        allocate(biom_matrix(0:int(max_height)))
        biom_matrix(:) = 0.0_r8
  
        !loop across cohorts to calculate canopy fuel load by 1m height bin
        currentCohort => currentPatch%tallest
        
        do while(associated(currentCohort))
          !zero cohort level variables
          woody_c            = 0.0_r8
          leaf_c             = 0.0_r8
          sapw_c             = 0.0_r8
          struct_c           = 0.0_r8
          crown_fuel_per_m   = 0.0_r8
          crown_depth        = 0.0_r8
          cbh_co             = 0.0_r8
          SF_val_CWD_frac_adj(ncwd) = 0.0_r8
  
          ! Calculate crown 1hr fuel biomass (leaf, twig sapwood, twig structural biomass)
          if ( int(prt_params%woody(currentCohort%pft)) == itrue) then !trees
            ! calculate canopy base height for cohort
            call CrownDepth(currentCohort%height,currentCohort%pft,crown_depth)
            cbh_co = currentCohort%height - crown_depth
        
            leaf_c   = currentCohort%prt%GetState(leaf_organ, carbon12_element)
            sapw_c   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
            struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
  
            woody_c  = currentCohort%n * (prt_params%allom_agb_frac(currentCohort%pft)* &
            (sapw_c + struct_c))
            leaf_c   = currentCohort%n * leaf_c
            
            call adjust_SF_CWD_frac(currentCohort%dbh, ncwd, SF_val_CWD_frac, SF_val_CWD_frac_adj)
            ! update canopy fuel load
            call currentPatch%fuel%CalculateCanopyFuelLoad(leaf_c, woody_c, SF_val_CWD_frac_adj)
            
            ! 1m biomass bin
            crown_fuel_per_m = (leaf_c + woody_c) / (carbon_2_biomass * crown_depth) !kg biomass / m
            ! sort crown fuel into bins from bottom to top of crown
            ! accumulate across cohorts to find density within canopy 1m sections
                 do h_idx = int(cbh_co), int(currentCohort%height)
                    biom_matrix(h_idx) = biom_matrix(h_idx) + crown_fuel_per_m
                 end do

          end if ! trees only
          currentCohort => currentCohort%shorter;
        end do ! end cohort loop

        biom_matrix(:) = biom_matrix(:) / currentPatch%area ! kg biomass / m3
        ! update canopy fuel bulk density
        call currentPatch%fuel%CalculateCanopyBulkDensity(biom_matrix, max_height)

        deallocate(biom_matrix)

      end if ! nocomp_bareground 
      currentPatch => currentPatch%younger;
    end do ! end patch loop

  end subroutine UpdateCanopyFuelCharacteristics

     
      

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
    real(r8), parameter :: alpha = 0.0035_r8  ! potential human ignition counts (alpha in Li et al. 2012) (#/person/month)

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
    if(currentPatch%nocomp_pft_label .eq. nocomp_bareground)then
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
      anthro_ignitions = alpha*6.8_r8*bc_in%pop_density(iofp)**0.43_r8/30.0_r8
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

    use SFParamsMod,    only : SF_val_miner_total, SF_val_part_dens
    use SFEquationsMod, only : OptimumPackingRatio, ReactionIntensity
    use SFEquationsMod, only : HeatofPreignition, EffectiveHeatingNumber
    use SFEquationsMod, only : WindFactor, PropagatingFlux
    use SFEquationsMod, only : ForwardRateOfSpread, BackwardRateOfSpread

    ! ARGUMENTS:
    type(ed_site_type), intent(in), target :: currentSite ! site object

    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch ! patch object 
    real(r8)                        :: beta         ! packing ratio [unitless]
    real(r8)                        :: beta_op      ! optimum packing ratio [unitless]
    real(r8)                        :: beta_ratio   ! relative packing ratio [unitless]
    real(r8)                        :: i_r          ! reaction intensity [kJ/m2/min]
    real(r8)                        :: xi           ! propagating flux ratio [unitless]
    real(r8)                        :: eps          ! effective heating number [unitless]
    real(r8)                        :: phi_wind     ! wind factor [unitless]
    real(r8)                        :: q_ig         ! heat of pre-ignition [kJ/kg]

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
      if (currentPatch%nocomp_pft_label /= nocomp_bareground .and.                       &
        currentPatch%fuel%non_trunk_loading > nearzero) then
        
        ! fraction of fuel array volume occupied by fuel, i.e. compactness of fuel bed [unitless]
        ! Rothermel 1972 Eq. 31
        beta = currentPatch%fuel%bulk_density_notrunks/SF_val_part_dens
        
        ! optimum packing ratio [unitless]
        beta_op = OptimumPackingRatio(currentPatch%fuel%SAV_notrunks)
        
        ! relative packing ratio [unitless]
        if (beta_op < nearzero) then 
          beta_ratio = 0.0_r8
        else
          beta_ratio = beta/beta_op 
        end if
        
        ! remove mineral content from fuel load per Thonicke 2010 
        currentPatch%fuel%non_trunk_loading = currentPatch%fuel%non_trunk_loading*(1.0_r8 - SF_val_miner_total) 
        
        ! reaction intensity [kJ/m2/min]
        i_r = ReactionIntensity(currentPatch%fuel%non_trunk_loading/0.45_r8,             &
          currentPatch%fuel%SAV_notrunks, beta_ratio,                                    &
          currentPatch%fuel%average_moisture_notrunks, currentPatch%fuel%MEF_notrunks)
   
        ! heat of preignition [kJ/kg] 
        q_ig = HeatofPreignition(currentPatch%fuel%average_moisture_notrunks)

        ! effective heating number [unitless]
        eps = EffectiveHeatingNumber(currentPatch%fuel%SAV_notrunks)
        
        ! wind factor [unitless]
        phi_wind = WindFactor(currentSite%fireWeather%effective_windspeed, beta_ratio,      &
         currentPatch%fuel%SAV_notrunks)

        ! propagating flux [unitless]       
        xi = PropagatingFlux(beta, currentPatch%fuel%SAV_notrunks)
        
        ! forward rate of spread [m/min]
        currentPatch%ROS_front = ForwardRateOfSpread(currentPatch%fuel%bulk_density_notrunks, &
         eps, q_ig, i_r, xi, phi_wind)

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
    !  Right now also calculates the area burnt...
    !
    use FatesConstantsMod, only : m2_per_km2
    use SFEquationsMod,    only : FireDuration, LengthToBreadth
    use SFEquationsMod,    only : AreaBurnt, FireSize, FireIntensity
    use SFParamsMod,       only : SF_val_fire_threshold, SF_val_rxfire_minthreshold, &
    SF_val_rxfire_maxthreshold, SF_val_rxfire_fuel_min, SF_val_rxfire_fuel_max
    use EDParamsMod,       only : rxfire_switch

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    
    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch                    ! patch object
    real(r8)                        :: fuel_consumed(num_fuel_classes) ! fuel consumed [kgC/m2]
    real(r8)                        :: tree_fraction_patch             ! treed fraction on patch [0-1]
    real(r8)                        :: length_to_breadth               ! length to breadth ratio of fire ellipse (unitless)
    real(r8)                        :: fire_size                       ! size of fire [m2]
    real(r8)                        :: area_burnt                      ! area burnt [m2/km2]
    logical                         :: is_rxfire                       ! is it a prescribed fire?
    logical                         :: rx_man                          ! prescribed fire use human ignition 
    logical                         :: rx_hyb                          ! prescribed fire due to both lightning strike and human ignition
    logical                         :: managed_wildfire                ! is it a wildfire with FI lower than the max rxfire intensity?[can either be Rx fire or wildfire]
    logical                         :: true_wildfire                   ! is it a wildfire that cannot be managed?
    logical                         :: is_wildfire                     ! combine both managed and true wildfire for now
    
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
        currentPatch%FI = 0.0_r8 
        currentPatch%fire = 0
        currentPatch%rxfire = 0
        currentPatch%rxfire_FI = 0.0_r8
        
        if (currentSite%NF > 0.0_r8 .or. currentSite%fireWeather%rx_flag .eq. itrue) then
          
          ! fire intensity [kW/m]
          currentPatch%FI = FireIntensity(currentPatch%TFC_ROS/0.45_r8, currentPatch%ROS_front/60.0_r8)

          ! Decide if prescribed fire or wildfire happen 
          ! prescribed fire and wildfire cannot happen on the same patch

          ! store some contion check here to simplify decision tree
          
          rx_man = (currentPatch%FI > SF_val_rxfire_minthreshold .and. &
          currentPatch%FI < SF_val_rxfire_maxthreshold .and. &
          currentSite%NF == 0.0_r8)

          rx_hyb = (currentPatch%FI < SF_val_fire_threshold .and. &
          currentPatch%FI > SF_val_rxfire_minthreshold .and. &
          currentPatch%FI < SF_val_rxfire_maxthreshold .and. &
          currentSite%NF > 0.0_r8)

          is_rxfire = (rx_man .or. rx_hyb)

          managed_wildfire = (currentSite%NF > 0.0_r8 .and. &
          currentPatch%FI > SF_val_fire_threshold .and. &
          currentPatch%FI < SF_val_rxfire_maxthreshold)

          true_wildfire = (currentSite%NF > 0.0_r8 .and. &
          currentPatch%FI > SF_val_fire_threshold .and. &
          currentPatch%FI > SF_val_rxfire_maxthreshold)

          is_wildfire = (managed_wildfire .or. true_wildfire)

          if (currentSite%fireWeather%rx_flag == itrue .and. & ! burn window check
          currentPatch%fuel%non_trunk_loading > SF_val_rxfire_fuel_min .and. & ! fuel load check 
          currentPatch%fuel%non_trunk_loading < SF_val_rxfire_fuel_max) then
            currentSite%rxfire_area_fuel = currentSite%rxfire_area_fuel + currentPatch%area ! record burnable area after fuel load check
            if (is_rxfire) then
              currentSite%rxfire_area_fi = currentSite%rxfire_area_fi + currentPatch%area ! record burnable area after FI check
              currentPatch%rxfire = 1
            else if (is_wildfire) then
              currentPatch%fire = 1
            end if

          else  ! not a patch suitable for conducting prescribed fire or rxfire is not even turned on
            ! track wildfires greater than kW/m energy threshold
            if (currentPatch%FI > SF_val_fire_threshold) then 
              currentPatch%fire = 1 
            end if

          end if

          if (currentPatch%fire == itrue) then
            currentSite%NF_successful = currentSite%NF_successful + &
            currentSite%NF*currentSite%FDI*currentPatch%area/area
          else if (currentPatch%rxfire == itrue) then
            currentPatch%rxfire_FI = currentPatch%FI
          end if
          
        end if
      end if

      currentPatch => currentPatch%younger
    end do    

  end subroutine CalculateSurfaceFireIntensity
   
  !---------------------------------------------------------------------------------------

  subroutine PassiveActiveCrownFireCheck(currentSite)
    !
    ! DESCRIPTION:
    ! 1) Calculate theoretical active crown fire ROS (R_active) using fuel model 10 and Rothermel's ROS model
    ! XLG: the calculation of R_active can also incluede canopy bulk density (CBD) and canopy water content (see EQ. 10),
    ! but those are currently ignored in Scott & Reinhardt 2001 due to contrasting effects of CBD on R_active in lit.,
    ! and the difficulty of obtaining a base-line FME (see discussion in Scott & Reinhardt 2001 pg12-13) 
    ! 2) Calculate the critical mass flow rate for sustaining a fully active crown fire given current
    ! canopy bulk density, which is R'_active in Scott & Reinhardt 2001; 
    ! 3) Calculate the critical open wind speed for sustaining a fully active crown fire, AKA crowning index (CI)
    ! using FM10 fuel characteristics 
    ! and the theoretical surface fire ROS (R'_SA, will be used in calculating final ROS for crown fire) using CI
    ! XLG: I'm not sure what SAV value they used in Scott & Reinhardt 2001 to calculate CI,
    ! Sam. L used SAV values from the BEHAVE model, but when use that to calculate some constants and comapre
    ! to EQ. 20 it doesn't seem right to me. This is something we have to discuss
    ! 4) Calculate the critical surface fire ROS (R'_initiation) for initiating a crown fire 
    ! 5) Determine passive or active crown fire by comparing current surface ROS_front to R'_initiation, and
    ! R_active to R'_active. Passive: ROS_front > R'_initiation but R_active < R'_active; 
    ! Active: ROS_front > R'_initiation and R_active > R'_active. Note: FI > FI_init is equivalent to ROS_front > R'_initiation
    ! 6) Calculate new ROS and FI and update ROS_front and FI   
    !
    use SFParamsMod,    only : SF_val_miner_total, SF_val_part_dens, SF_val_drying_ratio
    use EDParamsMod,    only : crown_fire_switch
    use EDTypesMod,     only : CalculateTreeGrassAreaSite
    use SFEquationsMod, only : OptimumPackingRatio, ReactionIntensity
    use SFEquationsMod, only : HeatofPreignition, EffectiveHeatingNumber
    use SFEquationsMod, only : WindFactor, PropagatingFlux
    use SFEquationsMod, only : ForwardRateOfSpread
    use SFEquationsMod, only : PassiveCrownFireIntensity, HeatReleasePerArea
    use SFEquationsMod, only : CrowningIndex, CrownFireIntensity

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite

    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch                    ! patch object


    real(r8)                        :: FI_init           ! critical surface fire intensity for initiating a crown fire [kW/m or kJ/m/s]
    real(r8)                        :: HPA               ! heat release per unit area [kW/m2]
    real(r8)                        :: ROS_init          ! critical surface ROS for initiating a crown fire [m/min]
    real(r8)                        :: ROS_active_min    ! critical ROS for sustaining a fully active crown fire [m/min]
    real(r8)                        :: canopy_frac_burnt ! canopy fraction burnt due to crown fire [fraction]
    real(r8)                        :: ROS_final         ! final ROS when a crown fire happens [m/min]
    real(r8)                        :: FI_final          ! final fireline intensity with crown consumption [kW/m or kJ/m/s]


   ! Local variables for calculating ROS_active 
    real(r8)                        :: ROS_active        ! theoretical crown fire rate of spread using FM 10 fuels [m/min]
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
  
                        
    ! local variables for calculating ROS_SA based on current patch fuel condition
    real(r8)                        :: ROS_SA            ! surface ROS calculated at open wind speed of CI [m/min]
    real(r8)                        :: beta         ! packing ratio of current patch [unitless]
    real(r8)                        :: beta_op      ! optimum packing ratio of current patch [unitless]
    real(r8)                        :: beta_ratio   ! relative packing ratio of current patch [unitless]
    real(r8)                        :: i_r          ! reaction intensity of current patch [kJ/m2/min]
    real(r8)                        :: xi           ! propagating flux ratio of current patch [unitless]
    real(r8)                        :: eps          ! effective heating number of current patch [unitless]
    real(r8)                        :: tree_fraction     ! site-level tree fraction [0-1]
    real(r8)                        :: grass_fraction    ! site-level grass fraction [0-1]
    real(r8)                        :: bare_fraction     ! site-level bare ground fraction [0-1]
    real(r8)                        :: CI                ! crowning index: open wind speed to sustain active crown fire [km/hr]
    real(r8)                        :: CI_effective      ! effective wind speed using CI [m/min]
    real(r8)                        :: phi_wind_SA       ! wind factor for calculating ROS_SA [unitless]
    real(r8)                        :: q_ig         ! heat of pre-ignition of current patch [kJ/kg]
   
                           

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
    real(r8),parameter  :: km_per_hr_to_m_per_min = 16.6667_r8         ! convert km/hour to m/min for wind speed
    real(r8), parameter :: wind_atten_tree = 0.4_r8                    ! wind attenuation factor for tree fraction
    real(r8), parameter :: wind_atten_grass = 0.6_r8                   ! wind attenuation factor for grass fraction

    
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
      ! initialize patch level variables
      currentPatch%passive_crown_fire = 0
      currentPatch%active_crown_fire = 0

      if (currentPatch%nocomp_pft_label /= nocomp_bareground .and.         &
      currentPatch%fire == itrue .and. crown_fire_switch) then
        ! calculate passive crown fire intensity, the minimum surface FI to initiate a crown fire
        FI_init = PassiveCrownFireIntensity(currentPatch%fuel%canopy_base_height)
        
        ! check if there is a crown fire 
        if (currentPatch%FI > FI_init ) then
          ! calculate ROS_active 
          fuel_1h     = fuel_1h_ton * tonnes_acre_to_kg_m2
          fuel_10h    = fuel_10h_ton * tonnes_acre_to_kg_m2
          fuel_100h   = fuel_100h_ton * tonnes_acre_to_kg_m2
          fuel_live   = fuel_live_ton * tonnes_acre_to_kg_m2
          total_fuel  = (fuel_1h + fuel_10h + fuel_100h + fuel_live) ! kg biomass /m2
          fuel_sav1h  = sav_1h_ft * sqft_cubicft_to_sqm_cubicm
          fuel_sav10h = sav_10h_ft * sqft_cubicft_to_sqm_cubicm
          fuel_sav100h = sav_100h_ft * sqft_cubicft_to_sqm_cubicm
          fuel_savlive  = sav_live_ft * sqft_cubicft_to_sqm_cubicm
          fuel_moist1h     = exp(-1.0_r8 * ((fuel_sav1h/SF_val_drying_ratio) * currentSite%fireWeather%fire_weather_index))
          fuel_moist10h    = exp(-1.0_r8 * ((fuel_sav10h/SF_val_drying_ratio) * currentSite%fireWeather%fire_weather_index))
          fuel_moist100h   = exp(-1.0_r8 * ((fuel_sav100h/SF_val_drying_ratio) * currentSite%fireWeather%fire_weather_index))
          fuel_moist_live  = exp(-1.0_r8 * ((fuel_savlive/SF_val_drying_ratio) * currentSite%fireWeather%fire_weather_index))
          fuel_depth       = fuel_depth_ft *0.3048_r8           !convert to meters
          fuel_bd          = total_fuel/fuel_depth              !fuel bulk density (kg biomass/m3)

          fuel_sav         = fuel_sav1h *(fuel_1h/total_fuel) + fuel_sav10h*(fuel_10h/total_fuel) + & 
                             fuel_sav100h*(fuel_100h/total_fuel) + fuel_savlive*(fuel_live/total_fuel)

          fuel_eff_moist   = fuel_moist1h *(fuel_1h/total_fuel) + fuel_moist10h*(fuel_10h/total_fuel) + & 
                             fuel_moist100h*(fuel_100h/total_fuel) + fuel_moist_live*(fuel_live/total_fuel)

          net_fuel  = total_fuel * (1.0_r8 - SF_val_miner_total)

          beta_fm10 = fuel_bd / SF_val_part_dens
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

          midflame_wind = currentSite%wind * 0.40_r8 ! Scott & Reinhardt 2001 use 40% of open wind speed as effective wind speed
                                                     ! XLG: should we scwitch to the way FATES calculateS effective wind speed?
          phi_wind_fm10 = WindFactor(midflame_wind, beta_ratio_fm10, fuel_sav)

          xi_fm10 = PropagatingFlux(beta_fm10, fuel_sav)
          
          ROS_active = ForwardRateOfSpread(fuel_bd, eps_fm10, q_ig_fm10, i_r_fm10, &
          xi_fm10, phi_wind_fm10)

          ! Calculate ROS_acitive_min, EQ 14 in Scott & Reinhardt 2001

          ROS_active_min = 3.0_r8 / currentPatch%fuel%canopy_bulk_density

          ! Calculate ROS_SA using current patch fuel conditions
          beta = currentPatch%fuel%bulk_density_notrunks/SF_val_part_dens
          beta_op = OptimumPackingRatio(currentPatch%fuel%SAV_notrunks)

          if (beta_op < nearzero) then 
            beta_ratio = 0.0_r8
          else
            beta_ratio = beta/beta_op 
          end if

          i_r = ReactionIntensity(currentPatch%fuel%non_trunk_loading/0.45_r8, &
            currentPatch%fuel%SAV_notrunks, beta_ratio,                        &
            currentPatch%fuel%average_moisture_notrunks, currentPatch%fuel%MEF_notrunks)
          q_ig = HeatofPreignition(currentPatch%fuel%average_moisture_notrunks)
          eps = EffectiveHeatingNumber(currentPatch%fuel%SAV_notrunks)

          ! Calculate crowning index, which is used for calculating phi_wind
          ! CI is also calculated using fuel characteristics of FM 10 
          CI = CrowningIndex(eps_fm10, q_ig_fm10, i_r_fm10, &
          currentPatch%fuel%canopy_bulk_density )
          CI = CI * km_per_hr_to_m_per_min  ! convert to m/min

          ! calculate effective wind speed at CI
          ! XLG: this is the disconnection from ROS_active calculation, where they used midflame wind speed 
          ! as effective wind speed. But it might also make sense as ROS_active is theoretical crown fire
          ! rate of spread, and ROS_SA is just the surface fire rate of spread at crowning index wind speed,
          ! so it makes senses to consider wind attenuation for grass cover as it's still surface fire??

          call CalculateTreeGrassAreaSite(currentSite,tree_fraction, grass_fraction, bare_fraction)
          CI_effective = CI * (tree_fraction*wind_atten_tree + &
          (grass_fraction + bare_fraction)*wind_atten_grass)
          ! phi_wind
          phi_wind_SA = WindFactor(CI_effective, beta_ratio,   &
              currentPatch%fuel%SAV_notrunks)
          xi = PropagatingFlux(beta, currentPatch%fuel%SAV_notrunks)
          ROS_SA = ForwardRateOfSpread(currentPatch%fuel%bulk_density_notrunks, &
              eps, q_ig, i_r, xi, phi_wind_SA)
          
          ! Calculate ROS_init, EQ. 12 in Scott & Reinhardt 2001
          ! first calculate heat release per unit area [kW/m2]
          HPA = HeatReleasePerArea(currentPatch%fuel%SAV_notrunks, i_r)
          ROS_init = (60.0_r8 * FI_init) / HPA 

          ! Now check if there is passive or active crown fire and calculate crown fraction burnt (CFB)
          ! XLG: there are alternative ways to calculate CFB, see pg 39-41 in Scott & Reinhardt 2001

          if (ROS_active >= ROS_active_min) then ! FI >= FI_init and ROS_active >= ROS_active_min
            currentPatch%active_crown_fire = 1 
            ! for active crown fire we set CFB to 1
            canopy_frac_burnt = 1.0_r8
          else if (ROS_active < ROS_active_min .and. &  ! FI >= FI_init but ROS_active < ROS_active_min
            currentPatch%ROS_front > ROS_init .and. &   ! XLG: it seems redudant when FI > FI_init is true, but calculation of ROS_init is 
            currentPatch%ROS_front < ROS_SA) then       ! different from ROS_front, let's check to be safe. I'm uncomfortable when calculations
              currentPatch%passive_crown_fire = 1       ! of all kinds of ROS are so disconnected from each other

              ! calculate crown fraction burnt EQ. 28 in Scott & Reinhardt 2001
              canopy_frac_burnt = min(1.0_r8, (currentPatch%ROS_front - ROS_init) / &
              (ROS_SA - ROS_init))
          else ! a condition where crown fire cessation happens??
            canopy_frac_burnt = 0.0_r8 
          end if
          ! calculate ROS_final EQ. 21 in Scott & Reinhardt 2001
          ROS_final = currentPatch%ROS_front + canopy_frac_burnt * &
                      ( ROS_active - currentPatch%ROS_front)
          ! update ROS_front with ROS_final
          currentPatch%ROS_front = ROS_final
          ! update fire intensity by accounting for burned canopy fuels
          ! EQ. 22 in Scott & Reinhardt 2001
          FI_final = CrownFireIntensity(HPA, currentPatch%fuel%canopy_fuel_load, &
          canopy_frac_burnt, ROS_final)

          if(write_SF == itrue)then
            if ( hlm_masterproc == itrue ) write(fates_log(),*) 'FI_final',FI_final
         endif

          ! only update FI when CFB > 0
          if (canopy_frac_burnt > 0.0_r8) then
            currentPatch%FI = FI_final
          end if

        end if ! end check if there is a crown fire 
      end if ! if there is a fire 

          
      currentPatch => currentPatch%younger 
    end do ! end patch loop

  end subroutine PassiveActiveCrownFireCheck

  !---------------------------------------------------------------------------------------
  
  subroutine CalculateAreaBurnt(currentSite)
    !
    !  DESCRIPTION:
    !  Calculates area burnt for each patch of a site
    !
    use FatesConstantsMod, only : m2_per_km2
    use SFEquationsMod,    only : FireDuration, LengthToBreadth
    use SFEquationsMod,    only : AreaBurnt, FireSize
    use SFParamsMod,       only : SF_val_fire_threshold

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    
    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch                    ! patch object
    real(r8)                        :: tree_fraction_patch             ! treed fraction on patch [0-1]
    real(r8)                        :: length_to_breadth               ! length to breadth ratio of fire ellipse (unitless)
    real(r8)                        :: fire_size                       ! size of fire [m2]
    real(r8)                        :: area_burnt                      ! area burnt [m2/km2]
    
    currentPatch => currentSite%oldest_patch 
    do while (associated(currentPatch))

      if (currentPatch%nocomp_pft_label /= nocomp_bareground) then

        ! initialize patch parameters to zero
        currentPatch%FD = 0.0_r8
        currentPatch%frac_burnt = 0.0_r8

        if (currentSite%NF > 0.0_r8 .and. currentPatch%FI > SF_val_fire_threshold) then

          ! fire duration [min]
          currentPatch%FD = FireDuration(currentSite%FDI)
          
          ! length-to-breadth ratio of fire ellipse [unitless]
          tree_fraction_patch  = currentPatch%total_tree_area/currentPatch%area
          length_to_breadth = LengthToBreadth(currentSite%fireWeather%effective_windspeed, tree_fraction_patch)

          ! fire size [m2]
          fire_size = FireSize(length_to_breadth, currentPatch%ROS_back,              &
              currentPatch%ROS_front, currentPatch%FD)

          ! area burnt [m2/km2]
          area_burnt = AreaBurnt(fire_size, currentSite%NF, currentSite%FDI)
          
          ! convert to area burned per area patch per day
          ! i.e., fraction of the patch burned on that day
          currentPatch%frac_burnt = min(0.99_r8, area_burnt/m2_per_km2)
          
        end if
      end if
      currentPatch => currentPatch%younger
    end do    

  end subroutine CalculateAreaBurnt
   
  !---------------------------------------------------------------------------------------

  !*****************************************************************
  subroutine CalculateRxfireAreaBurnt ( currentSite )
  !*****************************************************************

    !returns burned fraction for prescribed fire per patch by first checking
    !if total burnable fraction at site level is greater than user defined fraction of site area 
    !if yes, calculate burned fraction as (user defined frac / total burnable frac)

    use SFParamsMod,       only : SF_val_rxfire_AB !user defined prescribed fire area in fraction per day to reflect burning capacity

    ! ARGUMENTS
    type(ed_site_type), intent(inout), target :: currentSite

    !LOCALS
    type(fates_patch_type), pointer :: currentPatch  

    real(r8) :: total_burnable_frac        ! total fractional land area that can apply prescribed fire after condition checks at site level

    ! Testing parameters 
    real(r8), parameter :: min_frac_site = 0.1_r8 

    ! initialize site variables
    currentSite%rxfire_area_final = 0.0_r8 
    total_burnable_frac = 0.0_r8

    ! update total burnable fraction
    total_burnable_frac = currentSite%rxfire_area_fi / AREA
   
    currentPatch => currentSite%oldest_patch;

    do while(associated(currentPatch))

      if(currentPatch%nocomp_pft_label .ne. nocomp_bareground)then
        currentPatch%rxfire_frac_burnt = 0.0_r8
        if (currentPatch%rxfire .eq. itrue .and. & 
        total_burnable_frac .ge. min_frac_site ) then
          currentSite%rxfire_area_final = currentSite%rxfire_area_final + currentPatch%area ! the final burned total land area 
          currentPatch%rxfire_frac_burnt = min(0.99_r8, (SF_val_rxfire_AB / total_burnable_frac))
        else
          currentPatch%rxfire = 0 ! update rxfire occurence at patch 
          currentPatch%rxfire_FI = 0.0_r8
        end if
      end if

      currentPatch => currentPatch%younger;  
    end do ! end patch loop

  end subroutine CalculateRxfireAreaBurnt

  
!---------------------------------------------------------------------------------------

  
  !*****************************************************************
  subroutine  crown_scorching ( currentSite ) 
  !*****************************************************************

    !currentPatch%FI       average fire intensity of flaming front during day.  kW/m.
    !currentPatch%SH(pft)  scorch height for all cohorts of a given PFT on a given patch (m)

    type(ed_site_type), intent(in), target :: currentSite

    type(fates_patch_type), pointer  :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    real(r8) ::  tree_ag_biomass ! total amount of above-ground tree biomass in patch. kgC/m2
    real(r8) ::  leaf_c          ! leaf carbon      [kg]
    real(r8) ::  sapw_c          ! sapwood carbon   [kg]
    real(r8) ::  struct_c        ! structure carbon [kg]

    integer  ::  i_pft


    currentPatch => currentSite%oldest_patch;  
    do while(associated(currentPatch)) 

       if(currentPatch%nocomp_pft_label .ne. nocomp_bareground)then
       
       tree_ag_biomass = 0.0_r8
       if (currentPatch%fire == 1 .or. currentPatch%rxfire == 1) then
          currentCohort => currentPatch%tallest;
          do while(associated(currentCohort))  
             if ( prt_params%woody(currentCohort%pft) == itrue) then !trees only

                leaf_c = currentCohort%prt%GetState(leaf_organ, carbon12_element)
                sapw_c = currentCohort%prt%GetState(sapw_organ, carbon12_element)
                struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
                
                tree_ag_biomass = tree_ag_biomass + &
                      currentCohort%n * (leaf_c + & 
                      prt_params%allom_agb_frac(currentCohort%pft)*(sapw_c + struct_c))
             endif !trees only
             currentCohort=>currentCohort%shorter;
          enddo !end cohort loop

          do i_pft=1,numpft
             if (tree_ag_biomass > 0.0_r8  .and. prt_params%woody(i_pft) == itrue) then 
                
                !Equation 16 in Thonicke et al. 2010 !Van Wagner 1973 EQ8 !2/3 Byram (1959)
                currentPatch%Scorch_ht(i_pft) = EDPftvarcon_inst%fire_alpha_SH(i_pft) * (currentPatch%FI**0.667_r8)

                if(write_SF == itrue)then
                   if ( hlm_masterproc == itrue ) write(fates_log(),*) 'currentPatch%SH',currentPatch%Scorch_ht(i_pft)
                endif
             else
                currentPatch%Scorch_ht(i_pft) = 0.0_r8
             endif ! tree biomass
          end do

       endif !fire
       endif !nocomp_pft_label

       currentPatch => currentPatch%younger;  
    enddo !end patch loop

  end subroutine crown_scorching

  !*****************************************************************
  subroutine  crown_damage ( currentSite )
    !*****************************************************************

    !returns the updated currentCohort%fraction_crown_burned for each tree cohort within each patch.
    !currentCohort%fraction_crown_burned is the proportion of crown affected by fire

    type(ed_site_type), intent(in), target :: currentSite

    type(fates_patch_type) , pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort
    real(r8)                      :: crown_depth    ! Depth of crown in meters
    
    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch)) 

       if(currentPatch%nocomp_pft_label .ne. nocomp_bareground)then
       if (currentPatch%fire == 1 .or. currentPatch%rxfire == 1) then

          currentCohort=>currentPatch%tallest

          do while(associated(currentCohort))  
             currentCohort%fraction_crown_burned = 0.0_r8
             if ( prt_params%woody(currentCohort%pft) == itrue) then !trees only
                ! Flames lower than bottom of canopy. 
                ! c%height is height of cohort
                
                ! XLG: can we simply change this logic here to just use whethere 
                ! there is a passive or active crown fire? if no, crown fraction
                ! burnt is zero, if yes, we can calculate CFB given active or passive scenario?

                call CrownDepth(currentCohort%height,currentCohort%pft,crown_depth)
                
                if (currentPatch%Scorch_ht(currentCohort%pft) < &
                     (currentCohort%height-crown_depth)) then 
                   currentCohort%fraction_crown_burned = 0.0_r8
                else
                   ! Flames part of way up canopy. 
                   ! Equation 17 in Thonicke et al. 2010. 
                   ! flames over bottom of canopy, CFB depends on whether it's active or passive crown fire
                   if ((currentCohort%height > 0.0_r8).and.(currentPatch%Scorch_ht(currentCohort%pft) >=  &
                        (currentCohort%height-crown_depth))) then 
                          if (currentPatch%active_crown_fire == 1) then
                            currentCohort%fraction_crown_burned = 1.0_r8
                          else if (currentPatch%passive_crown_fire == 1) then 
                            ! XLG: should we use CFB calculated in the PassiveActiveCrownFireCheck routine??
                            ! since FI is calculated using that CFB, there will be mismatch between actual
                            ! biomass consumed using the CFB here (thus resulted fire intensity) and the 
                            ! calculated fire intenaity using the other CFB
                            currentCohort%fraction_crown_burned = (currentPatch%Scorch_ht(currentCohort%pft) - &
                             (currentCohort%height - crown_depth))/crown_depth
                          end if
                   endif

                endif
                ! Check for strange values. 
                currentCohort%fraction_crown_burned = min(1.0_r8, max(0.0_r8,currentCohort%fraction_crown_burned))              
             endif !trees only
             !shrink canopy to account for burnt section.     
             !currentCohort%canopy_trim = min(currentCohort%canopy_trim,(1.0_r8-currentCohort%fraction_crown_burned)) 

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop
       endif !fire?
       endif !nocomp_pft_label check

       currentPatch => currentPatch%younger;

    enddo !end patch loop

  end subroutine crown_damage

  !*****************************************************************
  subroutine  cambial_damage_kill ( currentSite ) 
    !*****************************************************************
    ! routine description.
    ! returns the probability that trees dies due to cambial char
    ! currentPatch%tau_l = duration of lethal stem heating (min). Calculated at patch level.

    type(ed_site_type), intent(in), target :: currentSite

    type(fates_patch_type) , pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    real(r8) :: tau_c !critical time taken to kill cambium (minutes) 
    real(r8) :: bt    !bark thickness in cm.

    currentPatch => currentSite%oldest_patch;  

    do while(associated(currentPatch)) 

       if(currentPatch%nocomp_pft_label .ne. nocomp_bareground)then

       if (currentPatch%fire == 1 .or. currentPatch%rxfire == 1) then
          currentCohort => currentPatch%tallest;
          do while(associated(currentCohort))  
             if ( prt_params%woody(currentCohort%pft) == itrue) then !trees only
                ! Equation 21 in Thonicke et al 2010
                bt = EDPftvarcon_inst%bark_scaler(currentCohort%pft)*currentCohort%dbh ! bark thickness. 
                ! Equation 20 in Thonicke et al. 2010. 
                tau_c = 2.9_r8*bt**2.0_r8 !calculate time it takes to kill cambium (min)
                ! Equation 19 in Thonicke et al. 2010
                if ((currentPatch%tau_l/tau_c) >= 2.0_r8) then
                   currentCohort%cambial_mort = 1.0_r8
                else
                   if ((currentPatch%tau_l/tau_c) > 0.22_r8) then
                      currentCohort%cambial_mort = (0.563_r8*(currentPatch%tau_l/tau_c)) - 0.125_r8
                   else
                      currentCohort%cambial_mort = 0.0_r8
                   endif
                endif
             endif !trees 

             currentCohort => currentCohort%shorter;

          enddo !end cohort loop
       endif !fire?
       endif !nocomp_pft_label check

       currentPatch=>currentPatch%younger;

    enddo !end patch loop

  end subroutine cambial_damage_kill

  !*****************************************************************
  subroutine  post_fire_mortality ( currentSite )
  !*****************************************************************

    !  returns the updated currentCohort%fire_mort value for each tree cohort within each patch.
    !  currentCohort%fraction_crown_burned is proportion of crown affected by fire
    !  currentCohort%crownfire_mort  probability of tree post-fire mortality due to crown scorch
    !  currentCohort%cambial_mort  probability of tree post-fire mortality due to cambial char
    !  currentCohort%fire_mort  post-fire mortality from cambial and crown damage assuming two are independent.

    type(ed_site_type), intent(in), target :: currentSite

    type(fates_patch_type),  pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    currentPatch => currentSite%oldest_patch

    do while(associated(currentPatch)) 

       if(currentPatch%nocomp_pft_label .ne. nocomp_bareground)then

       if (currentPatch%fire == 1 .or. currentPatch%rxfire == 1) then 
          currentCohort => currentPatch%tallest
          do while(associated(currentCohort))  
             currentCohort%fire_mort = 0.0_r8
             currentCohort%crownfire_mort = 0.0_r8
             currentCohort%rxfire_mort = 0.0_r8
             currentCohort%rxcrownfire_mort = 0.0_r8
             currentCohort%rxcambial_mort = 0.0_r8
             if ( prt_params%woody(currentCohort%pft) == itrue) then
                ! Equation 22 in Thonicke et al. 2010. 
                currentCohort%crownfire_mort = EDPftvarcon_inst%crown_kill(currentCohort%pft)*currentCohort%fraction_crown_burned**3.0_r8
                ! Equation 18 in Thonicke et al. 2010. 
                currentCohort%fire_mort = max(0._r8,min(1.0_r8,currentCohort%crownfire_mort+currentCohort%cambial_mort- &
                     (currentCohort%crownfire_mort*currentCohort%cambial_mort)))  !joint prob.   
             else
                currentCohort%fire_mort = 0.0_r8 !Set to zero. Grass mode of death is removal of leaves.
             endif !trees

             ! now decide which type of post-fire mortality, prescribed fire or wildfire?
             if (currentPatch%rxfire == itrue .and. currentPatch%fire == ifalse) then
              currentCohort%rxfire_mort = currentCohort%fire_mort
              currentCohort%rxcrownfire_mort = currentCohort%crownfire_mort
              currentCohort%rxcambial_mort = currentCohort%cambial_mort
              currentCohort%fire_mort = 0.0_r8
              currentCohort%crownfire_mort = 0.0_r8
              currentCohort%cambial_mort = 0.0_r8
             end if
             
             currentCohort => currentCohort%shorter

          enddo !end cohort loop
       endif !fire?
       endif !nocomp_pft_label check

       currentPatch => currentPatch%younger

    enddo !end patch loop

  end subroutine post_fire_mortality

  ! ============================================================================
end module SFMainMod
