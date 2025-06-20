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
      call CalculatePostFireMortality(currentSite)
    end if

  end subroutine DailyFireModel

  !---------------------------------------------------------------------------------------
 
  subroutine UpdateFireWeather(currentSite, bc_in)
    !
    !  DESCRIPTION:
    !  Updates the site's fire weather index and calculates effective windspeed based on 
    !   vegetation characteristics
    !
    !  Currently we use tree and grass fraction averaged over whole grid (site) to 
    !  prevent extreme divergence

    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : sec_per_day, sec_per_min
    use EDTypesMod,        only : CalculateTreeGrassAreaSite

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
    !
    use SFEquationsMod, only : FireIntensity
    use SFParamsMod,    only : SF_val_fire_threshold

    ! ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    
    ! LOCALS:
    type(fates_patch_type), pointer :: currentPatch                    ! patch object
    real(r8)                        :: fuel_consumed(num_fuel_classes) ! fuel consumed [kgC/m2]
    
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
        
        if (currentSite%NF > 0.0_r8) then
          
          ! fire intensity [kW/m]
          currentPatch%FI = FireIntensity(currentPatch%TFC_ROS/0.45_r8, currentPatch%ROS_front/60.0_r8)

          ! track fires greater than kW/m energy threshold
          if (currentPatch%FI > SF_val_fire_threshold) then 
            currentPatch%fire = 1 
            currentSite%NF_successful = currentSite%NF_successful + &
              currentSite%NF * currentSite%FDI*currentPatch%area / area
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
    use SFParamsMod,       only : SF_val_fire_threshold

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
        currentPatch%frac_burnt = 0.0_r8

        if (currentSite%NF > 0.0_r8 .and. currentPatch%FI > SF_val_fire_threshold) then

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
          currentPatch%frac_burnt = min(max_frac_burnt, area_burnt/m2_per_km2)
          
        end if
      end if
      currentPatch => currentPatch%younger
    end do    

  end subroutine CalculateAreaBurnt
   
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
