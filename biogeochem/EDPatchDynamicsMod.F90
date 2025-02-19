module EDPatchDynamicsMod
  ! ============================================================================
  ! Controls formation, creation, fusing and termination of patch level processes. 
  ! ============================================================================
  use FatesGlobals         , only : fates_log
  use FatesGlobals         , only : FatesWarn,N2S,A2S
  use FatesInterfaceTypesMod, only : hlm_freq_day
  use FatesInterfaceTypesMod, only : hlm_current_tod
  use EDPftvarcon          , only : EDPftvarcon_inst
  use EDPftvarcon          , only : GetDecompyFrac
  use PRTParametersMod      , only : prt_params
  use EDCohortDynamicsMod  , only : fuse_cohorts
  use EDTypesMod           , only : area_site => area
  use ChecksBalancesMod    , only : PatchMassStock
  use FatesLitterMod       , only : ncwd
  use FatesLitterMod       , only : ndcmpy
  use FatesLitterMod       , only : litter_type
  use FatesConstantsMod    , only : n_dbh_bins 
  use FatesLitterMod       , only : adjust_SF_CWD_frac
  use EDTypesMod           , only : homogenize_seed_pfts
  use EDTypesMod           , only : area
  use FatesConstantsMod    , only : patchfusion_dbhbin_loweredges
  use EDtypesMod           , only : force_patchfuse_min_biomass
  use EDTypesMod           , only : ed_site_type
  use FatesPatchMod        , only : fates_patch_type
  use EDTypesMod           , only : set_patchno 
  use FatesCohortMod       , only : fates_cohort_type
  use EDTypesMod           , only : site_massbal_type
  use EDTypesMod           , only : site_fluxdiags_type
  use EDTypesMod           , only : elem_diag_type
  use EDTypesMod           , only : min_patch_area
  use EDTypesMod           , only : min_patch_area_forced
  use EDParamsMod          , only : regeneration_model
  use FatesInterfaceTypesMod, only : numpft
  use FatesConstantsMod     , only : dtype_ifall
  use FatesConstantsMod     , only : dtype_ilog
  use FatesConstantsMod     , only : dtype_ifire
  use FatesConstantsMod     , only : dtype_ilandusechange
  use FatesConstantsMod    , only : ican_upper
  use PRTGenericMod        , only : num_elements
  use PRTGenericMod        , only : element_list
  use FatesFuelClassesMod  , only : fuel_classes
  use FatesConstantsMod    , only : N_DIST_TYPES
  use EDTypesMod           , only : AREA_INV
  use EDTypesMod           , only : dump_site
  use FatesConstantsMod    , only : rsnbl_math_prec
  use FatesConstantsMod    , only : fates_tiny
  use FatesConstantsMod    , only : nocomp_bareground
  use FatesInterfaceTypesMod    , only : hlm_use_planthydro
  use FatesInterfaceTypesMod    , only : bc_in_type
  use FatesInterfaceTypesMod    , only : numpft
  use FatesInterfaceTypesMod    , only : hlm_stepsize
  use FatesInterfaceTypesMod    , only : hlm_use_sp
  use FatesInterfaceTypesMod    , only : hlm_use_nocomp
  use FatesInterfaceTypesMod    , only : hlm_use_fixed_biogeog
  use FatesInterfaceTypesMod    , only : hlm_num_lu_harvest_cats
  use FatesInterfaceTypesMod    , only : hlm_use_luh
  use FatesInterfaceTypesMod    , only : hlm_num_luh2_states
  use FatesInterfaceTypesMod    , only : hlm_num_luh2_transitions
  use FatesGlobals         , only : endrun => fates_endrun
  use FatesConstantsMod    , only : r8 => fates_r8
  use FatesConstantsMod    , only : itrue, ifalse
  use FatesConstantsMod    , only : t_water_freeze_k_1atm
  use FatesConstantsMod    , only : TRS_regeneration
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use FatesPlantHydraulicsMod, only : AccumulateMortalityWaterStorage
  use FatesPlantHydraulicsMod, only : DeallocateHydrCohort
  use EDLoggingMortalityMod, only : logging_litter_fluxes 
  use EDLoggingMortalityMod, only : logging_time
  use EDLoggingMortalityMod, only : get_harvest_rate_area
  use EDLoggingMortalityMod, only : get_harvest_rate_carbon
  use EDLoggingMortalityMod, only : get_harvestable_carbon
  use EDLoggingMortalityMod, only : get_harvest_debt
  use FatesLandUseChangeMod, only : GetInitLanduseHarvestRate
  use EDParamsMod          , only : fates_mortality_disturbance_fraction
  use FatesAllometryMod    , only : carea_allom
  use FatesAllometryMod    , only : set_root_fraction
  use FatesConstantsMod    , only : g_per_kg
  use FatesConstantsMod    , only : ha_per_m2
  use FatesConstantsMod    , only : days_per_sec
  use FatesConstantsMod    , only : years_per_day
  use FatesConstantsMod    , only : nearzero
  use FatesConstantsMod    , only : primaryland, secondaryland, pastureland, rangeland, cropland
  use FatesConstantsMod    , only : nocomp_bareground_land
  use FatesConstantsMod    , only : n_landuse_cats
  use FatesLandUseChangeMod, only : GetLanduseTransitionRates
  use FatesLandUseChangeMod, only : GetInitLanduseTransitionRates
  use FatesLandUseChangeMod, only : GetLUHStatedata
  use FatesConstantsMod    , only : fates_unset_r8
  use FatesConstantsMod    , only : fates_unset_int
  use FatesConstantsMod    , only : hlm_harvest_carbon
  use EDCohortDynamicsMod  , only : InitPRTObject
  use ChecksBalancesMod,      only : SiteMassStock
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTLossFluxesMod,       only : PRTBurnLosses
  use FatesInterfaceTypesMod,      only : hlm_parteh_mode
  use PRTGenericMod,          only : prt_carbon_allom_hyp   
  use PRTGenericMod,          only : prt_cnp_flex_allom_hyp
  use SFParamsMod,            only : SF_VAL_CWD_FRAC
  use EDParamsMod,            only : logging_event_code
  use EDParamsMod,            only : logging_export_frac
  use EDParamsMod,            only : maxpatches_by_landuse
  use FatesRunningMeanMod,    only : ema_sdlng_mdd
  use FatesRunningMeanMod,    only : ema_sdlng_emerg_h2o, ema_sdlng_mort_par, ema_sdlng2sap_par
  use FatesRunningMeanMod,    only : ema_24hr, fixed_24hr, ema_lpa, ema_longterm
  use FatesRadiationMemMod,   only : num_swb

  ! CIME globals
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  
  !
  implicit none
  private
  
  public :: spawn_patches
  public :: fuse_patches
  public :: terminate_patches
  public :: patch_pft_size_profile
  public :: disturbance_rates
  public :: check_patch_area
  private:: fuse_2_patches

  character(len=*), parameter, private :: sourcefile = &
        __FILE__

  logical, parameter :: debug = .false.

  ! When creating new patches from other patches, we need to send some of the
  ! litter from the old patch to the new patch.  Likewise, when plants die
  ! from the disturbance, we need to send some amount from the old patch to
  ! the new patch.  If the plant matter falls straight down, then that
  ! would make a case for all of the litter going to the new patch. 
  ! This would be localization=1
  ! But if we think some of the plant matter drifts, or a tall tree lands
  ! outside of its gap, or are afraid of newly formed patches having 
  ! too much burnable material, then we drop the localization from 1 down
  ! a notch.
  ! Note that in all cases, a localization of 0, suggests that litter
  ! is dispensed randomly in space among the area of the new and old
  ! patch combined. A localization of 1 suggests that
  ! all litter is sent to the new patch.

  real(r8), parameter :: existing_litt_localization = 1.0_r8
  real(r8), parameter :: treefall_localization = 0.0_r8
  real(r8), parameter :: burn_localization = 0.0_r8
  real(r8), parameter :: landusechange_localization = 1.0_r8

  integer :: istat           ! return status code
  character(len=255) :: smsg ! Message string for deallocation errors
  character(len=512) :: msg  ! Message string for warnings and logging
  
  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

  ! ============================================================================
  subroutine disturbance_rates( site_in, bc_in)
    !
    ! !DESCRIPTION:
    ! Calculates the fire and mortality related disturbance rates for each patch,
    ! and then determines which is the larger at the patch scale (for now, there an only
    ! be one disturbance type for each timestep.  
    ! all disturbance rates here are per daily timestep. 
    
    ! 2016-2017
    ! Modify to add logging disturbance

    ! !USES:
    use EDMortalityFunctionsMod , only : mortality_rates
    use EDMortalityFunctionsMod , only : ExemptTreefallDist
    ! loging flux
    use EDLoggingMortalityMod , only : LoggingMortality_frac

  
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout) :: site_in
    type(bc_in_type) , intent(in) :: bc_in
    !
    ! !LOCAL VARIABLES:
    type (fates_patch_type) , pointer :: currentPatch
    type (fates_cohort_type), pointer :: currentCohort

    real(r8) :: cmort
    real(r8) :: bmort
    real(r8) :: hmort
    real(r8) :: frmort
    real(r8) :: smort
    real(r8) :: asmort
    real(r8) :: dgmort
    
    real(r8) :: lmort_direct
    real(r8) :: lmort_collateral
    real(r8) :: lmort_infra
    real(r8) :: l_degrad         ! fraction of trees that are not killed but suffer from forest 
                                 ! degradation (i.e. they are moved to newly-anthro-disturbed 
                                 ! secondary forest patch)
    real(r8) :: dist_rate_ldist_notharvested
    integer  :: threshold_sizeclass
    integer  :: i_dist
    integer  :: h_index
    real(r8) :: harvest_rate
    real(r8) :: tempsum
    real(r8) :: mean_temp
    real(r8) :: harvestable_forest_c(hlm_num_lu_harvest_cats)
    integer  :: harvest_tag(hlm_num_lu_harvest_cats)
    real(r8) :: current_fates_landuse_state_vector(n_landuse_cats)  ! [m2/m2]
    real(r8) :: state_vector(n_landuse_cats)
    real(r8), parameter :: max_daily_disturbance_rate = 0.999_r8
    logical  :: site_secondaryland_first_exceeding_min
    real(r8) :: secondary_young_fraction  ! what fraction of secondary land is young secondary land
    !----------------------------------------------------------------------------------------------
    ! Calculate Mortality Rates (these were previously calculated during growth derivatives)
    ! And the same rates in understory plants have already been applied to %dndt
    !----------------------------------------------------------------------------------------------
    
    ! first calculate the fraction of the site that is primary land
    current_fates_landuse_state_vector = site_in%get_current_landuse_statevector()

    ! and get the fraction of secondary land that is young secondary land
    secondary_young_fraction = site_in%get_secondary_young_fraction()

    ! check status of transition_landuse_from_off_to_on flag, and do some error checking on it
    if(site_in%transition_landuse_from_off_to_on) then
       if (sum(current_fates_landuse_state_vector(secondaryland:cropland)) .gt. nearzero) then
          write(fates_log(),*) 'flag for transition_landuse_from_off_to_on is set to true but site is not entirely primaryland'
          write(fates_log(), *) current_fates_landuse_state_vector
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
    endif

    ! get available biomass for harvest for all patches
    call get_harvestable_carbon(site_in, bc_in%site_area, bc_in%hlm_harvest_catnames, harvestable_forest_c)
 
    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))   

       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))        
          ! Mortality for trees in the understorey.
          !currentCohort%patchptr => currentPatch
          mean_temp = currentPatch%tveg24%GetMean()
          call mortality_rates(currentCohort,bc_in,currentPatch%btran_ft,      &
            mean_temp, cmort,hmort,bmort,frmort,smort,asmort,dgmort)
          currentCohort%dmort  = cmort+hmort+bmort+frmort+smort+asmort+dgmort
          call carea_allom(currentCohort%dbh,currentCohort%n,site_in%spread,currentCohort%pft, &
               currentCohort%crowndamage,currentCohort%c_area)

          ! Initialize diagnostic mortality rates
          currentCohort%cmort = cmort
          currentCohort%bmort = bmort
          currentCohort%hmort = hmort
          currentCohort%frmort = frmort
          currentCohort%smort = smort
          currentCohort%asmort = asmort
          currentCohort%dgmort = dgmort
          
          call LoggingMortality_frac(site_in, bc_in, currentCohort%pft, &
                currentCohort%dbh, currentCohort%canopy_layer, &
                lmort_direct,lmort_collateral,lmort_infra,l_degrad,&
                bc_in%hlm_harvest_rates, &
                bc_in%hlm_harvest_catnames, &
                bc_in%hlm_harvest_units, &
                currentPatch%land_use_label, &
                currentPatch%age_since_anthro_disturbance, &
                current_fates_landuse_state_vector, &
                harvestable_forest_c, &
                harvest_tag)
         
          currentCohort%lmort_direct     = lmort_direct
          currentCohort%lmort_collateral = lmort_collateral
          currentCohort%lmort_infra      = lmort_infra
          currentCohort%l_degrad         = l_degrad

          currentCohort => currentCohort%taller
       end do

       currentPatch => currentPatch%younger
    end do

    call get_harvest_debt(site_in, bc_in, harvest_tag)

    if ( hlm_use_luh .eq. itrue ) then
       if(.not. site_in%transition_landuse_from_off_to_on) then
          call GetLanduseTransitionRates(bc_in, site_in%min_allowed_landuse_fraction, &
               site_in%landuse_transition_matrix, site_in%landuse_vector_gt_min)
       else
          call GetInitLanduseTransitionRates(bc_in, site_in%min_allowed_landuse_fraction, &
               site_in%landuse_transition_matrix, site_in%landuse_vector_gt_min)
       endif
    else
       site_in%landuse_transition_matrix(:,:) = 0._r8
    endif

    ! ---------------------------------------------------------------------------------------------
    ! Calculate Disturbance Rates based on the mortality rates just calculated
    ! ---------------------------------------------------------------------------------------------

    ! Recalculate total canopy area prior to resolving the disturbance
    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))
       currentPatch%total_canopy_area = 0._r8
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))   
          if(currentCohort%canopy_layer==1)then
             currentPatch%total_canopy_area = currentPatch%total_canopy_area + currentCohort%c_area
          end if
          currentCohort => currentCohort%taller
       end do
       currentPatch => currentPatch%younger
    end do

    ! get some info needed to determine whether or not to apply land use change
    site_secondaryland_first_exceeding_min = .false.
    if (hlm_use_luh .eq. itrue) then
       call GetLUHStatedata(bc_in, state_vector)
       site_secondaryland_first_exceeding_min =  (state_vector(secondaryland) .gt. site_in%min_allowed_landuse_fraction) &
            .and. (.not. site_in%landuse_vector_gt_min(secondaryland))
    else
       state_vector = current_fates_landuse_state_vector
    end if

    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))   
       
       currentPatch%disturbance_rates(dtype_ifall) = 0.0_r8
       currentPatch%disturbance_rates(dtype_ilog)  = 0.0_r8
       currentPatch%disturbance_rates(dtype_ifire) = 0.0_r8

       dist_rate_ldist_notharvested = 0.0_r8

       ! transition matrix has units of area transitioned per unit area of the whole gridcell per time;
       ! need to change to area transitioned per unit area of that land-use type per time;
       ! because the land use state vector sums to one minus area bareground, need to also divide by that
       ! (or rather, multiply since it is in the denominator of the denominator)
       ! Avoid this calculation to avoid NaN due to division by zero result if luh is not used or applying
       ! to bare ground note that an alternative here might be to use what LUH thinks the state vector
       ! should be instead of what the FATES state vector is, in order to not amplify small deviations
       ! between the two...
       if (hlm_use_luh .eq. itrue .and. currentPatch%land_use_label .gt. nocomp_bareground_land) then
          currentPatch%landuse_transition_rates(1:n_landuse_cats) = min(1._r8, &
               site_in%landuse_transition_matrix(currentPatch%land_use_label,1:n_landuse_cats) &
               * (1._r8 - site_in%area_bareground) / &
               current_fates_landuse_state_vector(currentPatch%land_use_label))
       else
          currentPatch%landuse_transition_rates = 0.0_r8
       end if
       
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))   

          if(currentCohort%canopy_layer == 1)then

             ! Treefall Disturbance Rate.  Only count this for trees, not grasses
             if ( .not. ExemptTreefallDist(currentCohort) ) then
                currentPatch%disturbance_rates(dtype_ifall) = currentPatch%disturbance_rates(dtype_ifall) + &
                     fates_mortality_disturbance_fraction * &
                     min(1.0_r8,currentCohort%dmort)*hlm_freq_day*currentCohort%c_area/currentPatch%area
             end if

             ! Logging Disturbance Rate
             currentPatch%disturbance_rates(dtype_ilog) = currentPatch%disturbance_rates(dtype_ilog) + &
                   min(1.0_r8, currentCohort%lmort_direct +                          & 
                               currentCohort%lmort_collateral +                      &
                               currentCohort%lmort_infra +                           &
                               currentCohort%l_degrad ) *                            &
                               currentCohort%c_area/currentPatch%area

             if(currentPatch%disturbance_rates(dtype_ilog)>1.0) then
                 write(fates_log(),*) 'See luc mortalities:', currentCohort%lmort_direct, &
                     currentCohort%lmort_collateral, currentCohort%lmort_infra, currentCohort%l_degrad
             end if
             
             ! Non-harvested part of the logging disturbance rate
             dist_rate_ldist_notharvested = dist_rate_ldist_notharvested + currentCohort%l_degrad * &
                  currentCohort%c_area/currentPatch%area
             
          endif
          currentCohort => currentCohort%taller
       enddo !currentCohort

       ! for non-closed-canopy areas subject to logging, add an additional increment of area disturbed
       ! equivalent to the fraction logged to account for transfer of interstitial ground area to new secondary lands
       if ( (logging_time .or. site_in%transition_landuse_from_off_to_on .or. site_secondaryland_first_exceeding_min) .and. &
            (currentPatch%area - currentPatch%total_canopy_area) .gt. fates_tiny ) then
          ! The canopy is NOT closed. 

          if (.not. site_in%transition_landuse_from_off_to_on) then
             if(bc_in%hlm_harvest_units == hlm_harvest_carbon) then
                call get_harvest_rate_carbon (currentPatch%land_use_label, bc_in%hlm_harvest_catnames, &
                     bc_in%hlm_harvest_rates, currentPatch%age_since_anthro_disturbance, harvestable_forest_c, &
                     harvest_rate, harvest_tag)
             else
                call get_harvest_rate_area (currentPatch%land_use_label, bc_in%hlm_harvest_catnames, &
                     bc_in%hlm_harvest_rates, current_fates_landuse_state_vector, secondary_young_fraction, &
                     currentPatch%age_since_anthro_disturbance, harvest_rate)
             end if

             ! if the total intended area of secondary lands are less than what we can consider
             ! without having too-small patches, or if that was the case until just now,
             ! then there is special logic
             if (state_vector(secondaryland) .le. site_in%min_allowed_landuse_fraction) then
                harvest_rate = 0._r8
             else if (currentPatch%land_use_label .eq. primaryland .and. .not. &
                  site_in%landuse_vector_gt_min(secondaryland)) then
                harvest_rate = state_vector(secondaryland) / sum(state_vector(:))
             else
                harvest_rate = 0._r8
             end if
          else
             call GetInitLanduseHarvestRate(bc_in, site_in%min_allowed_landuse_fraction, &
                  harvest_rate, site_in%landuse_vector_gt_min)
          endif

          currentPatch%disturbance_rates(dtype_ilog) = currentPatch%disturbance_rates(dtype_ilog) + &
               (currentPatch%area - currentPatch%total_canopy_area) * harvest_rate / currentPatch%area

          ! Non-harvested part of the logging disturbance rate
          dist_rate_ldist_notharvested = dist_rate_ldist_notharvested + &
               (currentPatch%area - currentPatch%total_canopy_area) * harvest_rate / currentPatch%area
       endif

       ! fraction of the logging disturbance rate that is non-harvested
       if (currentPatch%disturbance_rates(dtype_ilog) .gt. nearzero) then
          currentPatch%fract_ldist_not_harvested = dist_rate_ldist_notharvested / &
               currentPatch%disturbance_rates(dtype_ilog)
       endif

       ! Fire Disturbance Rate
       currentPatch%disturbance_rates(dtype_ifire) = currentPatch%frac_burnt

       ! Fires can't burn the whole patch, as this causes /0 errors. 
       if (currentPatch%disturbance_rates(dtype_ifire) > 0.98_r8)then
          msg = 'very high fire areas'//trim(A2S(currentPatch%disturbance_rates(:)))//trim(N2S(currentPatch%frac_burnt))
          call FatesWarn(msg,index=2)
       endif

       ! if the sum of all disturbance rates is such that they will exceed total patch area on this day,
       ! then reduce them all proportionally.
       
       if ( (sum(currentPatch%disturbance_rates(:)) + sum(currentPatch%landuse_transition_rates(1:n_landuse_cats))) .gt. &
            max_daily_disturbance_rate ) then
          tempsum = sum(currentPatch%disturbance_rates(:)) + sum(currentPatch%landuse_transition_rates(1:n_landuse_cats))
          do i_dist = 1,N_DIST_TYPES
             currentPatch%disturbance_rates(i_dist) = max_daily_disturbance_rate * currentPatch%disturbance_rates(i_dist) &
                  / tempsum
          end do
          do i_dist = 1,n_landuse_cats
             currentPatch%landuse_transition_rates(i_dist) = max_daily_disturbance_rate * &
                  currentPatch%landuse_transition_rates(i_dist) / tempsum
          end do
       endif

       currentPatch => currentPatch%younger

    enddo !patch loop

    ! if the area of secondary land has just exceeded the minimum below which we ignore things,
    ! set the flag to keep track of that.
    if ( (state_vector(secondaryland) .gt. site_in%min_allowed_landuse_fraction) .and. &
         (.not. site_in%landuse_vector_gt_min(secondaryland)) ) then
       site_in%landuse_vector_gt_min(secondaryland) = .true.
       write(fates_log(),*) 'setting site_in%landuse_vector_gt_min(secondaryland) = .true.'
       if (debug) then
          currentPatch => site_in%oldest_patch
          do while (associated(currentPatch))
             write(fates_log(),*) 'cpatch area, LU, distrates(ilog): ', currentPatch%area, currentPatch%land_use_label, &
                  currentPatch%nocomp_pft_label, currentPatch%disturbance_rates(dtype_ilog), &
                  currentPatch%area - currentPatch%total_canopy_area
             currentPatch => currentPatch%younger
          end do
       end if
    end if

  end subroutine disturbance_rates

    ! ============================================================================

  subroutine spawn_patches( currentSite, bc_in)
    !
    ! !DESCRIPTION:
    ! In this subroutine, the following happens,
    ! all of which within a complex loop structure of (from outermost to innermost loop),
    ! nocomp-PFT, disturbance type, donor patch land use label, and receiver patch land use label:
    ! 1) the total area disturbed is calculated
    ! 2) a new patch is created
    ! 3) properties are averaged
    ! 4) litter fluxes from fire and mortality are added 
    ! 5) For mortality, plants in existing patch canopy are killed. 
    ! 6) For mortality, Plants in new and existing understorey are killed
    ! 7) For fire, burned plants are killed, and unburned plants are added to new patch. 
    ! 8) New cohorts are added to new patch and sorted. 
    ! 9) New patch is added into linked list
    ! 10) Area checked, and patchno recalculated. 
    !
    ! !USES:

    use EDParamsMod          , only : ED_val_understorey_death, logging_coll_under_frac
    use EDCohortDynamicsMod  , only : terminate_cohorts
    use FatesConstantsMod    , only : rsnbl_math_prec
    use FatesLandUseChangeMod, only : GetLanduseChangeRules
    !
    ! !ARGUMENTS:
    type (ed_site_type), intent(inout) :: currentSite
    type (bc_in_type), intent(in)      :: bc_in
    !
    ! !LOCAL VARIABLES:
    type (fates_patch_type) , pointer :: newPatch
    type (fates_patch_type) , pointer :: currentPatch
    type (fates_cohort_type), pointer :: currentCohort
    type (fates_cohort_type), pointer :: nc
    real(r8) :: site_areadis_primary         ! total area disturbed (to primary forest) in m2 per site per day
    real(r8) :: site_areadis_secondary       ! total area disturbed (to secondary forest) in m2 per site per day
    real(r8) :: patch_site_areadis           ! total area disturbed in m2 per patch per day
    real(r8) :: site_areadis                 ! total site area disturbed in m2 per day
    real(r8) :: age                          ! notional age of this patch in years
    integer  :: el                           ! element loop index
    integer  :: pft                          ! pft loop index
    integer  :: levcan                       ! canopy level
    real(r8) :: leaf_c                       ! leaf carbon [kg]
    real(r8) :: fnrt_c                       ! fineroot carbon [kg]
    real(r8) :: sapw_c                       ! sapwood carbon [kg]
    real(r8) :: store_c                      ! storage carbon [kg]
    real(r8) :: struct_c                     ! structure carbon [kg]
    real(r8) :: total_c                      ! total carbon of plant [kg]
    real(r8) :: leaf_burn_frac               ! fraction of leaves burned in fire
    ! for both woody and grass species
    real(r8) :: leaf_m                       ! leaf mass during partial burn calculations
    integer  :: min_nocomp_pft, max_nocomp_pft, i_nocomp_pft
    integer  :: i_disturbance_type, i_dist2  ! iterators for looping over disturbance types
    integer  :: i_landusechange_receiverpatchlabel  ! iterator for the land use change types
    integer  :: i_donorpatch_landuse_type    ! iterator for the land use change types donor patch
    integer  :: start_receiver_lulabel       ! starting bound for receiver landuse label type loop
    integer  :: end_receiver_lulabel         ! ending bound for receiver landuse label type loop
    real(r8) :: disturbance_rate             ! rate of disturbance being resolved [fraction of patch area / day]
    real(r8) :: oldarea                      ! old patch area prior to disturbance
    logical  :: clearing_matrix(n_landuse_cats,n_landuse_cats)  ! do we clear vegetation when transferring from one LU type to another?
    type (fates_patch_type) , pointer :: buffer_patch, temp_patch, copyPatch, previousPatch
    real(r8) :: nocomp_pft_area_vector(numpft)
    real(r8) :: nocomp_pft_area_vector_filled(numpft)
    real(r8) :: nocomp_pft_area_vector_alt(numpft)
    real(r8) :: newp_area_buffer_frac(numpft) 
    real(r8) :: newp_area_vector(numpft) 
    real(r8) :: max_val
    integer  :: i_land_use_label
    integer  :: i_pft
    real(r8) :: newp_area, area_to_keep, fraction_to_keep
    logical  :: buffer_patch_in_linked_list
    integer  :: n_pfts_by_landuse
    integer  :: which_pft_allowed
    logical  :: buffer_patch_used
    !---------------------------------------------------------------------

    if (hlm_use_nocomp .eq. itrue) then
       min_nocomp_pft = 0
       max_nocomp_pft = numpft
    else
       min_nocomp_pft = fates_unset_int
       max_nocomp_pft = fates_unset_int
    endif

    ! zero the diagnostic disturbance rate fields
    currentSite%disturbance_rates(:,:,:) = 0._r8

    ! get rules for vegetation clearing during land use change
    call GetLanduseChangeRules(clearing_matrix)
    
    ! in the nocomp cases, since every patch has a PFT identity, it can only receive patch area from patches
    ! that have the same identity. In order to allow this, we have this very high level loop over nocomp PFTs
    ! and only do the disturbance for any patches that have that nocomp PFT identity.
    ! If nocomp is not enabled, then this is not much of a loop, it only passes through once.
    nocomp_pft_loop: do i_nocomp_pft = min_nocomp_pft,max_nocomp_pft

       ! we want at the second-outermost loop to go through all disturbance types, because we resolve each of these separately
       disturbance_type_loop: do i_disturbance_type = 1,N_DIST_TYPES

          ! the next loop level is to go through patches that have a specific land-use type. the reason to do this
          ! is because the combination of
          ! disturbance type and donor land-use type uniquly define the land-use type of the receiver patch.
          landuse_donortype_loop: do i_donorpatch_landuse_type = 1, n_landuse_cats

             ! figure out what land use label(s) the receiver patch for disturbance from patches with
             ! this disturbance label and disturbance of this type will have, and set receiver label loop bounds accordingly.

             ! for fire and treefall disturbance, receiver land-use type is whatever the donor land-use type is.
             ! for logging disturbance, receiver land-use type is always secondary lands
             ! for land-use-change disturbance, we need to loop over all possible transition types for land-use-change from
             ! the current land-use type.

             select case(i_disturbance_type)
             case(dtype_ifire)
                start_receiver_lulabel = i_donorpatch_landuse_type
                end_receiver_lulabel = i_donorpatch_landuse_type
             case(dtype_ifall)
                start_receiver_lulabel = i_donorpatch_landuse_type
                end_receiver_lulabel = i_donorpatch_landuse_type
             case(dtype_ilog)
                start_receiver_lulabel = secondaryland
                end_receiver_lulabel = secondaryland
             case(dtype_ilandusechange)
                start_receiver_lulabel = 1  ! this could actually maybe be 2, as primaryland column of matrix should
                ! all be zeros, but leave as 1 for now
                end_receiver_lulabel = n_landuse_cats
             case default
                write(fates_log(),*) 'unknown disturbance mode?'
                write(fates_log(),*) 'i_disturbance_type: ',i_disturbance_type
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end select

             ! next loop level is the set of possible receiver patch land use types.
             ! for disturbance types other than land use change, this is sort of a dummy loop, per the above logic.
             landusechange_receiverpatchlabel_loop: do i_landusechange_receiverpatchlabel = start_receiver_lulabel, &
                  end_receiver_lulabel

                ! now we want to begin resolving all of the disturbance given the above categorical criteria of:
                ! nocomp-PFT, disturbance type, donor patch land use label, and receiver patch land use label.
                ! All of the disturbed area that meets these criteria (if any) will be put into a new patch whose area and
                ! properties are taken from one or more donor patches.

                ! calculate area of disturbed land that meets the above criteria, in this timestep, by summing contributions
                ! from each existing patch.
                currentPatch => currentSite%youngest_patch

                ! this variable site_areadis holds all the newly disturbed area from all patches for all disturbance being
                ! resolved now.
                site_areadis = 0.0_r8

                ! loop over all patches to figure out the total patch area generated as a result of all disturbance being
                ! resolved now.
                patchloop_areadis: do while(associated(currentPatch))

                   cp_nocomp_matches_1_if: if ( hlm_use_nocomp .eq. ifalse .or. &
                        currentPatch%nocomp_pft_label .eq. i_nocomp_pft ) then

                      patchlabel_matches_lutype_if_areadis: if (currentPatch%land_use_label .eq. i_donorpatch_landuse_type) then

                         disturbance_rate = 0.0_r8
                         if ( i_disturbance_type .ne. dtype_ilandusechange) then
                            disturbance_rate = currentPatch%disturbance_rates(i_disturbance_type)
                         else
                            disturbance_rate = currentPatch%landuse_transition_rates(i_landusechange_receiverpatchlabel)
                         endif
                         
                         if(disturbance_rate > (1.0_r8 + rsnbl_math_prec)) then
                            write(fates_log(),*) 'patch disturbance rate > 1 ?',disturbance_rate
                            call currentPatch%Dump()
                            call endrun(msg=errMsg(sourcefile, __LINE__))
                         else if (disturbance_rate > 1.0_r8) then
                            disturbance_rate = 1.0_r8
                         end if

                         ! Only create new patches that have non-negligible amount of land
                         if((currentPatch%area*disturbance_rate) > nearzero ) then

                            site_areadis = site_areadis + currentPatch%area * disturbance_rate

                            ! track disturbance rates to output to history
                            currentSite%disturbance_rates(i_disturbance_type,i_donorpatch_landuse_type,i_landusechange_receiverpatchlabel) = &
                                 currentSite%disturbance_rates(i_disturbance_type,i_donorpatch_landuse_type,i_landusechange_receiverpatchlabel) + &
                                 currentPatch%area * disturbance_rate * AREA_INV
                         end if

                      end if patchlabel_matches_lutype_if_areadis
                   end if cp_nocomp_matches_1_if
                   currentPatch => currentPatch%older
                enddo patchloop_areadis! end loop over patches. sum area disturbed for all patches.

                ! It is possible that no disturbance area was generated
                if ( site_areadis > nearzero) then

                   age = 0.0_r8

                   ! create an empty patch, to absorb newly disturbed area
                   allocate(newPatch)

                   call newPatch%Create(age, site_areadis, i_landusechange_receiverpatchlabel, i_nocomp_pft, &
                                         num_swb, numpft, currentSite%nlevsoil, hlm_current_tod,              &
                                         regeneration_model)

                   ! Initialize the litter pools to zero, these
                   ! pools will be populated by looping over the existing patches
                   ! and transfering in mass
                   do el=1,num_elements
                      call newPatch%litter(el)%InitConditions(init_leaf_fines=0._r8, &
                           init_root_fines=0._r8, &
                           init_ag_cwd=0._r8, &
                           init_bg_cwd=0._r8, &
                           init_seed=0._r8,   &
                           init_seed_germ=0._r8)
                   end do
                   newPatch%tallest  => null()
                   newPatch%shortest => null()

                endif

                ! we now have a new patch and know its area, but it is otherwise empty. Next, we
                ! loop round all the patches that contribute surviving individuals and litter
                ! pools to the new patch.  We only loop the pre-existing patches, so
                ! quit the loop if the current patch is null, and ignore the patch if the patch's categorical variables do not
                ! match those of the outermost set of loops (i.e. the patch's land-use label or nocomp-PFT label
                ! are not what we are resolving right now).

                currentPatch => currentSite%oldest_patch
                patchloop: do while(associated(currentPatch))

                   cp_nocomp_matches_2_if: if ( hlm_use_nocomp .eq. ifalse .or. &
                        currentPatch%nocomp_pft_label .eq. i_nocomp_pft ) then

                      patchlabel_matches_lutype_if: if (currentPatch%land_use_label .eq. i_donorpatch_landuse_type) then


                         ! disturbance_rate is the fraction of the patch's area that is disturbed and donated
                         disturbance_rate = 0.0_r8
                         if ( i_disturbance_type .ne. dtype_ilandusechange) then
                            disturbance_rate = currentPatch%disturbance_rates(i_disturbance_type)
                         else
                            disturbance_rate = currentPatch%landuse_transition_rates(i_landusechange_receiverpatchlabel)
                         endif

                         ! patch_site_areadis is the absolute amount of the patch's area that is disturbed and donated
                         patch_site_areadis = currentPatch%area * disturbance_rate
                         
                         areadis_gt_zero_if: if ( patch_site_areadis > nearzero ) then

                            if(.not.associated(newPatch))then
                               write(fates_log(),*) 'Patch spawning has attempted to point to'
                               write(fates_log(),*) 'an un-allocated patch'
                               call endrun(msg=errMsg(sourcefile, __LINE__))
                            end if

                            ! for the case where the donating patch is not primary, and
                            ! the current disturbance from this patch is non-anthropogenic,
                            ! then we need to average in the time-since-anthropogenic-disturbance
                            ! from the donor patch into that of the receiver patch
                            if ( currentPatch%land_use_label .gt. primaryland .and. &
                                 (i_disturbance_type .lt. dtype_ilog) ) then

                               newPatch%age_since_anthro_disturbance = newPatch%age_since_anthro_disturbance + &
                                     currentPatch%age_since_anthro_disturbance * (patch_site_areadis / site_areadis)
                                  
                            endif

                            ! Transfer the litter existing already in the donor patch to the new patch
                            ! This call will only transfer non-burned litter to new patch
                            ! and burned litter to atmosphere. Thus it is important to zero fuel%frac_burnt when
                            ! fire is not the current disturbance regime.

                            call CopyPatchMeansTimers(currentPatch, newPatch)

                            call TransLitterNewPatch( currentSite, currentPatch, newPatch, patch_site_areadis, i_disturbance_type)

                            ! Transfer in litter fluxes from plants in various contexts of death and destruction
                            select case(i_disturbance_type)
                            case (dtype_ilog)
                               call logging_litter_fluxes(currentSite, currentPatch, &
                                    newPatch, patch_site_areadis,bc_in)

                               ! if transitioning from primary to secondary, then may need to change nocomp pft,
                               ! so tag as having transitioned LU
                               if ( i_disturbance_type .eq. dtype_ilog .and. i_donorpatch_landuse_type .eq. primaryland) then
                                  newPatch%changed_landuse_this_ts = .true.
                               end if
                            case (dtype_ifire)
                               call fire_litter_fluxes(currentSite, currentPatch, &
                                    newPatch, patch_site_areadis,bc_in)
                            case (dtype_ifall)
                               call mortality_litter_fluxes(currentSite, currentPatch, &
                                    newPatch, patch_site_areadis,bc_in)
                            case (dtype_ilandusechange)
                               call landusechange_litter_fluxes(currentSite, currentPatch, &
                                    newPatch, patch_site_areadis,bc_in, &
                                    clearing_matrix(i_donorpatch_landuse_type,i_landusechange_receiverpatchlabel))

                               ! if land use change, then may need to change nocomp pft, so tag as having transitioned LU
                               newPatch%changed_landuse_this_ts = .true.
                            case default
                               write(fates_log(),*) 'unknown disturbance mode?'
                               write(fates_log(),*) 'i_disturbance_type: ',i_disturbance_type
                               call endrun(msg=errMsg(sourcefile, __LINE__))
                            end select

                            ! --------------------------------------------------------------------------
                            ! The newly formed patch from disturbance (newPatch), has now been given
                            ! some litter from dead plants and pre-existing litter from the donor patches.
                            !
                            ! Next, we loop through the cohorts in the donor patch, copy them with
                            ! area modified number density into the new patch, and apply survivorship.
                            ! -------------------------------------------------------------------------

                            currentCohort => currentPatch%shortest
                            cohortloop: do while(associated(currentCohort))

                               allocate(nc)
                               if(hlm_use_planthydro.eq.itrue) call InitHydrCohort(CurrentSite,nc)

                               ! Initialize the PARTEH object and point to the
                               ! correct boundary condition fields
                               nc%prt => null()
                               call InitPRTObject(nc%prt)
                               call nc%InitPRTBoundaryConditions()

                               !  (Keeping as an example)
                               ! Allocate running mean functions
                               !allocate(nc%tveg_lpa)
                               !call nc%tveg_lpa%InitRMean(ema_lpa,init_value=newPatch%tveg_lpa%GetMean())

                               call nc%ZeroValues()

                               ! nc is the new cohort that goes in the disturbed patch (newPatch)... currentCohort
                               ! is the curent cohort that stays in the donor patch (currentPatch)
                               call currentCohort%Copy(nc)

                               !this is the case as the new patch probably doesn't have a closed canopy, and
                               ! even if it does, that will be sorted out in canopy_structure.
                               nc%canopy_layer = 1
                               nc%canopy_layer_yesterday = 1._r8

                               sapw_c   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
                               struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
                               leaf_c   = currentCohort%prt%GetState(leaf_organ, carbon12_element)
                               fnrt_c   = currentCohort%prt%GetState(fnrt_organ, carbon12_element)
                               store_c  = currentCohort%prt%GetState(store_organ, carbon12_element)
                               total_c  = sapw_c + struct_c + leaf_c + fnrt_c + store_c

                               ! survivorship of plants in both the disturbed and undisturbed cohorts depends on what type of
                               ! disturbance is happening.

                               disttype_case: select case(i_disturbance_type)

                               ! treefall mortality is the current disturbance
                               case (dtype_ifall)
                               
                                  in_canopy_if_falldtype: if(currentCohort%canopy_layer == 1)then

                                     ! In the donor patch we are left with fewer trees because the area has decreased
                                     ! the plant density for large trees does not actually decrease in the donor patch
                                     ! because this is the part of the original patch where no trees have actually fallen
                                     ! The diagnostic cmort,bmort,hmort, and frmort  rates have already been saved

                                     currentCohort%n = currentCohort%n * (1.0_r8 - fates_mortality_disturbance_fraction * &
                                          min(1.0_r8,currentCohort%dmort * hlm_freq_day))

                                     nc%n = 0.0_r8      ! kill all of the trees who caused the disturbance.

                                     nc%cmort = nan     ! The mortality diagnostics are set to nan
                                     ! because the cohort should dissappear
                                     nc%hmort = nan
                                     nc%bmort = nan
                                     nc%frmort = nan
                                     nc%smort = nan
                                     nc%asmort = nan
                                     nc%dgmort = nan
                                     nc%lmort_direct     = nan
                                     nc%lmort_collateral = nan
                                     nc%lmort_infra      = nan
                                     nc%l_degrad         = nan

                                  else
                                     ! understory trees
                                     woody_if_falldtype: if( prt_params%woody(currentCohort%pft) == itrue)then


                                        ! Survivorship of undestory woody plants.  Two step process.
                                        ! Step 1:  Reduce current number of plants to reflect the
                                        !          change in area.
                                        !          The number density per square are doesn't change,
                                        !          but since the patch is smaller and cohort counts
                                        !          are absolute, reduce this number.

                                        nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                                        ! because the mortality rate due to impact for the cohorts which
                                        ! had been in the understory and are now in the newly-
                                        ! disturbed patch is very high, passing the imort directly to history
                                        ! results in large numerical errors, on account of the sharply
                                        ! reduced number densities.  so instead pass this info via a
                                        ! site-level diagnostic variable before reducing the number density.

                                        currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) = &
                                             currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) + &
                                             nc%n * ED_val_understorey_death / hlm_freq_day


                                        currentSite%imort_carbonflux(currentCohort%pft) = &
                                             currentSite%imort_carbonflux(currentCohort%pft) + &
                                             (nc%n * ED_val_understorey_death / hlm_freq_day ) * &
                                             total_c * g_per_kg * days_per_sec * years_per_day * ha_per_m2

                                        currentSite%imort_abg_flux(currentCohort%size_class, currentCohort%pft) = &
                                             currentSite%imort_abg_flux(currentCohort%size_class, currentCohort%pft) + &
                                             (nc%n * ED_val_understorey_death / hlm_freq_day ) * &
                                             ( (sapw_c + struct_c + store_c) * prt_params%allom_agb_frac(currentCohort%pft) + &
                                             leaf_c ) * &
                                             g_per_kg * days_per_sec * years_per_day * ha_per_m2


                                        ! Step 2:  Apply survivor ship function based on the understory death fraction
                                        ! remaining of understory plants of those that are knocked over
                                        ! by the overstorey trees dying...
                                        nc%n = nc%n * (1.0_r8 - ED_val_understorey_death)

                                        ! since the donor patch split and sent a fraction of its members
                                        ! to the new patch and a fraction to be preserved in itself,
                                        ! when reporting diagnostic rates, we must carry over the mortality rates from
                                        ! the donor that were applied before the patch split.  Remember this is only
                                        ! for diagnostics.  But think of it this way, the rates are weighted by
                                        ! number density in EDCLMLink, and the number density of this new patch is donated
                                        ! so with the number density must come the effective mortality rates.

                                        nc%cmort            = currentCohort%cmort
                                        nc%hmort            = currentCohort%hmort
                                        nc%bmort            = currentCohort%bmort
                                        nc%frmort           = currentCohort%frmort
                                        nc%smort            = currentCohort%smort
                                        nc%asmort           = currentCohort%asmort
                                        nc%dgmort           = currentCohort%dgmort
                                        nc%dmort            = currentCohort%dmort
                                        nc%lmort_direct     = currentCohort%lmort_direct
                                        nc%lmort_collateral = currentCohort%lmort_collateral
                                        nc%lmort_infra      = currentCohort%lmort_infra

                                        ! understory trees that might potentially be knocked over in the disturbance.
                                        ! The existing (donor) patch should not have any impact mortality, it should
                                        ! only lose cohorts due to the decrease in area.  This is not mortality.
                                        ! Besides, the current and newly created patch sum to unity

                                        currentCohort%n = currentCohort%n * (1._r8 -  patch_site_areadis/currentPatch%area)

                                     else
                                        ! grass is not killed by mortality disturbance events. Just move it into the new patch area.
                                        ! Just split the grass into the existing and new patch structures
                                        nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                                        ! Those remaining in the existing
                                        currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area)

                                        nc%cmort            = currentCohort%cmort
                                        nc%hmort            = currentCohort%hmort
                                        nc%bmort            = currentCohort%bmort
                                        nc%frmort           = currentCohort%frmort
                                        nc%smort            = currentCohort%smort
                                        nc%asmort           = currentCohort%asmort
                                        nc%dgmort           = currentCohort%dgmort
                                        nc%dmort            = currentCohort%dmort
                                        nc%lmort_direct    = currentCohort%lmort_direct
                                        nc%lmort_collateral = currentCohort%lmort_collateral
                                        nc%lmort_infra      = currentCohort%lmort_infra

                                     endif woody_if_falldtype
                                  endif in_canopy_if_falldtype

                               ! Fire is the current disturbance
                               case (dtype_ifire)

                                  ! Number of members in the new patch, before we impose fire survivorship
                                  nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                                  ! loss of individuals from source patch due to area shrinking
                                  currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area)

                                  levcan = currentCohort%canopy_layer

                                  if(levcan==ican_upper) then

                                     ! before changing number densities, track total rate of trees that died
                                     ! due to fire, as well as from each fire mortality term
                                     currentSite%fmort_rate_canopy(currentCohort%size_class, currentCohort%pft) = &
                                          currentSite%fmort_rate_canopy(currentCohort%size_class, currentCohort%pft) + &
                                          nc%n * currentCohort%fire_mort / hlm_freq_day

                                     currentSite%fmort_carbonflux_canopy(currentCohort%pft) = &
                                          currentSite%fmort_carbonflux_canopy(currentCohort%pft) + &
                                          (nc%n * currentCohort%fire_mort) * &
                                          total_c * g_per_kg * days_per_sec * ha_per_m2

                                  else
                                     ! understory
                                     currentSite%fmort_rate_ustory(currentCohort%size_class, currentCohort%pft) = &
                                          currentSite%fmort_rate_ustory(currentCohort%size_class, currentCohort%pft) + &
                                          nc%n * currentCohort%fire_mort / hlm_freq_day

                                     currentSite%fmort_carbonflux_ustory(currentCohort%pft) = &
                                          currentSite%fmort_carbonflux_ustory(currentCohort%pft) + &
                                          (nc%n * currentCohort%fire_mort) * &
                                          total_c * g_per_kg * days_per_sec * ha_per_m2
                                  end if

                                  currentSite%fmort_abg_flux(currentCohort%size_class, currentCohort%pft) = &
                                       currentSite%fmort_abg_flux(currentCohort%size_class, currentCohort%pft) + &
                                       (nc%n * currentCohort%fire_mort) * &
                                       ( (sapw_c + struct_c + store_c) * prt_params%allom_agb_frac(currentCohort%pft) + &
                                       leaf_c ) * &
                                       g_per_kg * days_per_sec * ha_per_m2


                                  currentSite%fmort_rate_cambial(currentCohort%size_class, currentCohort%pft) = &
                                       currentSite%fmort_rate_cambial(currentCohort%size_class, currentCohort%pft) + &
                                       nc%n * currentCohort%cambial_mort / hlm_freq_day
                                  currentSite%fmort_rate_crown(currentCohort%size_class, currentCohort%pft) = &
                                       currentSite%fmort_rate_crown(currentCohort%size_class, currentCohort%pft) + &
                                       nc%n * currentCohort%crownfire_mort / hlm_freq_day

                                  ! loss of individual from fire in new patch.
                                  nc%n = nc%n * (1.0_r8 - currentCohort%fire_mort)

                                  nc%cmort            = currentCohort%cmort
                                  nc%hmort            = currentCohort%hmort
                                  nc%bmort            = currentCohort%bmort
                                  nc%frmort           = currentCohort%frmort
                                  nc%smort            = currentCohort%smort
                                  nc%asmort           = currentCohort%asmort
                                  nc%dgmort           = currentCohort%dgmort
                                  nc%dmort            = currentCohort%dmort
                                  nc%lmort_direct     = currentCohort%lmort_direct
                                  nc%lmort_collateral = currentCohort%lmort_collateral
                                  nc%lmort_infra      = currentCohort%lmort_infra


                                  ! Some of of the leaf mass from living plants has been
                                  ! burned off.  Here, we remove that mass, and
                                  ! tally it in the flux we sent to the atmosphere
                                  if(prt_params%woody(currentCohort%pft) == itrue)then
                                     leaf_burn_frac = currentCohort%fraction_crown_burned
                                  else

                                     ! Grasses determine their fraction of leaves burned here

                                     leaf_burn_frac = currentPatch%fuel%frac_burnt(fuel_classes%live_grass())
                                  endif

                                  ! Perform a check to make sure that spitfire gave
                                  ! us reasonable mortality and burn fraction rates

                                  if( (leaf_burn_frac < 0._r8) .or. &
                                       (leaf_burn_frac > 1._r8) .or. &
                                       (currentCohort%fire_mort < 0._r8) .or. &
                                       (currentCohort%fire_mort > 1._r8)) then
                                     write(fates_log(),*) 'unexpected fire fractions'
                                     write(fates_log(),*) prt_params%woody(currentCohort%pft)
                                     write(fates_log(),*) leaf_burn_frac
                                     write(fates_log(),*) currentCohort%fire_mort
                                     call endrun(msg=errMsg(sourcefile, __LINE__))
                                  end if

                                  do el = 1,num_elements

                                     leaf_m = nc%prt%GetState(leaf_organ, element_list(el))
                                     ! for woody plants burn only leaves
                                     if(int(prt_params%woody(currentCohort%pft)) == itrue)then

                                        leaf_m = nc%prt%GetState(leaf_organ, element_list(el))

                                     else
                                        ! for grasses burn all aboveground tissues
                                        leaf_m = nc%prt%GetState(leaf_organ, element_list(el)) + &
                                             nc%prt%GetState(sapw_organ, element_list(el)) + &
                                             nc%prt%GetState(struct_organ, element_list(el))

                                     endif

                                     currentSite%mass_balance(el)%burn_flux_to_atm = &
                                          currentSite%mass_balance(el)%burn_flux_to_atm + &
                                          leaf_burn_frac * leaf_m * nc%n

                                     ! This diagnostic only tracks
                                     currentSite%flux_diags%elem(el)%burned_liveveg = &
                                          currentSite%flux_diags%elem(el)%burned_liveveg + & 
                                          leaf_burn_frac * leaf_m * nc%n * area_inv
                                     
                                  end do

                                  ! Here the mass is removed from the plant

                                  if(int(prt_params%woody(currentCohort%pft)) == itrue)then
                                     call PRTBurnLosses(nc%prt, leaf_organ, leaf_burn_frac)
                                  else
                                     call PRTBurnLosses(nc%prt, leaf_organ, leaf_burn_frac)
                                     call PRTBurnLosses(nc%prt, sapw_organ, leaf_burn_frac)
                                     call PRTBurnLosses(nc%prt, struct_organ, leaf_burn_frac)
                                  endif

                                  currentCohort%fraction_crown_burned = 0.0_r8
                                  nc%fraction_crown_burned            = 0.0_r8



                               ! Logging is the current disturbance
                               case (dtype_ilog)

                                  ! If this cohort is in the upper canopy. It generated
                                  in_canopy_if_logdtype: if(currentCohort%canopy_layer == 1)then

                                     ! calculate the survivorship of disturbed trees because non-harvested
                                     nc%n = currentCohort%n * currentCohort%l_degrad
                                     ! nc%n            = (currentCohort%l_degrad / (currentCohort%l_degrad + &
                                     !      currentCohort%lmort_direct + currentCohort%lmort_collateral +
                                     !   currentCohort%lmort_infra) ) * &
                                     !      currentCohort%n * patch_site_areadis/currentPatch%area

                                     ! Reduce counts in the existing/donor patch according to the logging rate
                                     currentCohort%n = currentCohort%n * &
                                          (1.0_r8 - min(1.0_r8,(currentCohort%lmort_direct +    &
                                          currentCohort%lmort_collateral + &
                                          currentCohort%lmort_infra + currentCohort%l_degrad)))

                                     nc%cmort            = currentCohort%cmort
                                     nc%hmort            = currentCohort%hmort
                                     nc%bmort            = currentCohort%bmort
                                     nc%frmort           = currentCohort%frmort
                                     nc%smort            = currentCohort%smort
                                     nc%asmort           = currentCohort%asmort
                                     nc%dgmort           = currentCohort%dgmort
                                     nc%dmort            = currentCohort%dmort

                                     ! since these are the ones that weren't logged,
                                     ! set the logging mortality rates as zero
                                     nc%lmort_direct     = 0._r8
                                     nc%lmort_collateral = 0._r8
                                     nc%lmort_infra      = 0._r8

                                  else

                                     ! What to do with cohorts in the understory of a logging generated
                                     ! disturbance patch?

                                     woody_if_logdtype: if(prt_params%woody(currentCohort%pft) == itrue)then


                                        ! Survivorship of undestory woody plants.  Two step process.
                                        ! Step 1:  Reduce current number of plants to reflect the
                                        !          change in area.
                                        !          The number density per square are doesn't change,
                                        !          but since the patch is smaller
                                        !          and cohort counts are absolute, reduce this number.
                                        nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                                        ! because the mortality rate due to impact for the cohorts which had
                                        ! been in the understory and are now in the newly-
                                        ! disturbed patch is very high, passing the imort directly to
                                        ! history results in large numerical errors, on account
                                        ! of the sharply reduced number densities.  so instead pass this info
                                        ! via a site-level diagnostic variable before reducing
                                        ! the number density.
                                        currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) = &
                                             currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) + &
                                             nc%n * currentPatch%fract_ldist_not_harvested * &
                                             logging_coll_under_frac / hlm_freq_day

                                        currentSite%imort_carbonflux(currentCohort%pft) = &
                                             currentSite%imort_carbonflux(currentCohort%pft) + &
                                             (nc%n * currentPatch%fract_ldist_not_harvested * &
                                             logging_coll_under_frac/ hlm_freq_day ) * &
                                             total_c * g_per_kg * days_per_sec * years_per_day * ha_per_m2

                                        currentSite%imort_abg_flux(currentCohort%size_class, currentCohort%pft) = &
                                             currentSite%imort_abg_flux(currentCohort%size_class, currentCohort%pft) + &
                                             (nc%n * currentPatch%fract_ldist_not_harvested * &
                                             logging_coll_under_frac/ hlm_freq_day ) * &
                                             ( ( sapw_c + struct_c + store_c) * prt_params%allom_agb_frac(currentCohort%pft) + &
                                             leaf_c ) * days_per_sec * years_per_day * ha_per_m2

                                        ! Step 2:  Apply survivor ship function based on the understory death fraction

                                        ! remaining of understory plants of those that are knocked
                                        ! over by the overstorey trees dying...
                                        ! LOGGING SURVIVORSHIP OF UNDERSTORY PLANTS IS SET AS A NEW PARAMETER
                                        ! in the fatesparameter files
                                        nc%n = nc%n * (1.0_r8 - &
                                             (1.0_r8-currentPatch%fract_ldist_not_harvested) * logging_coll_under_frac)

                                        ! Step 3: Reduce the number count of cohorts in the
                                        !         original/donor/non-disturbed patch to reflect the area change
                                        currentCohort%n = currentCohort%n * (1._r8 -  patch_site_areadis/currentPatch%area)

                                        nc%cmort            = currentCohort%cmort
                                        nc%hmort            = currentCohort%hmort
                                        nc%bmort            = currentCohort%bmort
                                        nc%frmort           = currentCohort%frmort
                                        nc%smort            = currentCohort%smort
                                        nc%asmort           = currentCohort%asmort
                                        nc%dgmort           = currentCohort%dgmort
                                        nc%dmort            = currentCohort%dmort
                                        nc%lmort_direct     = currentCohort%lmort_direct
                                        nc%lmort_collateral = currentCohort%lmort_collateral
                                        nc%lmort_infra      = currentCohort%lmort_infra

                                     else

                                        ! grass is not killed by mortality disturbance events.
                                        ! Just move it into the new patch area.
                                        ! Just split the grass into the existing and new patch structures
                                        nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                                        ! Those remaining in the existing
                                        currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area)

                                        ! No grass impact mortality imposed on the newly created patch
                                        nc%cmort            = currentCohort%cmort
                                        nc%hmort            = currentCohort%hmort
                                        nc%bmort            = currentCohort%bmort
                                        nc%frmort           = currentCohort%frmort
                                        nc%smort            = currentCohort%smort
                                        nc%asmort           = currentCohort%asmort
                                        nc%dgmort           = currentCohort%dgmort
                                        nc%dmort            = currentCohort%dmort
                                        nc%lmort_direct     = currentCohort%lmort_direct
                                        nc%lmort_collateral = currentCohort%lmort_collateral
                                        nc%lmort_infra      = currentCohort%lmort_infra

                                     endif woody_if_logdtype  ! is/is-not woody

                                  endif in_canopy_if_logdtype ! Select canopy layer

                               ! Land use change is the current disturbance type
                               case (dtype_ilandusechange)

                                  ! Number of members in the new patch, before we impose LUC survivorship
                                  nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                                  ! loss of individuals from source patch due to area shrinking
                                  currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area)

                                  ! now apply survivorship based on the type of landuse transition
                                  if ( clearing_matrix(i_donorpatch_landuse_type,i_landusechange_receiverpatchlabel) ) then
                                     ! kill everything
                                     nc%n = 0._r8
                                  end if

                               case default
                                  write(fates_log(),*) 'unknown disturbance mode?'
                                  write(fates_log(),*) 'i_disturbance_type: ',i_disturbance_type
                                  call endrun(msg=errMsg(sourcefile, __LINE__))
                               end select disttype_case    ! Select disturbance mode

                               ! if some plants in the new temporary cohort survived the transfer to the new patch,
                               ! then put the cohort into the linked list.
                               cohort_n_gt_zero: if (nc%n > 0.0_r8) then
                                 call newPatch%InsertCohort(nc)
                               else
                                  ! sadly, no plants in the cohort survived. on the bright side, we can deallocate their memory.
                                  call nc%FreeMemory()
                                  deallocate(nc, stat=istat, errmsg=smsg)
                                  if (istat/=0) then
                                     write(fates_log(),*) 'dealloc005: fail on deallocate(nc):'//trim(smsg)
                                     call endrun(msg=errMsg(sourcefile, __LINE__))
                                  endif
                               endif cohort_n_gt_zero

                               currentCohort => currentCohort%taller
                            enddo cohortloop
                            call newPatch%ValidateCohorts()

                            call currentPatch%SortCohorts()
                            call currentPatch%ValidateCohorts()

                            !update area of donor patch
                            oldarea = currentPatch%area
                            currentPatch%area = currentPatch%area - patch_site_areadis

                            ! for all disturbance rates that haven't been resolved yet, increase their amount so that
                            ! they are the same amount of gridcell-scale disturbance relative to the original patch size
                            if (i_disturbance_type .lt. N_DIST_TYPES) then
                               do i_dist2 = i_disturbance_type+1,N_DIST_TYPES-1
                                  currentPatch%disturbance_rates(i_dist2) = currentPatch%disturbance_rates(i_dist2) &
                                       * oldarea / currentPatch%area
                               end do
                               do i_dist2 = 1,n_landuse_cats
                                  currentPatch%landuse_transition_rates(i_dist2) = currentPatch%landuse_transition_rates(i_dist2) &
                                       * oldarea / currentPatch%area
                               end do
                            else
                               do i_dist2 = i_landusechange_receiverpatchlabel+1,n_landuse_cats
                                  currentPatch%landuse_transition_rates(i_dist2) = currentPatch%landuse_transition_rates(i_dist2) &
                                       * oldarea / currentPatch%area
                               end do
                            end if

                            ! sort out the cohorts, since some of them may be so small as to need removing.
                            ! the first call to terminate cohorts removes sparse number densities,
                            ! the second call removes for all other reasons (sparse culling must happen
                            ! before fusion)
                            call terminate_cohorts(currentSite, currentPatch, 1,16,bc_in)
                            call fuse_cohorts(currentSite,currentPatch, bc_in)
                            call terminate_cohorts(currentSite, currentPatch, 2,16,bc_in)
                            call currentPatch%SortCohorts()
                            call currentPatch%ValidateCohorts()

                         end if areadis_gt_zero_if   ! if ( newPatch%area > nearzero ) then

                      end if patchlabel_matches_lutype_if

                   end if cp_nocomp_matches_2_if
                   currentPatch => currentPatch%younger

                enddo patchloop ! currentPatch patch loop.

                !*************************/
                !**  INSERT NEW PATCH(ES) INTO LINKED LIST
                !*************************/

                if ( site_areadis .gt. nearzero) then

                   call InsertPatch(currentSite, newPatch)

                   ! sort out the cohorts, since some of them may be so small as to need removing.
                   ! the first call to terminate cohorts removes sparse number densities,
                   ! the second call removes for all other reasons (sparse culling must happen
                   ! before fusion)

                   call terminate_cohorts(currentSite, newPatch, 1,17, bc_in)
                   call fuse_cohorts(currentSite,newPatch, bc_in)
                   call terminate_cohorts(currentSite, newPatch, 2,17, bc_in)
                   call newPatch%SortCohorts()
                   call newPatch%ValidateCohorts()
                endif


                call check_patch_area(currentSite)
                call set_patchno(currentSite,.false.,0)

             end do landusechange_receiverpatchlabel_loop
          end do landuse_donortype_loop
       end do disturbance_type_loop

    end do nocomp_pft_loop

    nocomp_and_luh_if: if ( hlm_use_nocomp .eq. itrue .and. hlm_use_luh .eq. itrue ) then
       ! disturbance has just happened, and now the nocomp PFT identities of the newly-disturbed patches
       ! need to be remapped to those associated with the new land use type.

       ! logic:  loop over land use types. figure out the nocomp PFT fractions for all newly-disturbed patches that have become that land use type.
       ! if the

       lu_loop: do i_land_use_label = n_landuse_cats, 1, -1

          nocomp_pft_area_vector(:) = 0._r8
          nocomp_pft_area_vector_filled(:) = 0._r8
          
          currentPatch => currentSite%oldest_patch
          do while(associated(currentPatch))
             if (currentPatch%changed_landuse_this_ts .and. currentPatch%land_use_label .eq. i_land_use_label) then
                nocomp_pft_area_vector(currentPatch%nocomp_pft_label) = nocomp_pft_area_vector(currentPatch%nocomp_pft_label) + currentPatch%area
                copyPatch => currentPatch
             end if
             currentPatch => currentPatch%younger
          end do

          ! figure out how may PFTs on each land use type. if only 1, then the next calculation is much simpler: we just need to know which PFT is allowed.
          n_pfts_by_landuse = 0
          do i_pft = 1,numpft
             if ( currentSite%area_pft(i_pft,i_land_use_label) .gt. nearzero) then
                n_pfts_by_landuse = n_pfts_by_landuse + 1
                which_pft_allowed = i_pft
             end if
          end do
          if ( n_pfts_by_landuse .ne. 1) then
             which_pft_allowed = fates_unset_int
          endif

          patch_area_to_reallocate_if: if ( sum(nocomp_pft_area_vector(:)) .gt. nearzero ) then
             more_than_1_pft_to_handle_if: if ( n_pfts_by_landuse .ne. 1 ) then
                ! create buffer patch to put all of the pieces carved off of other patches
                allocate(buffer_patch)

                call buffer_patch%Create(0._r8, 0._r8, i_land_use_label, 0, &
                     num_swb, numpft, currentSite%nlevsoil, hlm_current_tod,              &
                     regeneration_model)

                ! Initialize the litter pools to zero
                do el=1,num_elements
                   call buffer_patch%litter(el)%InitConditions(init_leaf_fines=0._r8, &
                        init_root_fines=0._r8, &
                        init_ag_cwd=0._r8, &
                        init_bg_cwd=0._r8, &
                        init_seed=0._r8,   &
                        init_seed_germ=0._r8)
                end do
                buffer_patch%tallest  => null()
                buffer_patch%shortest => null()

                call CopyPatchMeansTimers(copyPatch, buffer_patch)

                ! make a note that this buffer patch has not been put into the linked list
                buffer_patch_in_linked_list = .false.
                buffer_patch_used = .false.

                currentPatch => currentSite%oldest_patch
                do while(associated(currentPatch))
                   if (currentPatch%changed_landuse_this_ts .and. currentPatch%land_use_label .eq. i_land_use_label) then

                      ! Calculate the areas to be given to potentially give to the buffer patch and those to keep in the current patch
                      area_to_keep = currentSite%area_pft(currentPatch%nocomp_pft_label,i_land_use_label) * sum(nocomp_pft_area_vector(:)) - & 
                                     nocomp_pft_area_vector_filled(currentPatch%nocomp_pft_label)
                      newp_area = currentPatch%area - area_to_keep
                      fraction_to_keep = area_to_keep / currentPatch%area

                      if (fraction_to_keep .le. nearzero .or. area_to_keep .lt. rsnbl_math_prec) then
                         ! we don't want any patch area with this PFT identity at all anymore. Fuse it into the buffer patch.

                         currentPatch%nocomp_pft_label = 0
                         if (associated(currentPatch%older)) then
                            previousPatch => currentPatch%older
                         else
                            previousPatch => currentPatch
                         endif

                         call fuse_2_patches(currentSite, currentPatch, buffer_patch)
                         currentPatch => previousPatch

                         buffer_patch_used = .true.

                      elseif ( area_to_keep .ge. rsnbl_math_prec .and. newp_area .ge. rsnbl_math_prec) then
                         ! we have more patch are of this PFT than we want, but we do want to keep some of it.
                         ! we want to split the patch into two here. leave one patch as-is, and put the rest into the buffer patch.

                         allocate(temp_patch)

                         call split_patch(currentSite, currentPatch, temp_patch, fraction_to_keep, newp_area)
                         !
                         temp_patch%nocomp_pft_label = 0

                         call fuse_2_patches(currentSite, temp_patch, buffer_patch)
                         !
                         nocomp_pft_area_vector_filled(currentPatch%nocomp_pft_label) = &
                              nocomp_pft_area_vector_filled(currentPatch%nocomp_pft_label) + currentPatch%area

                         currentPatch%changed_landuse_this_ts = .false.

                         buffer_patch_used = .true.
                      else
                         ! we want to keep all of this patch (and possibly more)

                         nocomp_pft_area_vector_filled(currentPatch%nocomp_pft_label) = &
                              nocomp_pft_area_vector_filled(currentPatch%nocomp_pft_label) + currentPatch%area
                         currentPatch%changed_landuse_this_ts = .false.

                      endif
                   end if

                   currentPatch => currentPatch%younger
                end do

                buffer_patch_used_if: if ( buffer_patch_used ) then
                   ! at this point, lets check that the total patch area remaining to be relabelled equals what we think that it is.
                   if (abs(sum(nocomp_pft_area_vector(:) - nocomp_pft_area_vector_filled(:)) - buffer_patch%area) .gt. rsnbl_math_prec) then
                      write(fates_log(),*) 'midway through patch reallocation and things are already not adding up.', i_land_use_label
                      write(fates_log(),*) currentSite%area_pft(:,i_land_use_label)
                      write(fates_log(),*) '-----'
                      write(fates_log(),*) nocomp_pft_area_vector_filled
                      write(fates_log(),*) '-----'
                      write(fates_log(),*) nocomp_pft_area_vector
                      write(fates_log(),*) '-----'
                      write(fates_log(),*) buffer_patch%area, buffer_patch%land_use_label, buffer_patch%nocomp_pft_label
                      write(fates_log(),*) sum(nocomp_pft_area_vector(:)),  sum(nocomp_pft_area_vector_filled(:)), buffer_patch%area
                      currentPatch => currentSite%oldest_patch
                      do while(associated(currentPatch))
                         write(fates_log(),*) currentPatch%area, currentPatch%land_use_label, currentPatch%nocomp_pft_label
                         currentPatch => currentPatch%younger
                      end do
                      call dump_site(currentSite)
                      call endrun(msg=errMsg(sourcefile, __LINE__))
                   end if

                   ! It's possible that we only need to move all of the buffer into one patch, so first determine what the new patch areas look
                   ! like and compare to the buffer patch area
                   newp_area_vector(:)= (currentSite%area_pft(:,i_land_use_label) * sum(nocomp_pft_area_vector(:))) - nocomp_pft_area_vector_filled(:)
                   newp_area_buffer_frac(:) = newp_area_vector(:) / buffer_patch%area

                   ! Find the maximum value of the vector
                   max_val = maxval(newp_area_buffer_frac)

                   ! If the max value is the only value in the array then loop through the array to find the max value pft index and insert buffer
                   if (abs(sum(newp_area_buffer_frac(:)) - max_val) .le. nearzero) then
                      i_pft = 1
                      do while(.not. buffer_patch_in_linked_list)
                         if (abs(newp_area_buffer_frac(i_pft) - max_val) .le. nearzero) then
                            
                            ! give the buffer patch the intended nocomp PFT label
                            buffer_patch%nocomp_pft_label = i_pft

                            ! track that we have added this patch area
                            nocomp_pft_area_vector_filled(i_pft) = nocomp_pft_area_vector_filled(i_pft) + buffer_patch%area

                            ! put the buffer patch directly into the linked list
                            call InsertPatch(currentSite, buffer_patch)
                            
                            ! Set flag to skip the next pft loop
                            buffer_patch_in_linked_list = .true.
                         end if
                         i_pft = i_pft + 1
                      end do
                   end if

                   ! Now we need to loop through the nocomp PFTs, and split the buffer patch into a set of patches to put back in the linked list
                   ! if not already done so above
                   nocomp_pft_loop_2: do i_pft = 1, numpft

                      ! Check the area fraction to makes sure that this pft should have area.  Also make sure that the buffer patch hasn't been 
                      ! added to the linked list already
                      if ( currentSite%area_pft(i_pft,i_land_use_label) .gt. nearzero .and. .not. buffer_patch_in_linked_list) then

                         ! Slightly complicated way of making sure that the same pfts are subtracted from each other which may help to avoid precision
                         ! errors due to differencing between very large and very small areas
                         nocomp_pft_area_vector_alt(:) = nocomp_pft_area_vector(:)
                         nocomp_pft_area_vector_alt(i_pft) = 0._r8
                         newp_area = (currentSite%area_pft(i_pft,i_land_use_label) * nocomp_pft_area_vector(i_pft)) - nocomp_pft_area_vector_filled(i_pft)
                         newp_area = newp_area + sum(currentSite%area_pft(i_pft,i_land_use_label)*nocomp_pft_area_vector_alt(:))

                         ! Compute area and fraction to keep in buffer
                         area_to_keep = buffer_patch%area - newp_area
                         fraction_to_keep = area_to_keep / buffer_patch%area

                         ! only bother doing this if the new new patch area needed is greater than some tiny amount
                         if ( newp_area .gt. rsnbl_math_prec * 0.01_r8) then

                            if (area_to_keep .gt. rsnbl_math_prec) then

                               ! split buffer patch in two, keeping the smaller buffer patch to put into new patches
                               allocate(temp_patch)

                               call split_patch(currentSite, buffer_patch, temp_patch, fraction_to_keep, newp_area)

                               ! give the new patch the intended nocomp PFT label
                               temp_patch%nocomp_pft_label = i_pft

                               ! track that we have added this patch area
                               nocomp_pft_area_vector_filled(i_pft) = nocomp_pft_area_vector_filled(i_pft) + temp_patch%area

                               ! put the new patch into the linked list
                               call InsertPatch(currentSite, temp_patch)

                            else
                               ! give the buffer patch the intended nocomp PFT label
                               buffer_patch%nocomp_pft_label = i_pft

                               ! track that we have added this patch area
                               nocomp_pft_area_vector_filled(i_pft) = nocomp_pft_area_vector_filled(i_pft) + buffer_patch%area

                               ! put the buffer patch directly into the linked list
                               call InsertPatch(currentSite, buffer_patch)

                               buffer_patch_in_linked_list = .true.

                            end if
                         end if
                      end if
                   end do nocomp_pft_loop_2

                   ! now we want to make sure that either the buffer_patch has zero area (presumably it was never used),
                   ! in which case it should be deallocated, or else it does have area but it has been put into the site
                   ! linked list. if either of those, that means everything worked properly, if not, then something has gone wrong.
                   if ( .not. buffer_patch_in_linked_list) then
                      if (buffer_patch%area .lt. rsnbl_math_prec) then
                         ! here we need to deallocate the buffer patch so that we don't get a memory leak.
                         call buffer_patch%FreeMemory(regeneration_model, numpft)
                         deallocate(buffer_patch, stat=istat, errmsg=smsg)
                         if (istat/=0) then
                            write(fates_log(),*) 'dealloc: fail on deallocate(dp):'//trim(smsg)
                            call endrun(msg=errMsg(sourcefile, __LINE__))
                         endif
                      else
                         write(fates_log(),*) 'Buffer patch still has area and it wasnt put into the linked list'
                         write(fates_log(),*) 'buffer_patch%area', buffer_patch%area
                         write(fates_log(),*) sum(nocomp_pft_area_vector_filled(:)), sum(nocomp_pft_area_vector(:))
                         write(fates_log(),*) sum(nocomp_pft_area_vector_filled(:) - nocomp_pft_area_vector(:))

                         call endrun(msg=errMsg(sourcefile, __LINE__))
                      end if
                   end if
                else
                   ! buffer patch was never even used. deallocate.
                   call buffer_patch%FreeMemory(regeneration_model, numpft)
                   deallocate(buffer_patch, stat=istat, errmsg=smsg)
                   if (istat/=0) then
                      write(fates_log(),*) 'dealloc: fail on deallocate(dp):'//trim(smsg)
                      call endrun(msg=errMsg(sourcefile, __LINE__))
                   endif
                end if buffer_patch_used_if

                ! check that the area we have added is the same as the area we have taken away. if not, crash.
                if ( abs(sum(nocomp_pft_area_vector_filled(:) - nocomp_pft_area_vector(:))) .gt. rsnbl_math_prec) then
                   write(fates_log(),*) 'patch reallocation logic doesnt add up. difference is: ', sum(nocomp_pft_area_vector_filled(:) - nocomp_pft_area_vector(:))
                   write(fates_log(),*) nocomp_pft_area_vector_filled
                   write(fates_log(),*) nocomp_pft_area_vector
                   write(fates_log(),*) i_land_use_label
                   currentPatch => currentSite%oldest_patch
                   do while(associated(currentPatch))
                      write(fates_log(),*) currentPatch%area, currentPatch%land_use_label, currentPatch%nocomp_pft_label
                      currentPatch => currentPatch%younger
                   end do
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
             else
                ! if there is only one PFT allowed on this land use type, then all we need to do is relabel all of the patches that just changed
                ! land use type and let patch fusion take care of the rest.
                currentPatch => currentSite%oldest_patch
                do while(associated(currentPatch))
                   if (currentPatch%changed_landuse_this_ts .and. currentPatch%land_use_label .eq. i_land_use_label) then
                      currentPatch%nocomp_pft_label = which_pft_allowed
                      currentPatch%changed_landuse_this_ts = .false.
                   end if
                   currentPatch => currentPatch%younger
                end do
             endif more_than_1_pft_to_handle_if
          end if patch_area_to_reallocate_if
          call check_patch_area(currentSite)
       end do lu_loop
    else
       ! if not using a configuration where the changed_landuse_this_ts is relevant, loop through all patches and reset it
       currentPatch => currentSite%oldest_patch
       do while(associated(currentPatch))
          currentPatch%changed_landuse_this_ts = .false.
          currentPatch => currentPatch%younger
       end do
    endif nocomp_and_luh_if

    !zero disturbance rate trackers on all patches
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       currentPatch%disturbance_rates(:) = 0._r8
       currentPatch%fract_ldist_not_harvested = 0._r8
       currentPatch => currentPatch%younger
    end do

    return
  end subroutine spawn_patches

  ! -----------------------------------------------------------------------------------------

  subroutine split_patch(currentSite, currentPatch, new_patch, fraction_to_keep, area_to_remove)
    !
    ! !DESCRIPTION:
    !  Split a patch into two patches that are identical except in their areas
    !
    ! !ARGUMENTS:
    type(ed_site_type),intent(inout) :: currentSite
    type(fates_patch_type) , intent(inout), pointer :: currentPatch      ! Donor Patch
    type(fates_patch_type) , intent(inout), pointer :: new_patch         ! New Patch
    real(r8), intent(in)                            :: fraction_to_keep  ! fraction of currentPatch to keep, the rest goes to newpatch
    real(r8), intent(in), optional                  :: area_to_remove     ! area of currentPatch to remove, the rest goes to newpatch
    !
    ! !LOCAL VARIABLES:
    integer  :: el                           ! element loop index
    type (fates_cohort_type), pointer :: nc
    type (fates_cohort_type), pointer :: currentCohort
    integer  :: pft
    real(r8) :: temp_area

    temp_area = 0._r8
    if (present(area_to_remove)) then
       temp_area = area_to_remove
    else
       temp_area = currentPatch%area - (currentPatch%area * fraction_to_keep)
    end if

    ! first we need to make the new patch
    call new_patch%Create(0._r8, temp_area, &
         currentPatch%land_use_label, currentPatch%nocomp_pft_label, &
         num_swb, numpft, currentSite%nlevsoil, hlm_current_tod, &
         regeneration_model)

    ! Initialize the litter pools to zero, these
    ! pools will be populated shortly
    do el=1,num_elements
       call new_patch%litter(el)%InitConditions(init_leaf_fines=0._r8, &
            init_root_fines=0._r8, &
            init_ag_cwd=0._r8, &
            init_bg_cwd=0._r8, &
            init_seed=0._r8,   &
            init_seed_germ=0._r8)
    end do
    new_patch%tallest  => null()
    new_patch%shortest => null()

    call CopyPatchMeansTimers(currentPatch, new_patch)

    call TransLitterNewPatch( currentSite, currentPatch, new_patch, temp_area, 0)

    ! Next, we loop through the cohorts in the donor patch, copy them with
    ! area modified number density into the new-patch, and apply survivorship.
    ! -------------------------------------------------------------------------

    currentCohort => currentPatch%shortest
    do while(associated(currentCohort))

       allocate(nc)
       if(hlm_use_planthydro.eq.itrue) call InitHydrCohort(CurrentSite,nc)

       ! Initialize the PARTEH object and point to the
       ! correct boundary condition fields
       nc%prt => null()
       call InitPRTObject(nc%prt)
       call nc%InitPRTBoundaryConditions()

       !  (Keeping as an example)
       ! Allocate running mean functions
       !allocate(nc%tveg_lpa)
       !call nc%tveg_lpa%InitRMean(ema_lpa,init_value=new_patch%tveg_lpa%GetMean())

       call nc%ZeroValues()

       ! nc is the new cohort that goes in the disturbed patch (new_patch)... currentCohort
       ! is the curent cohort that stays in the donor patch (currentPatch)
       call currentCohort%Copy(nc)

       ! Number of members in the new patch
       nc%n = currentCohort%n * (1._r8 - fraction_to_keep)

       ! loss of individuals from source patch due to area shrinking
       currentCohort%n = currentCohort%n * fraction_to_keep

       call new_patch%InsertCohort(nc)

       currentCohort => currentCohort%taller
    enddo ! currentCohort
    call new_patch%ValidateCohorts()

    call currentPatch%SortCohorts()
    call currentPatch%ValidateCohorts()

    !update area of donor patch
    currentPatch%area = currentPatch%area - temp_area

  end subroutine split_patch

  ! ============================================================================

  subroutine check_patch_area( currentSite )
    !
    ! !DESCRIPTION:
    !  Check to see that total area is not exceeded.  
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout) :: currentSite
    !
    ! !LOCAL VARIABLES:
    real(r8)                     :: areatot
    type(fates_patch_type), pointer :: currentPatch 
    type(fates_patch_type), pointer :: largestPatch
    real(r8)                     :: largest_area
    integer                      :: el
    real(r8)                     :: live_stock
    real(r8)                     :: seed_stock
    real(r8)                     :: litter_stock
    real(r8)                     :: mass_gain
    real(r8), parameter          :: area_error_fail = 1.0e-6_r8
    !---------------------------------------------------------------------

    areatot = 0._r8
    largest_area = 0._r8
    largestPatch => null()
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       areatot = areatot + currentPatch%area
       
       if(currentPatch%area>largest_area) then
          largestPatch => currentPatch
          largest_area = currentPatch%area
       end if
       
       currentPatch => currentPatch%younger
    end do
    
    if ( abs( areatot - area_site ) > nearzero ) then 
       
       if ( abs(areatot-area_site) > area_error_fail ) then
          write(fates_log(),*) 'Patch areas do not sum to 10000 within tolerance'
          write(fates_log(),*) 'Total area: ',areatot,'absolute error: ',areatot-area_site

          currentPatch => currentSite%oldest_patch
          do while(associated(currentPatch))
             write(fates_log(),*) 'area, LU, PFT', currentPatch%area, currentPatch%land_use_label, currentPatch%nocomp_pft_label
             currentPatch => currentPatch%younger
          end do

          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       if(debug) then
          write(fates_log(),*) 'Total patch area precision being fixed, adjusting',(areatot-area_site)
          write(fates_log(),*) 'largest patch. This may have slight impacts on carbon balance.'
       end if
       
       do el = 1,num_elements
           ! This returns the total mass on the patch for the current area [kg]
           call PatchMassStock(largestPatch,el,live_stock,seed_stock,litter_stock)
           
           ! Then we scale the total mass by the added area
           mass_gain = (seed_stock+litter_stock) * &
                 (area_site-areatot)/largestPatch%area

           currentSite%mass_balance(el)%patch_resize_err = &
                 currentSite%mass_balance(el)%patch_resize_err + mass_gain

       end do
       
       largestPatch%area = largestPatch%area + (area_site-areatot)
       
    endif

    return
  end subroutine check_patch_area

  ! ============================================================================

  subroutine TransLitterNewPatch(currentSite,        &
                                 currentPatch,       &
                                 newPatch,           &
                                 patch_site_areadis, &
                                 dist_type)

    ! -----------------------------------------------------------------------------------
    ! 
    ! This routine transfers litter fluxes and rates from a donor patch "currentPatch" into 
    ! the new patch. 
    ! This may include the transfer of existing litter from a patch that burned.
    ! This ROUTINE DOES TRANSFER PARTIALLY BURNED LITTER
    !
    ! Also, note we are not transfering in diagnostics that were calculated
    ! prior to disturbance, because those diagnostics we applied to the patch
    ! before it split, so the diagnostics should reflect those ages and areas.
    !
    ! We do transfer fragmentation fluxes, because we need maintain mass conservation.
    !
    ! We do transfer the seed pool, because we don't currently burn seeds.
    ! Note the seed-pool can decay into the litter pool, where
    ! it can burn.
    !
    ! The "newPatch" is the newly created patch. This patch has already been given
    ! an area which is the sum of disturbed area from a list of donors.  
    ! At this point in the call sequence, we are looping over a list of
    ! donor patches, and transferring over their litter pools which is in units 
    ! kg/m2, we need to make sure that we are conserving mass.
    !
    ! AD = Area of each Donor    [m2]
    ! LD = Litter of each Donor  [kg/m2] 
    !
    ! SUM( AD * LD)  / SUM (AD)    =   SUM( AD*LD/SUM(AD) )
    !
    ! newPatch%area = SUM(AD)   the sum of all donor areas has already been given to newPatch
    ! patch_site_areadis = AD   this is the currently donated area
    !
    ! The fragmentation/decomposition flux from donor patches has 
    ! already occured in existing patches.  However some of their area 
    ! has been carved out for this new patch which is receiving donations.
    ! Lets maintain conservation on that pre-existing mass flux in these 
    ! newly disturbed patches.  Include only the fragmentation flux.
    ! -----------------------------------------------------------------------------------
     
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type)  , intent(in)    :: currentSite        ! site
    type(fates_patch_type) , intent(in)    :: currentPatch       ! Donor patch
    type(fates_patch_type) , intent(inout) :: newPatch           ! New patch
    real(r8)            , intent(in)    :: patch_site_areadis ! Area being donated
                                                              ! by current patch
    integer,              intent(in)    :: dist_type          ! disturbance type

    
    ! locals
    type(site_massbal_type), pointer :: site_mass
    type(litter_type),pointer :: curr_litt  ! litter object for current patch
    type(litter_type),pointer :: new_litt  ! litter object for the new patch
    real(r8) :: remainder_area             ! amount of area remaining in patch after donation
    real(r8) :: burned_mass                ! the mass of litter that was supposed to be provided
                                           ! by the donor, but was burned [kg] 
    real(r8) :: donatable_mass             ! mass of donatable litter [kg]
    real(r8) :: donate_frac                ! the fraction of litter mass sent to the new patch
    real(r8) :: retain_frac                ! the fraction of litter mass retained by the donor patch
    real(r8) :: donate_m2                  ! area normalization for litter mass destined to new patch [m-2]
    real(r8) :: retain_m2                  ! area normalization for litter mass destined to old patch [m-2]
    integer  :: el                         ! element loop counter
    integer  :: c                          ! CWD loop counter
    integer  :: pft                        ! PFT loop counter
    integer  :: dcmpy                      ! Decomposibility loop counter
    integer  :: sl                         ! soil layer loop counter
    real(r8) :: litter_stock0,litter_stock1
    real(r8) :: burn_flux0,burn_flux1
    real(r8) :: error
    real(r8) :: frac_burnt                 ! fraction burnt of current fuel type [0-1]

    do el = 1,num_elements

       site_mass => currentSite%mass_balance(el)
       curr_litt  => currentPatch%litter(el)
       new_litt  => newPatch%litter(el)

       ! Distribute the fragmentation litter flux rates. This is only used for diagnostics
       ! at this point.  Litter fragmentation has already been passed to the output
       ! boundary flux arrays.

       do c = 1,ncwd 
          new_litt%ag_cwd_frag(c) = new_litt%ag_cwd_frag(c) + &
               curr_litt%ag_cwd_frag(c) * patch_site_areadis/newPatch%area
          
          do sl=1,currentSite%nlevsoil
             new_litt%bg_cwd_frag(c,sl) = new_litt%bg_cwd_frag(c,sl) + &
                   curr_litt%bg_cwd_frag(c,sl) * patch_site_areadis/newPatch%area
          end do
       enddo
       
       do dcmpy = 1,ndcmpy

          new_litt%leaf_fines_frag(dcmpy) = new_litt%leaf_fines_frag(dcmpy) + &
               curr_litt%leaf_fines_frag(dcmpy) * patch_site_areadis/newPatch%area
          
          do sl=1,currentSite%nlevsoil
             new_litt%root_fines_frag(dcmpy,sl) = new_litt%root_fines_frag(dcmpy,sl) + &
                   curr_litt%root_fines_frag(dcmpy,sl) * patch_site_areadis/newPatch%area
          end do
          
       enddo

       do pft = 1,numpft
          
          new_litt%seed_decay(pft) = new_litt%seed_decay(pft) + &
               curr_litt%seed_decay(pft)*patch_site_areadis/newPatch%area

          new_litt%seed_germ_decay(pft) = new_litt%seed_germ_decay(pft) + &
               curr_litt%seed_germ_decay(pft)*patch_site_areadis/newPatch%area
          
       end do

       ! -----------------------------------------------------------------------------
       ! Distribute the existing litter that was already in place on the donor
       ! patch.  Some of this burns and is sent to the atmosphere, and some goes to the 
       ! litter stocks of the newly created patch. ALso, some may be retained in the 
       ! donor patch.
       !
       ! This routine handles litter transfer for all types. Note that some of the
       ! transfer may burn. If this routine is being called for a tree-fall
       ! or logging disturbance, it is assumed that the burned_masses should come
       ! out to zero.
       ! -----------------------------------------------------------------------------

       ! If/when sending litter fluxes to the old patch, we divide the total 
       ! mass sent to that patch, by the area it will have remaining
       ! after it donates area.
       ! i.e. subtract the area it is donating.
       
       remainder_area = currentPatch%area - patch_site_areadis

       ! Calculate the fraction of litter to be retained versus donated
       ! vis-a-vis the new and donor patch

       retain_frac = (1.0_r8-existing_litt_localization) * &
             remainder_area/(newPatch%area+remainder_area)
       donate_frac = 1.0_r8-retain_frac
        
       if(remainder_area > rsnbl_math_prec) then
           retain_m2 = retain_frac/remainder_area
           donate_m2 = (1.0_r8-retain_frac)/newPatch%area
       else
           retain_m2 = 0._r8
           donate_m2 = 1.0_r8/newPatch%area
       end if


       if (debug) then
          burn_flux0    = site_mass%burn_flux_to_atm
          litter_stock0 = curr_litt%GetTotalLitterMass()*currentPatch%area + & 
                          new_litt%GetTotalLitterMass()*newPatch%area
       end if
       
       do c = 1,ncwd
         frac_burnt = 0.0_r8
         if (dist_type == dtype_ifire .and. currentPatch%fire == 1) then
            frac_burnt = currentPatch%fuel%frac_burnt(c)
         end if 
             
          ! Transfer above ground CWD
          donatable_mass     = curr_litt%ag_cwd(c) * patch_site_areadis * &
                               (1._r8 - frac_burnt)

          burned_mass        = curr_litt%ag_cwd(c) * patch_site_areadis * &
                               frac_burnt
 
          new_litt%ag_cwd(c) = new_litt%ag_cwd(c) + donatable_mass*donate_m2
          curr_litt%ag_cwd(c) = curr_litt%ag_cwd(c) + donatable_mass*retain_m2

          site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass
             
          ! Transfer below ground CWD (none burns)
          
          do sl = 1,currentSite%nlevsoil
             donatable_mass         = curr_litt%bg_cwd(c,sl) * patch_site_areadis
             new_litt%bg_cwd(c,sl)  = new_litt%bg_cwd(c,sl) + donatable_mass*donate_m2
             curr_litt%bg_cwd(c,sl) = curr_litt%bg_cwd(c,sl) + donatable_mass*retain_m2
          end do
          
       enddo
       
       frac_burnt = 0.0_r8
       if (dist_type == dtype_ifire .and. currentPatch%fire == 1) then
         frac_burnt = currentPatch%fuel%frac_burnt(fuel_classes%dead_leaves())
      end if 
             
       do dcmpy=1,ndcmpy

           ! Transfer leaf fines
           donatable_mass           = curr_litt%leaf_fines(dcmpy) * patch_site_areadis * &
                                      (1._r8 - frac_burnt)

           burned_mass              = curr_litt%leaf_fines(dcmpy) * patch_site_areadis * &
                                       frac_burnt

           new_litt%leaf_fines(dcmpy) = new_litt%leaf_fines(dcmpy) + donatable_mass*donate_m2
           curr_litt%leaf_fines(dcmpy) = curr_litt%leaf_fines(dcmpy) + donatable_mass*retain_m2
           
           site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass
           
           ! Transfer root fines (none burns)
           do sl = 1,currentSite%nlevsoil
               donatable_mass = curr_litt%root_fines(dcmpy,sl) * patch_site_areadis             
               new_litt%root_fines(dcmpy,sl) = new_litt%root_fines(dcmpy,sl) + donatable_mass*donate_m2
               curr_litt%root_fines(dcmpy,sl) = curr_litt%root_fines(dcmpy,sl) + donatable_mass*retain_m2
          end do
          
       end do
             
       do pft = 1,numpft

          ! Transfer seeds (currently we don't burn seeds)
          donatable_mass = curr_litt%seed(pft) * patch_site_areadis

          new_litt%seed(pft) = new_litt%seed(pft) + donatable_mass * donate_m2
          curr_litt%seed(pft) = curr_litt%seed(pft) + donatable_mass * retain_m2
          
          donatable_mass = curr_litt%seed_germ(pft) * patch_site_areadis

          new_litt%seed_germ(pft) = new_litt%seed_germ(pft) + donatable_mass * donate_m2
          curr_litt%seed_germ(pft) = curr_litt%seed_germ(pft) + donatable_mass * retain_m2
          
       enddo

       ! --------------------------------------------------------------------------
       ! Mass conservation check, set debug=.true. if mass imbalances in 
       ! EDMainMod start triggering.
       ! --------------------------------------------------------------------------
       if (debug) then
          burn_flux1    = site_mass%burn_flux_to_atm
          litter_stock1 = curr_litt%GetTotalLitterMass()*remainder_area + & 
                          new_litt%GetTotalLitterMass()*newPatch%area
          error = (litter_stock1 - litter_stock0) + (burn_flux1-burn_flux0)
          if(abs(error)>1.e-8_r8) then
             write(fates_log(),*) 'non trivial carbon mass balance error in litter transfer'
             write(fates_log(),*) 'abs error: ',error
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       end if


    end do

    return
  end subroutine TransLitterNewPatch

  ! ============================================================================

  subroutine fire_litter_fluxes(currentSite, currentPatch, &
       newPatch, patch_site_areadis, bc_in)
    !
    ! !DESCRIPTION:
    !  CWD pool burned by a fire. 
    !  Carbon going from burned trees into CWD pool
    !  Burn parts of trees that don't die in fire
    !  Burn live grasses and kill them. 
    !  Note: The number density of living plants in the donating patch (currentPatch)
    !        has not been scaled down by area yet. That happens after this routine.

    !
    ! !USES:
    use SFParamsMod,          only : SF_VAL_CWD_FRAC
    !
    ! !ARGUMENTS:
    type(ed_site_type)  , intent(inout), target :: currentSite
    type(fates_patch_type) , intent(inout), target :: currentPatch   ! Donor Patch
    type(fates_patch_type) , intent(inout), target :: newPatch   ! New Patch
    real(r8)            , intent(in)            :: patch_site_areadis ! Area being donated
    type(bc_in_type)    , intent(in)            :: bc_in
    
    !
    ! !LOCAL VARIABLES:

    type(fates_cohort_type), pointer      :: currentCohort
    type(litter_type), pointer         :: new_litt
    type(litter_type), pointer         :: curr_litt
    type(site_massbal_type), pointer   :: site_mass
    type(elem_diag_type), pointer      :: elflux_diags

    real(r8) :: donatable_mass       ! non-burned litter mass provided by the donor [kg]
                                     ! some may or may not be retained by the donor
    real(r8) :: burned_mass          ! the mass of litter that was supposed to be provided
                                     ! by the donor, but was burned [kg]
    real(r8) :: remainder_area       ! current patch's remaining area after donation [m2]
    real(r8) :: retain_frac          ! the fraction of litter mass retained by the donor patch
    real(r8) :: bcroot               ! amount of below ground coarse root per cohort kg
    real(r8) :: bstem                ! amount of above ground stem biomass per cohort kg
    real(r8) :: leaf_burn_frac       ! fraction of leaves burned 
    real(r8) :: leaf_m               ! leaf mass [kg]
    real(r8) :: fnrt_m               ! fineroot mass [kg]
    real(r8) :: sapw_m               ! sapwood mass [kg]
    real(r8) :: store_m              ! storage mass [kg]
    real(r8) :: struct_m             ! structure mass [kg]
    real(r8) :: repro_m              ! Reproductive mass (seeds/flowers) [kg]
    real(r8) :: num_dead_trees       ! total number of dead trees passed in with the burn area
    real(r8) :: num_live_trees       ! total number of live trees passed in with the burn area
    real(r8) :: donate_m2            ! area normalization for litter mass destined to new patch [m-2]
    real(r8) :: retain_m2            ! area normalization for litter mass staying in donor patch [m-2]
    real(r8) :: dcmpy_frac           ! fraction of mass going to each decomposability partition
    integer  :: el                   ! element loop index
    integer  :: sl                   ! soil layer index
    integer  :: c                    ! loop index for coarse woody debris pools
    integer  :: pft                  ! loop index for plant functional types
    integer  :: dcmpy                ! loop index for decomposability pool
    integer  :: element_id           ! parteh compatible global element index
    real(r8) :: SF_val_CWD_frac_adj(4) !Updated wood partitioning to CWD based on dbh
    !---------------------------------------------------------------------

    ! Only do this if there was a fire in this actual patch. 
    if ( currentPatch%fire  ==  ifalse ) return

    ! If plant hydraulics are turned on, account for water leaving the plant-soil
    ! mass balance through the dead trees
    if (hlm_use_planthydro == itrue) then
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))
          num_dead_trees  = (currentCohort%fire_mort * &
                currentCohort%n*patch_site_areadis/currentPatch%area)
          call AccumulateMortalityWaterStorage(currentSite,currentCohort,num_dead_trees)
          currentCohort => currentCohort%taller
       end do
    end if


    ! If/when sending litter fluxes to the donor patch, we divide the total 
    ! mass sent to that patch, by the area it will have remaining
    ! after it donates area.
    ! i.e. subtract the area it is donating.
    
    remainder_area = currentPatch%area - patch_site_areadis
   
    ! Calculate the fraction of litter to be retained versus donated
    ! vis-a-vis the new and donor patch (if the area remaining
    ! in the original donor patch is small, don't bother 
    ! retaining anything.)
    retain_frac = (1.0_r8-burn_localization) * &
          remainder_area/(newPatch%area+remainder_area)

    if(remainder_area > rsnbl_math_prec) then
        retain_m2 = retain_frac/remainder_area
        donate_m2 = (1.0_r8-retain_frac)/newPatch%area
    else
        retain_m2 = 0._r8
        donate_m2 = 1.0_r8/newPatch%area
    end if

    do el = 1,num_elements
       
       element_id = element_list(el)
       site_mass  => currentSite%mass_balance(el)
       elflux_diags => currentSite%flux_diags%elem(el)
       curr_litt  => currentPatch%litter(el)      ! Litter pool of "current" patch
       new_litt   => newPatch%litter(el)          ! Litter pool of "new" patch
       
       ! -----------------------------------------------------------------------------
       ! PART 1) Handle mass fluxes associated with plants that died in the fire. This
       ! includes transfer of non burned plant material to litter, and the burned
       ! part to the atmosphere.
       ! ------------------------------------------------------------------------------

       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))
          
             pft = currentCohort%pft

             ! Number of trees that died because of the fire, per m2 of ground. 
             ! Divide their litter into the four litter streams, and spread 
             ! across ground surface. 
             ! -----------------------------------------------------------------------

             fnrt_m   = currentCohort%prt%GetState(fnrt_organ, element_id)
             store_m  = currentCohort%prt%GetState(store_organ, element_id)
             repro_m  = currentCohort%prt%GetState(repro_organ, element_id)
          
             if (prt_params%woody(currentCohort%pft) == itrue) then
                ! Assumption: for woody plants fluxes from deadwood and sapwood go together in CWD pool
                leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id)
                sapw_m          = currentCohort%prt%GetState(sapw_organ,element_id)
                struct_m        = currentCohort%prt%GetState(struct_organ,element_id)
             else
                ! for non-woody plants all stem fluxes go into the same leaf litter pool
                leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id) + &
                     currentCohort%prt%GetState(sapw_organ,element_id) + &
                     currentCohort%prt%GetState(struct_organ,element_id)
                sapw_m          = 0._r8
                struct_m        = 0._r8
             end if


             ! Absolute number of dead trees being transfered in with the donated area
             num_dead_trees = (currentCohort%fire_mort*currentCohort%n * &
                               patch_site_areadis/currentPatch%area)

             ! Contribution of dead trees to leaf litter
             donatable_mass = num_dead_trees * (leaf_m+repro_m) * &
                              (1.0_r8-currentCohort%fraction_crown_burned)

             ! Contribution of dead trees to leaf burn-flux
             burned_mass  = num_dead_trees * (leaf_m+repro_m) * currentCohort%fraction_crown_burned

             do dcmpy=1,ndcmpy
                 dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)
                 new_litt%leaf_fines(dcmpy) = new_litt%leaf_fines(dcmpy) + &
                                              donatable_mass*donate_m2*dcmpy_frac
                 curr_litt%leaf_fines(dcmpy) = curr_litt%leaf_fines(dcmpy) + &
                                               donatable_mass*retain_m2*dcmpy_frac
             end do

             site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass

             
             
             call set_root_fraction(currentSite%rootfrac_scr, pft, currentSite%zi_soil, &
                  bc_in%max_rooting_depth_index_col)

             ! Contribution of dead trees to root litter (no root burn flux to atm)
             do dcmpy=1,ndcmpy
                 dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
                 do sl = 1,currentSite%nlevsoil
                     donatable_mass = num_dead_trees * (fnrt_m+store_m) * currentSite%rootfrac_scr(sl)
                     new_litt%root_fines(dcmpy,sl) = new_litt%root_fines(dcmpy,sl) + &
                                                     donatable_mass*donate_m2*dcmpy_frac
                     curr_litt%root_fines(dcmpy,sl) = curr_litt%root_fines(dcmpy,sl) + &
                                                      donatable_mass*retain_m2*dcmpy_frac
                 end do
             end do

             ! Track as diagnostic fluxes
             elflux_diags%surf_fine_litter_input(pft) = &
                  elflux_diags%surf_fine_litter_input(pft) + &
                  num_dead_trees * (leaf_m+repro_m) * (1.0_r8-currentCohort%fraction_crown_burned)

             elflux_diags%root_litter_input(pft) = &
                  elflux_diags%root_litter_input(pft) + &
                  (fnrt_m + store_m) * num_dead_trees

             ! coarse root biomass per tree
             bcroot = (sapw_m + struct_m) * (1.0_r8 - prt_params%allom_agb_frac(pft) )

             ! below ground coarse woody debris from burned trees

             !adjust the how wood is partitioned between the cwd classes based on cohort dbh
	     call adjust_SF_CWD_frac(currentCohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)

             do c = 1,ncwd
                do sl = 1,currentSite%nlevsoil
                   donatable_mass =  num_dead_trees * SF_val_CWD_frac_adj(c) * &
                         bcroot * currentSite%rootfrac_scr(sl)

                   new_litt%bg_cwd(c,sl) = new_litt%bg_cwd(c,sl) + &
                         donatable_mass * donate_m2
                   curr_litt%bg_cwd(c,sl) = curr_litt%bg_cwd(c,sl) + &
                         donatable_mass * retain_m2

                   ! track diagnostics
                   elflux_diags%cwd_bg_input(c) = &
                        elflux_diags%cwd_bg_input(c) + &
                        donatable_mass
                enddo
             end do

             ! stem biomass per tree
             bstem  = (sapw_m + struct_m) * prt_params%allom_agb_frac(pft)

             ! Above ground coarse woody debris from twigs and small branches
             ! a portion of this pool may burn
             do c = 1,ncwd
                 donatable_mass = num_dead_trees * SF_val_CWD_frac_adj(c) * bstem
                 if (c == 1 .or. c == 2) then
                      donatable_mass = donatable_mass * (1.0_r8-currentCohort%fraction_crown_burned)
                      burned_mass = num_dead_trees * SF_val_CWD_frac_adj(c) * bstem * &
                      currentCohort%fraction_crown_burned
                      site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass
                endif
                new_litt%ag_cwd(c) = new_litt%ag_cwd(c) + donatable_mass * donate_m2
                curr_litt%ag_cwd(c) = curr_litt%ag_cwd(c) + donatable_mass * retain_m2
                elflux_diags%cwd_ag_input(c) = elflux_diags%cwd_ag_input(c) + donatable_mass
             enddo


            currentCohort => currentCohort%taller
        enddo
    end do
    
    return
  end subroutine fire_litter_fluxes

  ! ============================================================================

  subroutine mortality_litter_fluxes(currentSite, currentPatch, &
       newPatch, patch_site_areadis,bc_in)
    !
    ! !DESCRIPTION:
    ! Carbon going from mortality associated with disturbance into CWD pools. 
    ! By "associated with disturbance", this includes tree death that
    ! forms gaps, as well as tree death due to impacts from those trees.
    !
    ! We calculate the fraction of litter to be retained versus donated
    ! vis-a-vis the new and donor patch. At this step, we have not
    ! yet removed the area from the pre-existing patch (currentPatch),
    ! so we pre-compute "remainder_area", which is the soon-to-be
    ! area of the patch once disturbance is completed.
    !
    ! !USES:
    use EDParamsMod,  only : ED_val_understorey_death
    use SFParamsMod,  only : SF_val_cwd_frac
    !
    ! !ARGUMENTS:
    type(ed_site_type)  , intent(inout), target :: currentSite 
    type(fates_patch_type) , intent(inout), target :: currentPatch
    type(fates_patch_type) , intent(inout), target :: newPatch
    real(r8)            , intent(in)            :: patch_site_areadis
    type(bc_in_type)    , intent(in)            :: bc_in
    !
    ! !LOCAL VARIABLES:
    type(fates_cohort_type), pointer      :: currentCohort
    type(litter_type), pointer         :: new_litt
    type(litter_type), pointer         :: curr_litt
    type(site_massbal_type), pointer   :: site_mass
    type(elem_diag_type), pointer      :: elflux_diags

    real(r8) :: remainder_area       ! amount of area remaining in patch after donation
    real(r8) :: num_dead
    real(r8) :: donatable_mass       ! mass of donatable litter [kg]
    real(r8) :: leaf_m               ! leaf mass [kg]
    real(r8) :: fnrt_m               ! fineroot mass [kg]
    real(r8) :: sapw_m               ! sapwood mass [kg]
    real(r8) :: store_m              ! storage mass [kg]
    real(r8) :: struct_m             ! structure mass [kg]
    real(r8) :: repro_m              ! reproductive mass [kg]
    real(r8) :: retain_frac          ! Fraction of mass to be retained
    real(r8) :: donate_frac          ! Fraction of mass to be donated
    real(r8) :: donate_m2            ! area normalization for litter mass destined to new patch [m-2]
    real(r8) :: retain_m2            ! area normalization for litter mass destined to old patch [m-2]
    real(r8) :: ag_wood              ! Total above ground mass in wood [kg]
    real(r8) :: bg_wood              ! Total bg mass in wood [kg]
    real(r8) :: seed_mass            ! Total seed mass generated from storage death [kg]
    integer  :: pft                  ! plant functional type index
    integer  :: dcmpy                ! decomposability index
    integer  :: c                    ! coarse woody debris pool index
    integer  :: el                   ! element loop index
    integer  :: sl                   ! soil layer index
    integer  :: element_id           ! parteh compatible global element index
    real(r8) :: dcmpy_frac           ! decomposability fraction
    real(r8) :: SF_val_CWD_frac_adj(4) !Updated wood partitioning to CWD based on dbh
    !---------------------------------------------------------------------

    remainder_area = currentPatch%area - patch_site_areadis
    
    retain_frac = (1.0_r8-treefall_localization) * &
         remainder_area/(newPatch%area+remainder_area)
    donate_frac = 1.0_r8-retain_frac
    
    if(remainder_area > rsnbl_math_prec) then
       retain_m2 = retain_frac/remainder_area
       donate_m2 = (1.0_r8-retain_frac)/newPatch%area
    else
       retain_m2 = 0._r8
       donate_m2 = 1._r8/newPatch%area
    end if


    do el = 1,num_elements
       
       element_id = element_list(el)
       site_mass  => currentSite%mass_balance(el)
       elflux_diags => currentSite%flux_diags%elem(el)
       curr_litt  => currentPatch%litter(el)   ! Litter pool of "current" patch
       new_litt   => newPatch%litter(el)       ! Litter pool of "new" patch

       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))       

          pft = currentCohort%pft
   
          fnrt_m   = currentCohort%prt%GetState(fnrt_organ, element_id)
          store_m  = currentCohort%prt%GetState(store_organ, element_id)
          repro_m  = currentCohort%prt%GetState(repro_organ, element_id)

          if (prt_params%woody(currentCohort%pft) == itrue) then
             ! Assumption: for woody plants fluxes from deadwood and sapwood go together in CWD pool
             leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id)
             sapw_m          = currentCohort%prt%GetState(sapw_organ,element_id)
             struct_m        = currentCohort%prt%GetState(struct_organ,element_id)
          else
             ! for non-woody plants all stem fluxes go into the same leaf litter pool
             leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id) + &
                     currentCohort%prt%GetState(sapw_organ,element_id) + &
                     currentCohort%prt%GetState(struct_organ,element_id)
             sapw_m          = 0._r8
             struct_m        = 0._r8
          end if

          if(currentCohort%canopy_layer == 1)then

             ! Upper canopy trees. The total dead is based on their disturbance
             ! generating mortality rate.
             
             num_dead = currentCohort%n * min(1.0_r8,currentCohort%dmort * &
                   hlm_freq_day * fates_mortality_disturbance_fraction)
             
          elseif(prt_params%woody(pft) == itrue) then
             
             ! Understorey trees. The total dead is based on their survivorship
             ! function, and the total area of disturbance.
             
             num_dead = ED_val_understorey_death * currentCohort%n * &
                   (patch_site_areadis/currentPatch%area) 

          else
             
             ! The only thing left is uderstory grasses. These guys aren't
             ! killed by tree-fall disturbance events.

             num_dead = 0._r8
             
          end if

          ! Update water balance by removing dead plant water
          ! but only do this once (use the carbon element id)
          if( (element_id == carbon12_element) .and. &
              (hlm_use_planthydro == itrue) ) then
              call AccumulateMortalityWaterStorage(currentSite,currentCohort, num_dead)
          end if
          
          ! Transfer leaves of dying trees to leaf litter (includes seeds too)
          do dcmpy=1,ndcmpy
              dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)
              new_litt%leaf_fines(dcmpy) = new_litt%leaf_fines(dcmpy) + &
                    num_dead*(leaf_m+repro_m)*donate_m2*dcmpy_frac
              
              curr_litt%leaf_fines(dcmpy) = curr_litt%leaf_fines(dcmpy) + &
                    num_dead*(leaf_m+repro_m)*retain_m2*dcmpy_frac
          end do
                 
          ! Pre-calculate Structural and sapwood, below and above ground, total mass [kg]
          ag_wood = num_dead * (struct_m + sapw_m) * prt_params%allom_agb_frac(pft)
          bg_wood = num_dead * (struct_m + sapw_m) * (1.0_r8-prt_params%allom_agb_frac(pft))
          
          call set_root_fraction(currentSite%rootfrac_scr, pft, currentSite%zi_soil, &
               bc_in%max_rooting_depth_index_col)

          ! Adjust how wood is partitioned between the cwd classes based on cohort dbh
	  call adjust_SF_CWD_frac(currentCohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)

          do c=1,ncwd
             
             ! Transfer wood of dying trees to AG CWD pools
             new_litt%ag_cwd(c) = new_litt%ag_cwd(c) + ag_wood * &
                    SF_val_CWD_frac_adj(c) * donate_m2

             curr_litt%ag_cwd(c) = curr_litt%ag_cwd(c) + ag_wood * &
                   SF_val_CWD_frac_adj(c) * retain_m2
             
             ! Transfer wood of dying trees to BG CWD pools
             do sl = 1,currentSite%nlevsoil
                new_litt%bg_cwd(c,sl) = new_litt%bg_cwd(c,sl) + bg_wood * &
                       currentSite%rootfrac_scr(sl) * SF_val_CWD_frac_adj(c) * &
                       donate_m2

                curr_litt%bg_cwd(c,sl) = curr_litt%bg_cwd(c,sl) + bg_wood * &
                      currentSite%rootfrac_scr(sl) * SF_val_CWD_frac_adj(c) * &
                      retain_m2
             end do
          end do

          ! Transfer fine roots of dying trees to below ground litter pools
          do dcmpy=1,ndcmpy
              dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
              do sl=1,currentSite%nlevsoil
                  new_litt%root_fines(dcmpy,sl) = new_litt%root_fines(dcmpy,sl) + &
                        num_dead * currentSite%rootfrac_scr(sl) * &
                        (fnrt_m + store_m*(1.0_r8-EDPftvarcon_inst%allom_frbstor_repro(pft))) * &
                        donate_m2 * dcmpy_frac
                  
                  curr_litt%root_fines(dcmpy,sl) = curr_litt%root_fines(dcmpy,sl) + &
                        num_dead * currentSite%rootfrac_scr(sl) * &
                        (fnrt_m + store_m*(1.0_r8-EDPftvarcon_inst%allom_frbstor_repro(pft))) * &
                        retain_m2 * dcmpy_frac
              end do
          end do
              
          ! Transfer some of the storage that is shunted to reproduction
          ! upon death, to the seed-pool. This is was designed for grasses,
          ! but it is possible that some trees may utilize this behavior too

          seed_mass = num_dead * store_m * EDPftvarcon_inst%allom_frbstor_repro(pft)

          ! SEED DISTRIBUTION IS BREAKING MASS CONSERVATION RIGHT NOW...
!          call DistributeSeeds(currentSite,seed_mass,el,pft)

          new_litt%seed(pft) = new_litt%seed(pft) + seed_mass * donate_m2
          curr_litt%seed(pft) = curr_litt%seed(pft) + seed_mass * retain_m2
          
          ! track diagnostic fluxes
          do c=1,ncwd
             elflux_diags%cwd_ag_input(c) = & 
                  elflux_diags%cwd_ag_input(c) + SF_val_CWD_frac_adj(c) * ag_wood
             
             elflux_diags%cwd_bg_input(c) = &
                  elflux_diags%cwd_bg_input(c) + SF_val_CWD_frac_adj(c) * bg_wood
          end do

          elflux_diags%surf_fine_litter_input(pft) = elflux_diags%surf_fine_litter_input(pft) + &
               num_dead*(leaf_m + repro_m)

          elflux_diags%root_litter_input(pft) = elflux_diags%root_litter_input(pft) + & 
               num_dead * (fnrt_m + store_m*(1.0_r8-EDPftvarcon_inst%allom_frbstor_repro(pft)))
          

          
          currentCohort => currentCohort%taller      
       enddo !currentCohort         

    enddo


    return
  end subroutine mortality_litter_fluxes

    ! ============================================================================

  subroutine landusechange_litter_fluxes(currentSite, currentPatch, &
       newPatch, patch_site_areadis, bc_in, &
       clearing_matrix_element)
    !
    ! !DESCRIPTION:
    !  CWD pool from land use change.
    !  Carbon going from felled trees into CWD pool
    !  Either kill everything or nothing on disturbed land, depending on clearing matrix
    !
    ! !USES:
    use SFParamsMod,          only : SF_VAL_CWD_FRAC
    !
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout), target :: currentSite
    type(fates_patch_type) , intent(inout), target :: currentPatch   ! Donor Patch
    type(fates_patch_type) , intent(inout), target :: newPatch   ! New Patch
    real(r8)               , intent(in)            :: patch_site_areadis ! Area being donated
    type(bc_in_type)       , intent(in)            :: bc_in
    logical                , intent(in)            :: clearing_matrix_element ! whether or not to clear vegetation

    !
    ! !LOCAL VARIABLES:

    type(fates_cohort_type), pointer      :: currentCohort
    type(litter_type), pointer         :: new_litt
    type(litter_type), pointer         :: curr_litt
    type(site_massbal_type), pointer   :: site_mass
    type(elem_diag_type), pointer      :: elflux_diags

    real(r8) :: donatable_mass       ! non-burned litter mass provided by the donor [kg]
                                     ! some may or may not be retained by the donor
    real(r8) :: burned_mass          ! the mass of litter that was supposed to be provided
                                     ! by the donor, but was burned [kg]
    real(r8) :: remainder_area       ! current patch's remaining area after donation [m2]
    real(r8) :: retain_frac          ! the fraction of litter mass retained by the donor patch
    real(r8) :: bcroot               ! amount of below ground coarse root per cohort kg
    real(r8) :: bstem                ! amount of above ground stem biomass per cohort kg
    real(r8) :: leaf_burn_frac       ! fraction of leaves burned
    real(r8) :: leaf_m               ! leaf mass [kg]
    real(r8) :: fnrt_m               ! fineroot mass [kg]
    real(r8) :: sapw_m               ! sapwood mass [kg]
    real(r8) :: store_m              ! storage mass [kg]
    real(r8) :: struct_m             ! structure mass [kg]
    real(r8) :: repro_m              ! Reproductive mass (seeds/flowers) [kg]
    real(r8) :: num_dead_trees       ! total number of dead trees passed in with the burn area
    real(r8) :: num_live_trees       ! total number of live trees passed in with the burn area
    real(r8) :: donate_m2            ! area normalization for litter mass destined to new patch [m-2]
    real(r8) :: retain_m2            ! area normalization for litter mass staying in donor patch [m-2]
    real(r8) :: dcmpy_frac           ! fraction of mass going to each decomposability partition
    integer  :: el                   ! element loop index
    integer  :: sl                   ! soil layer index
    integer  :: c                    ! loop index for coarse woody debris pools
    integer  :: pft                  ! loop index for plant functional types
    integer  :: dcmpy                ! loop index for decomposability pool
    integer  :: element_id           ! parteh compatible global element index
    real(r8) :: trunk_product_site   ! flux of carbon in trunk products exported off site      [ kgC/site ]
                                     ! (note we are accumulating over the patch, but scale is site level)
    real(r8) :: woodproduct_mass     ! mass that ends up in wood products [kg]

    !---------------------------------------------------------------------

    clear_veg_if: if (clearing_matrix_element) then

       ! If plant hydraulics are turned on, account for water leaving the plant-soil
       ! mass balance through the dead trees
       if (hlm_use_planthydro == itrue) then
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             num_dead_trees  = (currentCohort%n*patch_site_areadis/currentPatch%area)
             call AccumulateMortalityWaterStorage(currentSite,currentCohort,num_dead_trees)
             currentCohort => currentCohort%taller
          end do
       end if


       ! If/when sending litter fluxes to the donor patch, we divide the total 
       ! mass sent to that patch, by the area it will have remaining
       ! after it donates area.
       ! i.e. subtract the area it is donating.

       remainder_area = currentPatch%area - patch_site_areadis

       ! Calculate the fraction of litter to be retained versus donated
       ! vis-a-vis the new and donor patch (if the area remaining
       ! in the original donor patch is small, don't bother
       ! retaining anything.)
       retain_frac = (1.0_r8-landusechange_localization) * &
            remainder_area/(newPatch%area+remainder_area)

       if(remainder_area > rsnbl_math_prec) then
          retain_m2 = retain_frac/remainder_area
          donate_m2 = (1.0_r8-retain_frac)/newPatch%area
       else
          retain_m2 = 0._r8
          donate_m2 = 1.0_r8/newPatch%area
       end if

       do el = 1,num_elements

          ! Zero some site level accumulator diagnsotics
          trunk_product_site  = 0.0_r8

          element_id = element_list(el)
          site_mass  => currentSite%mass_balance(el)
          elflux_diags => currentSite%flux_diags%elem(el)
          curr_litt  => currentPatch%litter(el)      ! Litter pool of "current" patch
          new_litt   => newPatch%litter(el)          ! Litter pool of "new" patch

          ! -----------------------------------------------------------------------------
          ! PART 1) Handle mass fluxes associated with plants that died in the land use transition
          ! ------------------------------------------------------------------------------

          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))

             pft = currentCohort%pft

             ! Number of trees that died because of the land use transition, per m2 of ground.
             ! Divide their litter into the four litter streams, and spread
             ! across ground surface.
             ! -----------------------------------------------------------------------

             fnrt_m   = currentCohort%prt%GetState(fnrt_organ, element_id)
             store_m  = currentCohort%prt%GetState(store_organ, element_id)
             repro_m  = currentCohort%prt%GetState(repro_organ, element_id)

             if (prt_params%woody(currentCohort%pft) == itrue) then
                ! Assumption: for woody plants fluxes from deadwood and sapwood go together in CWD pool
                leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id)
                sapw_m          = currentCohort%prt%GetState(sapw_organ,element_id)
                struct_m        = currentCohort%prt%GetState(struct_organ,element_id)
             else
                ! for non-woody plants all stem fluxes go into the same leaf litter pool
                leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id) + &
                     currentCohort%prt%GetState(sapw_organ,element_id) + &
                     currentCohort%prt%GetState(struct_organ,element_id)
                sapw_m          = 0._r8
                struct_m        = 0._r8
             end if


             ! Absolute number of dead trees being transfered in with the donated area
             num_dead_trees = (currentCohort%n * &
                  patch_site_areadis/currentPatch%area)

             ! Contribution of dead trees to leaf litter
             donatable_mass = num_dead_trees * (leaf_m+repro_m) * &
                  (1.0_r8-EDPftvarcon_inst%landusechange_frac_burned(pft))

             ! Contribution of dead trees to leaf burn-flux
             burned_mass  = num_dead_trees * (leaf_m+repro_m) * EDPftvarcon_inst%landusechange_frac_burned(pft)

             do dcmpy=1,ndcmpy
                dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)
                new_litt%leaf_fines(dcmpy) = new_litt%leaf_fines(dcmpy) + &
                     donatable_mass*donate_m2*dcmpy_frac
                curr_litt%leaf_fines(dcmpy) = curr_litt%leaf_fines(dcmpy) + &
                     donatable_mass*retain_m2*dcmpy_frac
             end do

             site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass
             
             call set_root_fraction(currentSite%rootfrac_scr, pft, currentSite%zi_soil, &
                  bc_in%max_rooting_depth_index_col)

             ! Contribution of dead trees to root litter (no root burn flux to atm)
             do dcmpy=1,ndcmpy
                dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
                do sl = 1,currentSite%nlevsoil
                   donatable_mass = num_dead_trees * (fnrt_m+store_m) * currentSite%rootfrac_scr(sl)
                   new_litt%root_fines(dcmpy,sl) = new_litt%root_fines(dcmpy,sl) + &
                        donatable_mass*donate_m2*dcmpy_frac
                   curr_litt%root_fines(dcmpy,sl) = curr_litt%root_fines(dcmpy,sl) + &
                        donatable_mass*retain_m2*dcmpy_frac
                end do
             end do

             ! Track as diagnostic fluxes
             elflux_diags%surf_fine_litter_input(pft) = &
                  elflux_diags%surf_fine_litter_input(pft) + &
                  num_dead_trees * (leaf_m+repro_m) * (1.0_r8-EDPftvarcon_inst%landusechange_frac_burned(pft))

             elflux_diags%root_litter_input(pft) = &
                  elflux_diags%root_litter_input(pft) + &
                  (fnrt_m + store_m) * num_dead_trees

             ! coarse root biomass per tree
             bcroot = (sapw_m + struct_m) * (1.0_r8 - prt_params%allom_agb_frac(pft) )

             ! below ground coarse woody debris from felled trees
             do c = 1,ncwd
                do sl = 1,currentSite%nlevsoil
                   donatable_mass =  num_dead_trees * SF_val_CWD_frac(c) * &
                        bcroot * currentSite%rootfrac_scr(sl)

                   new_litt%bg_cwd(c,sl) = new_litt%bg_cwd(c,sl) + &
                        donatable_mass * donate_m2
                   curr_litt%bg_cwd(c,sl) = curr_litt%bg_cwd(c,sl) + &
                        donatable_mass * retain_m2

                   ! track diagnostics
                   elflux_diags%cwd_bg_input(c) = &
                        elflux_diags%cwd_bg_input(c) + &
                        donatable_mass
                enddo
             end do

             ! stem biomass per tree
             bstem  = (sapw_m + struct_m) * prt_params%allom_agb_frac(pft)

             ! Above ground coarse woody debris from twigs and small branches
             ! a portion of this pool may burn
             ! a portion may also be carried offsite as wood product
             do c = 1,ncwd
                donatable_mass = num_dead_trees * SF_val_CWD_frac(c) * bstem
                if (c == 1 .or. c == 2) then  ! these pools can burn
                   donatable_mass = donatable_mass * (1.0_r8-EDPftvarcon_inst%landusechange_frac_burned(pft))
                   burned_mass = num_dead_trees * SF_val_CWD_frac(c) * bstem * &
                        EDPftvarcon_inst%landusechange_frac_burned(pft)

                   site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass
                else ! all other pools can end up as timber products or burn or go to litter
                   donatable_mass = donatable_mass * (1.0_r8-EDPftvarcon_inst%landusechange_frac_exported(pft)) * &
                        (1.0_r8-EDPftvarcon_inst%landusechange_frac_burned(pft))

                   burned_mass = num_dead_trees * SF_val_CWD_frac(c) * bstem * &
                        (1.0_r8-EDPftvarcon_inst%landusechange_frac_exported(pft)) * &
                        EDPftvarcon_inst%landusechange_frac_burned(pft)

                   woodproduct_mass = num_dead_trees * SF_val_CWD_frac(c) * bstem * &
                        EDPftvarcon_inst%landusechange_frac_exported(pft)

                   site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass

                   trunk_product_site = trunk_product_site + &
                        woodproduct_mass

                   ! Amount of trunk mass exported off site [kg/m2]
                   elflux_diags%exported_harvest = elflux_diags%exported_harvest + &
                        woodproduct_mass * area_inv

                   site_mass%wood_product_landusechange(pft) = site_mass%wood_product_landusechange(pft) + &
                        woodproduct_mass
                   
                endif
                new_litt%ag_cwd(c) = new_litt%ag_cwd(c) + donatable_mass * donate_m2
                curr_litt%ag_cwd(c) = curr_litt%ag_cwd(c) + donatable_mass * retain_m2
                elflux_diags%cwd_ag_input(c) = elflux_diags%cwd_ag_input(c) + donatable_mass
             enddo

             currentCohort => currentCohort%taller
          enddo

          ! Update the amount of carbon exported from the site through logging.

          if(element_id .eq. carbon12_element) then
             currentSite%resources_management%trunk_product_site  = &
                  currentSite%resources_management%trunk_product_site + &
                  trunk_product_site
          end if


       end do

    end if clear_veg_if
    return
  end subroutine landusechange_litter_fluxes

  ! ============================================================================

  subroutine fuse_patches( csite, bc_in )
    !
    ! !DESCRIPTION:
    !  Decide to fuse patches if their cohort structures are similar           
    !
    ! !USES:
    use EDParamsMod , only : ED_val_patch_fusion_tol
    use EDTypesMod , only : patch_fusion_tolerance_relaxation_increment
    use EDTypesMod , only : max_age_of_second_oldest_patch
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target  :: csite
    type(bc_in_type), intent(in)               :: bc_in
    !
    ! !LOCAL VARIABLES:
    type(ed_site_type) , pointer :: currentSite
    type(fates_patch_type), pointer :: currentPatch,tpp,tmpptr
    integer  :: ft,z        !counters for pft and height class
    real(r8) :: norm        !normalized difference between biomass profiles
    real(r8) :: profiletol  !tolerance of patch fusion routine. Starts off high and is reduced if there are too many patches.
    integer  :: nopatches(n_landuse_cats)   !number of patches presently in gridcell
    integer  :: iterate     !switch of patch reduction iteration scheme. 1 to keep going, 0 to stop
    integer  :: fuse_flag   !do patches get fused (1) or not (0).
    integer  :: i_lulabel    !iterator over anthropogenic disturbance categories
    integer  :: i_pftlabel  !nocomp pft iterator
    real(r8) :: primary_land_fraction_beforefusion,primary_land_fraction_afterfusion
    integer  :: pftlabelmin, pftlabelmax
    integer  :: num_bareground_patches
    integer  :: i
    !
    !---------------------------------------------------------------------

    currentSite => csite 

    profiletol = ED_val_patch_fusion_tol

    primary_land_fraction_beforefusion = 0._r8
    primary_land_fraction_afterfusion = 0._r8

    nopatches(1:n_landuse_cats) = 0
    num_bareground_patches = 0
    
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       if ( currentPatch%land_use_label .gt. nocomp_bareground_land) then
          nopatches(currentPatch%land_use_label) = &
               nopatches(currentPatch%land_use_label) + 1
       
          if (currentPatch%land_use_label .eq. primaryland) then
             primary_land_fraction_beforefusion = primary_land_fraction_beforefusion + &
                  currentPatch%area * AREA_INV
          endif
       else
          num_bareground_patches = num_bareground_patches + 1
       endif
       currentPatch => currentPatch%older
    enddo

    if (num_bareground_patches .gt. 1 ) then
       write(fates_log(),*) 'somehow there is more than one bare ground patch. this shouldnt have happened.'
       call endrun(msg=errMsg(sourcefile, __LINE__))                
    endif
    
    pftlabelmin = 0
    if ( hlm_use_nocomp .eq. itrue ) then
       pftlabelmax = numpft
    else
       pftlabelmax = 0
    endif

    !---------------------------------------------------------------------!
    ! iterate over land use categories
    !---------------------------------------------------------------------!    

    lulabel_loop: do i_lulabel = 1, n_landuse_cats

       !---------------------------------------------------------------------!
       !  We only really care about fusing patches if nopatches > 1          !
       !---------------------------------------------------------------------!

       iterate = 1

       !---------------------------------------------------------------------!
       !  Keep doing this until nopatches <= maxpatches_by_landuse(i_lulabel)!
       !---------------------------------------------------------------------!

       iterate_eq_1_loop: do while(iterate == 1)

        !---------------------------------------------------------------------!
        ! iterate over nocomp pft labels (if nocomp is false, then this isn't much of a loop)
        !---------------------------------------------------------------------!

        pftlabel_loop: do i_pftlabel = pftlabelmin, pftlabelmax

          !---------------------------------------------------------------------!
          ! Calculate the biomass profile of each patch                         !
          !---------------------------------------------------------------------!  
          currentPatch => currentSite%youngest_patch
          do while(associated(currentPatch))
             call patch_pft_size_profile(currentPatch)
             currentPatch => currentPatch%older
          enddo

          !-------------------------------------------------------------------------------!
          ! Loop round current & target (currentPatch,tpp) patches to assess combinations !
          !-------------------------------------------------------------------------------!   
          currentPatch => currentSite%youngest_patch
          currentpatch_loop: do while(associated(currentPatch))      
             tpp => currentSite%youngest_patch
             tpp_loop: do while(associated(tpp))

                both_associated_if: if(associated(tpp).and.associated(currentPatch))then
                   !--------------------------------------------------------------------!
                   ! only fuse patches whose anthropogenic disturbance category matches !
                   ! that of the outer loop that we are in                              !
                   !--------------------------------------------------------------------!
                   landuse_labels_match_if: if ( tpp%land_use_label .eq. i_lulabel .and. &
                        currentPatch%land_use_label .eq. i_lulabel) then

                    nocomp_pft_labels_match_if: if (hlm_use_nocomp .eq. ifalse .or. &
                         (tpp%nocomp_pft_label .eq. i_pftlabel .and. &
                         currentPatch%nocomp_pft_label .eq. i_pftlabel)) then

                      !--------------------------------------------------------------------------------------------
                      ! The default is to fuse the patches, unless some criteria is met which keeps them separated.
                      ! there are multiple criteria which all need to be met to keep them distinct:
                      ! (a) one of them is younger than the max age at which we force fusion;
                      ! (b) there is more than a threshold (tiny) amount of biomass in at least one of the patches;
                      ! (c) for at least one pft x size class, where there is biomass in that class in at least one patch,
                      ! and the normalized difference between the patches exceeds a threshold.
                      !--------------------------------------------------------------------------------------------

                      fuse_flag = 1
                      different_patches_if: if(currentPatch%patchno /= tpp%patchno) then   !these should be the same patch

                         !-----------------------------------------------------------------------------------
                         ! check to see if both patches are older than the age at which we force them to fuse
                         !-----------------------------------------------------------------------------------

                         maxage_if: if ( tpp%age .le. max_age_of_second_oldest_patch .or. &
                              currentPatch%age .le. max_age_of_second_oldest_patch ) then


                            !------------------------------------------------------------
                            ! the next bit of logic forces fusion of two patches which 
                            ! both have tiny biomass densities. without this,
                            ! fates gives a bunch of really young patches which all have 
                            ! almost no biomass and so don't need to be distinguished 
                            ! from each other. but if force_patchfuse_min_biomass is too big,
                            ! it takes too long for the youngest patch to build up enough 
                            ! biomass to be its own distinct entity, which leads to large 
                            ! oscillations in the patch dynamics and dependent variables.
                            !------------------------------------------------------------

                            patchfuse_min_biomass_if: if &
                                 (sum(currentPatch%pft_agb_profile(:,:)) > force_patchfuse_min_biomass .or. &
                                 sum(tpp%pft_agb_profile(:,:)) > force_patchfuse_min_biomass ) then

                               !---------------------------------------------------------------------!
                               ! Calculate the difference criteria for each pft and dbh class        !
                               !---------------------------------------------------------------------!   

                               pft_loop: do ft = 1,numpft        ! loop over pfts
                                  hgt_bin_loop: do z = 1,n_dbh_bins      ! loop over hgt bins 

                                     !----------------------------------
                                     ! is there biomass in this category?
                                     !----------------------------------

                                     agbprof_gt_zero_if: if &
                                          (currentPatch%pft_agb_profile(ft,z)  > 0.0_r8 .or.  &
                                          tpp%pft_agb_profile(ft,z) > 0.0_r8)then 

                                        !---------------------------------------------------------------------!
                                        ! what is the relative difference in biomass in this category between
                                        ! the two patches?
                                        !---------------------------------------------------------------------!

                                        norm = abs(currentPatch%pft_agb_profile(ft,z) - &
                                             tpp%pft_agb_profile(ft,z))/(0.5_r8 * &
                                             &(currentPatch%pft_agb_profile(ft,z) + tpp%pft_agb_profile(ft,z)))

                                        !---------------------------------------------------------------------!
                                        ! Look for differences in profile biomass, above the minimum biomass  !
                                        !---------------------------------------------------------------------!

                                        if(norm  > profiletol)then

                                           fuse_flag = 0 !do not fuse  - keep apart. 

                                        endif
                                     endif agbprof_gt_zero_if
                                  enddo hgt_bin_loop
                               enddo pft_loop
                            endif patchfuse_min_biomass_if
                         endif maxage_if

                         !-------------------------------------------------------------------------!
                         ! Call the patch fusion routine if there is not a meaningful difference   !
                         ! any of the pft x height categories                                      !
                         ! or both are older than forced fusion age                                !
                         !-------------------------------------------------------------------------!

                         fuseflagset_if: if(fuse_flag  ==  1)then
                            
                            !-----------------------!
                            ! fuse the two patches  !
                            !-----------------------!
                            
                            tmpptr => currentPatch%older       
                            call fuse_2_patches(csite, currentPatch, tpp)
                            call fuse_cohorts(csite,tpp, bc_in)
                            call tpp%SortCohorts()
                            call tpp%ValidateCohorts()
                            currentPatch => tmpptr

                            !------------------------------------------------------------------------!
                            ! since we've just fused two patches, but are still in the midst of      !
                            ! a patch x patch loop, reset the patch fusion tolerance to the starting !
                            ! value so that any subsequent fusions in this loop are done with that   !
                            ! value. otherwise we can end up in a situation where we've loosened the !
                            ! fusion tolerance to get nopatches <= maxpatches_by_landuse(i_lulabel), but then,      !
                            ! having accomplished that, we continue through all the patch x patch    !
                            ! combinations and then all the patches get fused, ending up with        !
                            ! nopatches << maxpatches_by_landuse(i_lulabel) and losing all heterogeneity.           !
                            !------------------------------------------------------------------------!

                            profiletol = ED_val_patch_fusion_tol

                         endif fuseflagset_if
                      endif different_patches_if
                    endif nocomp_pft_labels_match_if
                   endif landuse_labels_match_if
                endif both_associated_if

                tpp => tpp%older
             enddo tpp_loop

             if(associated(currentPatch))then 
                currentPatch => currentPatch%older 
             else
                currentPatch => null()
             endif !associated currentPatch

          enddo currentpatch_loop

        end do pftlabel_loop

          !---------------------------------------------------------------------!
          ! Is the number of patches larger than the maximum?                   !
          !---------------------------------------------------------------------!   
          nopatches(i_lulabel) = 0
          currentPatch => currentSite%youngest_patch
          do while(associated(currentPatch))
             if (currentPatch%land_use_label .eq. i_lulabel) then
                nopatches(i_lulabel) = nopatches(i_lulabel) +1
             endif
             currentPatch => currentPatch%older
          enddo

          if(nopatches(i_lulabel) > maxpatches_by_landuse(i_lulabel))then
             iterate = 1
             profiletol = profiletol * patch_fusion_tolerance_relaxation_increment

             !---------------------------------------------------------------------!
             ! Making profile tolerance larger means that more fusion will happen  !
             !---------------------------------------------------------------------!        

             ! its possible that there are too many categorical patch types and the tolerances
             ! will never allow patch fusion to occur.  In this case crash and let the user know.
             ! the 100 is sort of a random number, in principle since profile tolerance is compared 
             ! against relative biomass size, it shoudnt ever get above 2 (which would mean fusing 
             ! a zero with a nonzero biomass in a given category)
             if (profiletol .gt. 100._r8) then
                write(fates_log(),*) 'profile tolerance is too big, this shouldnt happen.'
                write(fates_log(),*) 'probably this means there are too many distinct categorical '
                write(fates_log(),*) 'patch types for the maximum number of patches'
                call dump_site(currentSite)
                write(fates_log(),*) 'currentSite%area_bareground', currentSite%area_bareground
                do i = 1, n_landuse_cats
                   write(fates_log(),*) 'i, currentSite%area_pft(:,i)',i, currentSite%area_pft(:,i)
                end do
                tmpptr => currentSite%youngest_patch
                do while(associated(tmpptr))
                   write(fates_log(),*) tmpptr%area, tmpptr%nocomp_pft_label, tmpptr%land_use_label
                   tmpptr => tmpptr%older
                end do
                
                call endrun(msg=errMsg(sourcefile, __LINE__))                
             endif
          else
             iterate = 0
          endif

       enddo iterate_eq_1_loop ! iterate .eq. 1 ==> nopatches>maxpatches_by_landuse(i_lulabel)

    end do lulabel_loop

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))

       if (currentPatch%land_use_label .eq. primaryland) then
          primary_land_fraction_afterfusion = primary_land_fraction_afterfusion + &
               currentPatch%area * AREA_INV
       endif

       currentPatch => currentPatch%older
    enddo

    currentSite%primary_land_patchfusion_error = primary_land_fraction_afterfusion - primary_land_fraction_beforefusion
 
  end subroutine fuse_patches

  ! ============================================================================

  subroutine fuse_2_patches(csite, dp, rp)
    !
    ! !DESCRIPTION:
    ! This function fuses the two patches specified in the argument.
    ! It fuses the first patch in the argument (the "donor") into the second
    ! patch in the argument (the "recipient"), and frees the memory 
    ! associated with the secnd patch
    !
    ! !USES:
    use FatesSizeAgeTypeIndicesMod, only: get_age_class_index
    !
    ! !ARGUMENTS:
    type (ed_site_type), intent(inout),target :: csite  ! Current site 
    type (fates_patch_type) , pointer :: dp                ! Donor Patch
    type (fates_patch_type) , target, intent(inout) :: rp  ! Recipient Patch

    !
    ! !LOCAL VARIABLES:
    type (fates_cohort_type), pointer :: currentCohort ! Current Cohort
    type (fates_cohort_type), pointer :: nextc         ! Remembers next cohort in list 
    integer                        :: c,p          !counters for pft and litter size class. 
    integer                        :: el           ! loop counting index for elements
    integer                        :: pft          ! loop counter for pfts
    type(fates_patch_type), pointer   :: youngerp     ! pointer to the patch younger than donor
    type(fates_patch_type), pointer   :: olderp       ! pointer to the patch older than donor
    real(r8)                       :: inv_sum_area ! Inverse of the sum of the two patches areas
    !-----------------------------------------------------------------------------------------------

    ! Generate a litany of area weighted averages

    inv_sum_area = 1.0_r8/(dp%area + rp%area)
    
    rp%age = (dp%age * dp%area + rp%age * rp%area) * inv_sum_area
    rp%age_since_anthro_disturbance = (dp%age_since_anthro_disturbance * dp%area &
         + rp%age_since_anthro_disturbance * rp%area) * inv_sum_area

    rp%age_class = get_age_class_index(rp%age)

    do el = 1,num_elements
       call rp%litter(el)%FuseLitter(rp%area,dp%area,dp%litter(el))
    end do
    
    call rp%fuel%Fuse(rp%area, dp%area, dp%fuel)

    if ( rp%land_use_label .ne. dp%land_use_label) then
       write(fates_log(),*) 'trying to fuse patches with different land_use_label values'
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    if ( hlm_use_nocomp .eq. itrue .and. rp%nocomp_pft_label .ne. dp%nocomp_pft_label) then
       write(fates_log(),*) 'trying to fuse patches with different nocomp_pft_label values'
       write(fates_log(),*) 'rp%nocomp_pft_label, dp%nocomp_pft_label',rp%nocomp_pft_label, dp%nocomp_pft_label
       write(fates_log(),*) 'rp%area, dp%area',rp%area, dp%area
       write(fates_log(),*) 'sum(rp%pft_agb_profile(:,:), sum(dp%pft_agb_profile(:,:)',sum(rp%pft_agb_profile(:,:)), sum(dp%pft_agb_profile(:,:))
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! Weighted mean of the running means
    call rp%tveg24%FuseRMean(dp%tveg24,rp%area*inv_sum_area)
    call rp%tveg_lpa%FuseRMean(dp%tveg_lpa,rp%area*inv_sum_area)

    if ( regeneration_model == TRS_regeneration ) then
       call rp%seedling_layer_par24%FuseRMean(dp%seedling_layer_par24,rp%area*inv_sum_area)
       call rp%sdlng_mort_par%FuseRMean(dp%sdlng_mort_par,rp%area*inv_sum_area)
       call rp%sdlng2sap_par%FuseRMean(dp%sdlng2sap_par,rp%area*inv_sum_area)
       do pft = 1,numpft
          call rp%sdlng_emerg_smp(pft)%p%FuseRMean(dp%sdlng_emerg_smp(pft)%p,rp%area*inv_sum_area)
          call rp%sdlng_mdd(pft)%p%FuseRMean(dp%sdlng_mdd(pft)%p,rp%area*inv_sum_area)
       enddo
    end if
    
    call rp%tveg_longterm%FuseRMean(dp%tveg_longterm,rp%area*inv_sum_area)

    rp%livegrass               = (dp%livegrass*dp%area + rp%livegrass*rp%area) * inv_sum_area
    rp%ros_front               = (dp%ros_front*dp%area + rp%ros_front*rp%area) * inv_sum_area
    rp%tau_l                   = (dp%tau_l*dp%area + rp%tau_l*rp%area) * inv_sum_area
    rp%tfc_ros              = (dp%tfc_ros*dp%area + rp%tfc_ros*rp%area) * inv_sum_area
    rp%fi                   = (dp%fi*dp%area + rp%fi*rp%area) * inv_sum_area
    rp%fd                   = (dp%fd*dp%area + rp%fd*rp%area) * inv_sum_area
    rp%ros_back             = (dp%ros_back*dp%area + rp%ros_back*rp%area) * inv_sum_area
    rp%scorch_ht(:)         = (dp%scorch_ht(:)*dp%area + rp%scorch_ht(:)*rp%area) * inv_sum_area
    rp%frac_burnt           = (dp%frac_burnt*dp%area + rp%frac_burnt*rp%area) * inv_sum_area
    rp%btran_ft(:)          = (dp%btran_ft(:)*dp%area + rp%btran_ft(:)*rp%area) * inv_sum_area
    rp%zstar                = (dp%zstar*dp%area + rp%zstar*rp%area) * inv_sum_area
    rp%c_stomata            = (dp%c_stomata*dp%area + rp%c_stomata*rp%area) * inv_sum_area
    rp%c_lblayer            = (dp%c_lblayer*dp%area + rp%c_lblayer*rp%area) * inv_sum_area
    rp%rad_error(1)         = (dp%rad_error(1)*dp%area + rp%rad_error(1)*rp%area) * inv_sum_area
    rp%rad_error(2)         = (dp%rad_error(2)*dp%area + rp%rad_error(2)*rp%area) * inv_sum_area
    
    rp%area = rp%area + dp%area !THIS MUST COME AT THE END!

    !insert donor cohorts into recipient patch
    if(associated(dp%shortest))then

       currentCohort => dp%shortest
       if(associated(currentCohort)) then
          nextc => currentCohort%taller
       endif

       do while(associated(dp%shortest))
       
         call rp%InsertCohort(currentCohort)

          currentCohort => nextc

          dp%shortest => currentCohort

          if(associated(currentCohort)) then
             nextc => currentCohort%taller
          endif

       enddo !cohort
       call rp%ValidateCohorts()
    endif !are there any cohorts?

    call patch_pft_size_profile(rp) ! Recalculate the patch size profile for the resulting patch

    ! Define some aliases for the donor patches younger and older neighbors
    ! which may or may not exist.  After we set them, we will remove the donor
    ! And then we will go about re-setting the map.
    if(associated(dp%older))then
       olderp => dp%older
    else
       olderp => null()
    end if
    if(associated(dp%younger))then
       youngerp => dp%younger
    else
       youngerp => null()
    end if

    ! We have no need for the dp pointer anymore, we have passed on it's legacy
    call dp%FreeMemory(regeneration_model, numpft)
    deallocate(dp, stat=istat, errmsg=smsg)
    if (istat/=0) then
       write(fates_log(),*) 'dealloc006: fail on deallocate(dp):'//trim(smsg)
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    ! if neither youngerp nor olderp are associated, that means that the patch we are fusing into
    ! is not part of the linked-list structure, and so no further action needs to be taken.
    if(associated(youngerp) .or. associated(olderp))then

       if(associated(youngerp))then
          ! Update the younger patch's new older patch (because it isn't dp anymore)
          youngerp%older => olderp
       else
          ! There was no younger patch than dp, so the head of the young order needs
          ! to be set, and it is set as the patch older than dp.  That patch
          ! already knows it's older patch (so no need to set or change it)
          csite%youngest_patch => olderp
          olderp%younger       => null()
       end if

       if(associated(olderp))then
          ! Update the older patch's new younger patch (becuase it isn't dp anymore)
          olderp%younger => youngerp
       else
          ! There was no patch older than dp, so the head of the old patch order needs
          ! to be set, and it is set as the patch younger than dp.  That patch already
          ! knows it's younger patch, no need to set
          csite%oldest_patch => youngerp
          youngerp%older     => null()
       end if
    end if

  end subroutine fuse_2_patches

  ! ============================================================================

  subroutine terminate_patches(currentSite, bc_in)
    !
    ! !DESCRIPTION:
    !  Terminate Patches if they  are too small                          
    !
    !
    ! !ARGUMENTS:
    type(ed_site_type), target, intent(inout) :: currentSite
    type(bc_in_type), intent(in)               :: bc_in
    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type), pointer :: currentPatch
    type(fates_patch_type), pointer :: olderPatch
    type(fates_patch_type), pointer :: youngerPatch
    type(fates_patch_type), pointer :: patchpointer
    type(fates_patch_type), pointer :: largest_patch
    integer, parameter           :: max_cycles = 10  ! After 10 loops through
                                                     ! You should had fused
    integer                      :: count_cycles
    logical                      :: gotfused
    logical                      :: current_patch_is_youngest_lutype
    integer                      :: i_landuse, i_pft
    integer                      :: land_use_type_to_remove

    real(r8) areatot ! variable for checking whether the total patch area is wrong.
    real(r8) :: state_vector_driver(n_landuse_cats)  ! [m2/m2]
    real(r8) :: state_vector_internal(n_landuse_cats)  ! [m2/m2]
    !---------------------------------------------------------------------

    ! Initialize the count cycles
    count_cycles = 0

    ! Start at the youngest patch in the list and assume that the largest patch is this patch
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       lessthan_min_patcharea_if: if(currentPatch%area <= min_patch_area)then

          nocomp_if: if (hlm_use_nocomp .eq. itrue) then

             gotfused = .false.
             patchpointer => currentSite%youngest_patch
             do while(associated(patchpointer))
                if ( .not.associated(currentPatch,patchpointer) .and. &
                     patchpointer%nocomp_pft_label .eq. currentPatch%nocomp_pft_label .and. &
                     patchpointer%land_use_label .eq. currentPatch%land_use_label .and. &
                     .not. gotfused) then

                   call fuse_2_patches(currentSite, patchpointer, currentPatch)
                   
                   gotfused = .true.
                else
                   patchpointer => patchpointer%older
                endif
             end do

             if ( .not. gotfused ) then
                !! somehow didn't find a patch to fuse with.
                write(fates_log(),*) 'Warning. small nocomp patch wasnt able to find another patch to fuse with.', &
                     currentPatch%nocomp_pft_label, currentPatch%land_use_label, currentPatch%area
             endif

          else nocomp_if

             ! Even if the patch area is small, avoid fusing it into its neighbor
             ! if it is the youngest of all patches. We do this in attempts to maintain
             ! a discrete patch for very young patches.
             ! However, if the patch to be fused is excessively small, then fuse
             ! at all costs.

             notyoungest_if: if ( .not.associated(currentPatch,currentSite%youngest_patch) .or. &
               currentPatch%area <= min_patch_area_forced ) then

                gotfused = .false.

                associated_older_if: if(associated(currentPatch%older)) then

                   if(debug) &
                        write(fates_log(),*) 'fusing to older patch because this one is too small',&
                        currentPatch%area, &
                        currentPatch%older%area

                   ! We set a pointer to this patch, because
                   ! it will be returned by the subroutine as de-referenced

                   olderPatch => currentPatch%older

                   ! If the older patch has the same landuse label fuse the older (donor) patch into the current patch
                   distlabel_1_if: if (currentPatch%land_use_label .eq. olderPatch%land_use_label) then

                      if(debug) &
                           write(fates_log(),*) 'terminate: fused to older patch, same label: ', currentPatch%land_use_label, olderPatch%land_use_label

                      call fuse_2_patches(currentSite, olderPatch, currentPatch)

                      ! The fusion process has updated the "older" pointer on currentPatch
                      ! for us.

                      ! This logic checks to make sure that the younger patch is not the youngest
                      ! patch. As mentioned earlier, we try not to fuse it.

                      gotfused = .true.
                   else distlabel_1_if !i.e. anthro labels of two patches are not the same
                      countcycles_if: if (count_cycles .gt. 0) then
                         ! if we're having an incredibly hard time fusing patches because of their differing anthropogenic disturbance labels,
                         ! since the size is so small, let's sweep the problem under the rug and change the tiny patch's label to that of its older sibling
                         ! and then allow them to fuse together.
                         currentPatch%land_use_label = olderPatch%land_use_label
                         currentPatch%age_since_anthro_disturbance = olderPatch%age_since_anthro_disturbance
                         call fuse_2_patches(currentSite, olderPatch, currentPatch)
                         gotfused = .true.
                      endif countcycles_if
                   endif distlabel_1_if
                endif associated_older_if

                not_gotfused_if: if( .not. gotfused .and. associated(currentPatch%younger) ) then

                   if(debug) &
                        write(fates_log(),*) 'fusing to younger patch because oldest one is too small', &
                        currentPatch%area

                   youngerPatch => currentPatch%younger

                   distlabel_2_if: if (currentPatch%land_use_label .eq. youngerPatch%land_use_label) then

                      call fuse_2_patches(currentSite, youngerPatch, currentPatch)

                      ! The fusion process has updated the "younger" pointer on currentPatch

                      gotfused = .true.

                   else distlabel_2_if
                      if (count_cycles .gt. 0) then
                         ! if we're having an incredibly hard time fusing patches because of their differing anthropogenic disturbance labels,
                         ! since the size is so small, let's sweep the problem under the rug and change the tiny patch's label to that of its younger sibling
                         currentPatch%land_use_label = youngerPatch%land_use_label
                         currentPatch%age_since_anthro_disturbance = youngerPatch%age_since_anthro_disturbance
                         call fuse_2_patches(currentSite, youngerPatch, currentPatch)
                         gotfused = .true.
                      endif ! count cycles
                   endif distlabel_2_if     ! anthro labels
                endif not_gotfused_if ! has an older patch
             endif notyoungest_if ! is not the youngest patch
        endif nocomp_if
        endif lessthan_min_patcharea_if ! very small patch

       ! It is possible that an incredibly small patch just fused into another incredibly
       ! small patch, resulting in an incredibly small patch.  It is also possible that this
       ! resulting incredibly small patch is the oldest patch.  If this was true than
       ! we would had been at the end of the loop, and left with an incredibly small patch.
       ! Think this is impossible? No, this really happens, especially when we have fires.
       ! So, we don't move forward until we have merged enough area into this thing.

        if(currentPatch%area > min_patch_area_forced)then
          currentPatch => currentPatch%older
         
          count_cycles = 0
       else
          count_cycles = count_cycles + 1
       end if

       if(count_cycles > max_cycles) then
          write(fates_log(),*) 'FATES is having difficulties fusing very small patches.'
          write(fates_log(),*) 'It is possible that a either a secondary or primary'
          write(fates_log(),*) 'patch has become the only patch of its kind, and it is'
          write(fates_log(),*) 'is very very small. You can test your luck by'
          write(fates_log(),*) 'disabling the endrun statement following this message.'
          write(fates_log(),*) 'FATES may or may not continue to operate within error'
          write(fates_log(),*) 'tolerances, but will generate another fail if it does not.'
          write(fates_log(),*) 'otherwise, dumping some diagnostics.'
          write(fates_log(),*) currentPatch%area, currentPatch%nocomp_pft_label, currentPatch%land_use_label
          call dump_site(currentSite)

          write(fates_log(),*) 'currentSite%area_bareground', currentSite%area_bareground
          write(fates_log(),*) 'currentSite%area_pft(:,:)', currentSite%area_pft(:,:)
          patchpointer => currentSite%youngest_patch
          do while(associated(patchpointer))
             write(fates_log(),*) patchpointer%area, patchpointer%nocomp_pft_label, patchpointer%land_use_label
             patchpointer => patchpointer%older
          end do
          state_vector_internal = currentSite%get_current_landuse_statevector()
          write(fates_log(),*) 'current landuse state vector: ', state_vector_internal
          write(fates_log(),*) 'current landuse state vector (not including bare gruond): ', state_vector_internal/(1._r8-currentSite%area_bareground)
          call GetLUHStatedata(bc_in, state_vector_driver)
          write(fates_log(),*) 'driver data landuse state vector: ', state_vector_driver
          write(fates_log(),*) 'min_allowed_landuse_fraction: ', currentSite%min_allowed_landuse_fraction
          write(fates_log(),*) 'landuse_vector_gt_min: ', currentSite%landuse_vector_gt_min
          do i_landuse = 1, n_landuse_cats
             write(fates_log(),*) 'trans matrix from: ', i_landuse, currentSite%landuse_transition_matrix(i_landuse,:)
          end do

          if ( (state_vector_driver(currentPatch%land_use_label) .lt. currentSite%min_allowed_landuse_fraction ) .or. &
               (state_vector_internal(currentPatch%land_use_label) .lt. currentSite%min_allowed_landuse_fraction ) ) then

             ! try fusing all of the patches with this land use label into the largest patch on the site.
             land_use_type_to_remove = currentPatch%land_use_label

             write(fates_log(),*) 'removing all patches with land use type ',land_use_type_to_remove

             ! first find the largest patch on the site
             patchpointer => currentSite%youngest_patch
             largest_patch => currentSite%youngest_patch
             do while(associated(patchpointer))
                if (patchpointer%area .gt. largest_patch%area .and. patchpointer%nocomp_pft_label .ne. nocomp_bareground) then
                   largest_patch => patchpointer
                endif
                patchpointer => patchpointer%older
             end do

             ! now go and fuse all patches that have the land use type we are removing into that patch
             patchpointer => currentSite%youngest_patch
             do while(associated(patchpointer))
                if ( patchpointer%land_use_label .eq. land_use_type_to_remove ) then

                   write(fates_log(),*) 'fusing into patch with types, age, and size of:', largest_patch%land_use_label, &
                        largest_patch%nocomp_pft_label, largest_patch%age, largest_patch%area

                   write(fates_log(),*) 'fusing away patch with types, age, and size of:', patchpointer%land_use_label, &
                        patchpointer%nocomp_pft_label, patchpointer%age, patchpointer%area

                   ! reset the categorical properties of the patch and fuse it into the largest patch
                   patchpointer%land_use_label = largest_patch%land_use_label
                   patchpointer%nocomp_pft_label = largest_patch%nocomp_pft_label
                   patchpointer%age_since_anthro_disturbance = largest_patch%age_since_anthro_disturbance
                   call fuse_2_patches(currentSite, patchpointer, largest_patch)

                   ! start over in the loop to make sure we are removing every patch with the targeted land use type
                   patchpointer => currentSite%youngest_patch

                else
                   patchpointer => patchpointer%older
                endif
             end do

             write(fates_log(),*) 'resetting currentSite%landuse_vector_gt_min(i) to .false.'
             ! now reset the allowed land use vector element so that we don't make any more such patches unless they exceed the min area
             currentSite%landuse_vector_gt_min(land_use_type_to_remove) = .false.
             count_cycles = 0
             currentPatch => currentSite%youngest_patch
          else
             write(fates_log(),*) 'this isnt because the land use was less than allowed'

             call endrun(msg=errMsg(sourcefile, __LINE__))
          
             ! Note to user. If you DO decide to remove the end-run above this line
             ! Make sure that you keep the pointer below this line, or you will get
             ! an infinite loop.
             currentPatch => currentPatch%older
             count_cycles = 0
          endif
       end if  !count cycles
       
    enddo ! current patch loop
    
    !check area is not exceeded
    call check_patch_area( currentSite )

    call set_patchno( currentSite, .false., 0)
    
    return
  end subroutine terminate_patches

  ! =====================================================================================

  subroutine DistributeSeeds(currentSite,seed_mass,el,pft)
      
      ! !ARGUMENTS:
      type(ed_site_type), target, intent(inout) :: currentSite  !
      real(r8), intent(in)                      :: seed_mass    ! mass of seed input [kg]
      integer, intent(in)                       :: el           ! element index
      integer, intent(in)                       :: pft          ! pft index

      ! !LOCAL VARIABLES:
      type(fates_patch_type), pointer              :: currentPatch
      type(litter_type), pointer                :: litt

      
      currentPatch => currentSite%oldest_patch
      do while(associated(currentPatch)) 
          litt => currentPatch%litter(el)
          
          if(homogenize_seed_pfts) then
              litt%seed(:) = litt%seed(:) + seed_mass/(area_site*real(numpft,r8))
          else
              litt%seed(pft) = litt%seed(pft) + seed_mass/area_site
          end if
          
          currentPatch => currentPatch%younger
      end do
          

      return
  end subroutine DistributeSeeds

  ! ============================================================================
  subroutine patch_pft_size_profile(cp_pnt)
    !
    ! !DESCRIPTION:
    !  Binned patch size profiles generated for patch fusion routine        
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(fates_patch_type), target, intent(inout) :: cp_pnt
    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type) , pointer  :: currentPatch
    type(fates_cohort_type), pointer  :: currentCohort
    real(r8) :: mind(N_DBH_BINS) ! Bottom of DBH bin 
    real(r8) :: maxd(N_DBH_BINS) ! Top of DBH bin
    real(r8) :: delta_dbh   ! Size of DBH bin
    integer  :: p    ! Counter for PFT 
    integer  :: j    ! Counter for DBH bins 
    real(r8), parameter :: gigantictrees = 1.e8_r8
    !---------------------------------------------------------------------

    currentPatch => cp_pnt

    currentPatch%pft_agb_profile(:,:) = 0.0_r8

    do j = 1,N_DBH_BINS   
        if (j == N_DBH_BINS) then
           mind(j) = patchfusion_dbhbin_loweredges(j)
           maxd(j) = gigantictrees
        else 
           mind(j) = patchfusion_dbhbin_loweredges(j)
           maxd(j) = patchfusion_dbhbin_loweredges(j+1)
        endif
    enddo

    currentCohort => currentPatch%shortest
    do while(associated(currentCohort))    
       do j = 1,N_DBH_BINS   
          if((currentCohort%dbh  >  mind(j)) .AND. (currentCohort%dbh  <=  maxd(j)))then

             currentPatch%pft_agb_profile(currentCohort%pft,j) = &
                  currentPatch%pft_agb_profile(currentCohort%pft,j) + &
                  currentCohort%prt%GetState(struct_organ, carbon12_element) * &
                  currentCohort%n/currentPatch%area

          endif
       enddo ! dbh bins

       currentCohort => currentCohort%taller

    enddo !currentCohort 
   
  end subroutine patch_pft_size_profile

  ! =====================================================================================
  function countPatches( nsites, sites ) result ( totNumPatches ) 
    !
    ! !DESCRIPTION:
    !  Loop over all Patches to count how many there are
    !
    ! !USES:
    use EDTypesMod , only : ed_site_type
    !
    ! !ARGUMENTS:
    integer,             intent(in)            :: nsites
    type(ed_site_type) , intent(inout), target :: sites(nsites)
    !
    ! !LOCAL VARIABLES:
    type (fates_patch_type), pointer :: currentPatch
    integer :: totNumPatches  ! total number of patches.  
    integer :: s
    !---------------------------------------------------------------------

    totNumPatches = 0

    do s = 1,nsites
       currentPatch => sites(s)%oldest_patch
       do while(associated(currentPatch))
          totNumPatches = totNumPatches + 1
          currentPatch => currentPatch%younger
       enddo
    enddo

   end function countPatches

  ! =====================================================================================

 subroutine InsertPatch(currentSite, newPatch)

    ! !DESCRIPTION:
    ! Insert patch into linked list
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type (ed_site_type), intent(inout)              :: currentSite
    type (fates_patch_type), intent(inout), pointer :: newPatch

    ! !LOCAL VARIABLES:
    type (fates_patch_type), pointer :: currentPatch
    logical                       :: patch_inserted

    ! The goal here is to have patches ordered in a specific way. That way is to have them
    ! looped as the following, where LU refers to the land use label, PFT refers to the
    ! nocomp PFT label, and Y and O refer to continuous patch ages.
    !
    !       LU1     ----     LU2       ----     LU3         --  etc
    !    /       \         /       \        /        \
    !  PFT1 --- PFT2  |  PFT1 --- PFT2  |  PFT1 --- PFT2    --  etc
    !  / \      /  \     / \      /  \     / \      /  \  
    ! O - Y    O -  Y   O - Y    O  - Y   O - Y    O -  Y   --  etc

    ! I.e. to treat land use as the outermost loop element, then nocomp PFT as next loop element,
    ! and then age as the innermost loop element.  Visualizing the above as a linked list patches:

    ! LU1/PFT1/O <-> LU1/PFT1/Y <-> LU1/PFT2/O <- ... -> LU3/PFT2/O <-> LU3/PFT2/Y

    ! Mapping this setup onto the existing "older/younger" scheme means that lower number
    ! land use and pft labels are considered "older".  Note that this matches the current
    ! initialization scheme in which patches are created and linked in increasing pft
    ! numerical order starting from 1.  This also aligns with the current set_patchno scheme
    ! in which patches are given an indexable number for the API iteration loops.
    
    ! The way to accomplsh this most simply is to define a pseudo-age that includes all of the
    ! above info and sort the patches based on the pseudo-age. i.e. take some number larger
    ! than any patch will ever reach in actual age. Then take the LU, multiply it by the big
    ! number squared, add it to the pft number multiplied by the big number, and add to the age.
    ! And lastly to sort using that instead of the actual age.

    ! If land use is turned off or nocomp is turned off, then this should devolve to the prior
    ! behavior of just age sorting.

    patch_inserted = .false.
    
    if (GetPseudoPatchAge(newPatch) .le. GetPseudoPatchAge(currentSite%youngest_patch)) then

       ! insert new patch at the head of the linked list
       newPatch%older   => currentSite%youngest_patch
       newPatch%younger => null()
       currentSite%youngest_patch%younger => newPatch
       currentSite%youngest_patch => newPatch

       patch_inserted = .true.
    else if (GetPseudoPatchAge(newPatch) .ge. GetPseudoPatchAge(currentSite%oldest_patch)) then

       ! insert new patch at the end of the linked list
       newPatch%younger   => currentSite%oldest_patch
       newPatch%older     => null()
       currentSite%oldest_patch%older => newPatch
       currentSite%oldest_patch => newPatch

       patch_inserted = .true.
    else
       ! new patch has a pseudo-age somewhere within the linked list. find the first patch which
       ! has a pseudo age older than it, and put it ahead of that patch
       currentPatch => currentSite%youngest_patch
       do while (associated(currentPatch) .and. ( .not. patch_inserted) )   
          if (GetPseudoPatchAge(newPatch) .lt. GetPseudoPatchAge(currentPatch)) then
             newPatch%older => currentPatch
             newPatch%younger => currentPatch%younger
             currentPatch%younger%older => newPatch
             currentPatch%younger => newPatch

             patch_inserted = .true.
          endif
          currentPatch => currentPatch%older
       end do
    end if

    if ( .not. patch_inserted) then
       ! something has gone wrong. abort.
       write(fates_log(),*) 'something has gone wrong in the patch insertion, because no place to put the new patch was found' 
          call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

  end subroutine InsertPatch

  ! =====================================================================================

  function GetPseudoPatchAge(CurrentPatch) result(pseudo_age)
    
    ! Purpose: we want to sort the patches in a way that takes into account both their
    ! continuous and categorical variables. Calculate a pseudo age that does this, by taking
    ! the integer labels, multiplying these by large numbers, and adding to the continuous age.
    ! Note to ensure that lower integer land use label and pft label numbers are considered
    ! "younger" (i.e higher index patchno) in the linked list, they are summed and multiplied by
    ! negative one.  The patch age is still added normally to this negative pseudoage calculation
    ! as a higher age will result in a less negative number correlating with an "older" patch.

    type (fates_patch_type), intent(in), pointer :: CurrentPatch
    real(r8)            :: pseudo_age    
    real(r8), parameter :: max_actual_age = 1.e4  ! hard to imagine a patch older than 10,000 years
    real(r8), parameter :: max_actual_age_squared = 1.e8

    pseudo_age = -1.0_r8 * (real(CurrentPatch%land_use_label,r8) * max_actual_age_squared + &
         real(CurrentPatch%nocomp_pft_label,r8) * max_actual_age) + CurrentPatch%age
    
  end function GetPseudoPatchAge

  ! =====================================================================================

 subroutine CopyPatchMeansTimers(dp, rp)

    ! !DESCRIPTION:
    ! Copy any means or timers from the original patch to the new patch
    ! These values will inherit all info from the original patch
    ! --------------------------------------------------------------------------
    !
    ! !ARGUMENTS:
    type (fates_patch_type), intent(in)    :: dp  ! Donor Patch
    type (fates_patch_type), intent(inout) :: rp  ! Recipient Patch

    ! LOCAL:
    integer :: ipft   ! pft index

    call rp%tveg24%CopyFromDonor(dp%tveg24)
    call rp%tveg_lpa%CopyFromDonor(dp%tveg_lpa)
    call rp%tveg_longterm%CopyFromDonor(dp%tveg_longterm)

    if ( regeneration_model == TRS_regeneration ) then
       call rp%seedling_layer_par24%CopyFromDonor(dp%seedling_layer_par24)
       call rp%sdlng_mort_par%CopyFromDonor(dp%sdlng_mort_par)
       call rp%sdlng2sap_par%CopyFromDonor(dp%sdlng2sap_par)
       do ipft = 1,numpft
          call rp%sdlng_emerg_smp(ipft)%p%CopyFromDonor(dp%sdlng_emerg_smp(ipft)%p)
          call rp%sdlng_mdd(ipft)%p%CopyFromDonor(dp%sdlng_mdd(ipft)%p)
       enddo
    end if

 end subroutine CopyPatchMeansTimers

 end module EDPatchDynamicsMod
