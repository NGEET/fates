module EDPatchDynamicsMod

  ! ============================================================================
  ! Controls formation, creation, fusing and termination of patch level processes. 
  ! ============================================================================
  use FatesGlobals         , only : fates_log 
  use FatesInterfaceTypesMod    , only : hlm_freq_day
  use EDPftvarcon          , only : EDPftvarcon_inst
  use EDPftvarcon          , only : GetDecompyFrac
  use EDCohortDynamicsMod  , only : fuse_cohorts, sort_cohorts, insert_cohort
  use EDCohortDynamicsMod  , only : DeallocateCohort
  use EDTypesMod           , only : area_site => area
  use ChecksBalancesMod    , only : PatchMassStock
  use FatesLitterMod       , only : ncwd
  use FatesLitterMod       , only : ndcmpy
  use FatesLitterMod       , only : litter_type
  use EDTypesMod           , only : homogenize_seed_pfts
  use EDTypesMod           , only : n_dbh_bins, area, patchfusion_dbhbin_loweredges
  use EDtypesMod           , only : force_patchfuse_min_biomass
  use EDTypesMod           , only : maxPatchesPerSite
  use EDTypesMod           , only : maxPatchesPerSite_by_disttype  
  use EDTypesMod           , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod           , only : site_massbal_type
  use EDTypesMod           , only : site_fluxdiags_type
  use EDTypesMod           , only : min_patch_area
  use EDTypesMod           , only : min_patch_area_forced
  use EDTypesMod           , only : nclmax
  use EDTypesMod           , only : maxpft
  use EDTypesMod           , only : dtype_ifall
  use EDTypesMod           , only : dtype_ilog
  use EDTypesMod           , only : dtype_ifire
  use EDTypesMod           , only : ican_upper
  use EDTypesMod           , only : num_elements
  use EDTypesMod           , only : element_list
  use EDTypesMod           , only : element_pos
  use EDTypesMod           , only : lg_sf
  use EDTypesMod           , only : dl_sf
  use EDTypesMod           , only : dump_patch
  use FatesConstantsMod    , only : rsnbl_math_prec
  use FatesInterfaceTypesMod    , only : hlm_use_planthydro
  use FatesInterfaceTypesMod    , only : hlm_numSWb
  use FatesInterfaceTypesMod    , only : bc_in_type
  use FatesInterfaceTypesMod    , only : hlm_days_per_year
  use FatesInterfaceTypesMod    , only : numpft
  use FatesGlobals         , only : endrun => fates_endrun
  use FatesConstantsMod    , only : r8 => fates_r8
  use FatesConstantsMod    , only : itrue, ifalse
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use FatesPlantHydraulicsMod, only : AccumulateMortalityWaterStorage
  use FatesPlantHydraulicsMod, only : DeallocateHydrCohort
  use EDLoggingMortalityMod, only : logging_litter_fluxes 
  use EDLoggingMortalityMod, only : logging_time
  use EDParamsMod          , only : fates_mortality_disturbance_fraction
  use FatesAllometryMod    , only : carea_allom
  use FatesAllometryMod    , only : set_root_fraction
  use FatesConstantsMod    , only : g_per_kg
  use FatesConstantsMod    , only : ha_per_m2
  use FatesConstantsMod    , only : days_per_sec
  use FatesConstantsMod    , only : years_per_day
  use FatesConstantsMod    , only : nearzero
  use FatesConstantsMod    , only : primaryforest, secondaryforest
  use FatesConstantsMod    , only : n_anthro_disturbance_categories
  use FatesConstantsMod    , only : fates_unset_r8
  use FatesConstantsMod    , only : fates_unset_int
  use EDCohortDynamicsMod  , only : InitPRTObject
  use EDCohortDynamicsMod  , only : InitPRTBoundaryConditions
  use ChecksBalancesMod,      only : SiteMassStock
  use PRTGenericMod,          only : all_carbon_elements
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

  ! CIME globals
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  
  !
  implicit none
  private
  !
  public :: create_patch
  public :: spawn_patches
  public :: zero_patch
  public :: fuse_patches
  public :: terminate_patches
  public :: patch_pft_size_profile
  public :: disturbance_rates
  public :: check_patch_area
  public :: set_patchno
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
    ! loging flux
    use EDLoggingMortalityMod , only : LoggingMortality_frac

  
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout), target :: site_in
    type(bc_in_type) , intent(in) :: bc_in
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort

    real(r8) :: cmort
    real(r8) :: bmort
    real(r8) :: hmort
    real(r8) :: frmort
    real(r8) :: smort
    real(r8) :: asmort

    real(r8) :: lmort_direct
    real(r8) :: lmort_collateral
    real(r8) :: lmort_infra
    real(r8) :: l_degrad         ! fraction of trees that are not killed but suffer from forest 
                                 ! degradation (i.e. they are moved to newly-anthro-disturbed 
                                 ! secondary forest patch)
    real(r8) :: dist_rate_ldist_notharvested
    integer  :: threshold_sizeclass

    !----------------------------------------------------------------------------------------------
    ! Calculate Mortality Rates (these were previously calculated during growth derivatives)
    ! And the same rates in understory plants have already been applied to %dndt
    !----------------------------------------------------------------------------------------------
    
    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))   

       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))        
          ! Mortality for trees in the understorey.
          currentCohort%patchptr => currentPatch

          call mortality_rates(currentCohort,bc_in,cmort,hmort,bmort,frmort,smort,asmort)
          currentCohort%dmort  = cmort+hmort+bmort+frmort+smort+asmort
          call carea_allom(currentCohort%dbh,currentCohort%n,site_in%spread,currentCohort%pft, &
               currentCohort%c_area)

          ! Initialize diagnostic mortality rates
          currentCohort%cmort = cmort
          currentCohort%bmort = bmort
          currentCohort%hmort = hmort
          currentCohort%frmort = frmort
          currentCohort%smort = smort
          currentCohort%asmort = asmort

          call LoggingMortality_frac(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_layer, &
                lmort_direct,lmort_collateral,lmort_infra,l_degrad )
         
          currentCohort%lmort_direct     = lmort_direct
          currentCohort%lmort_collateral = lmort_collateral
          currentCohort%lmort_infra      = lmort_infra
          currentCohort%l_degrad         = l_degrad
          
          currentCohort => currentCohort%taller
       end do
       currentPatch%disturbance_mode = fates_unset_int
       currentPatch => currentPatch%younger
    end do

    ! ---------------------------------------------------------------------------------------------
    ! Calculate Disturbance Rates based on the mortality rates just calculated
    ! ---------------------------------------------------------------------------------------------

    currentPatch => site_in%oldest_patch
    do while (associated(currentPatch))   
       
       currentPatch%disturbance_rates(dtype_ifall) = 0.0_r8
       currentPatch%disturbance_rates(dtype_ilog)  = 0.0_r8
       currentPatch%disturbance_rates(dtype_ifire) = 0.0_r8

       dist_rate_ldist_notharvested = 0.0_r8
       
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))   

          if(currentCohort%canopy_layer == 1)then

             ! Treefall Disturbance Rate
             currentPatch%disturbance_rates(dtype_ifall) = currentPatch%disturbance_rates(dtype_ifall) + &
                  fates_mortality_disturbance_fraction * &
                  min(1.0_r8,currentCohort%dmort)*hlm_freq_day*currentCohort%c_area/currentPatch%area

             ! Logging Disturbance Rate
             currentPatch%disturbance_rates(dtype_ilog) = currentPatch%disturbance_rates(dtype_ilog) + &
                   min(1.0_r8, currentCohort%lmort_direct +                          & 
                               currentCohort%lmort_collateral +                      &
                               currentCohort%lmort_infra +                           &
                               currentCohort%l_degrad ) *                            &
                               currentCohort%c_area/currentPatch%area
             
             ! Non-harvested part of the logging disturbance rate
             dist_rate_ldist_notharvested = dist_rate_ldist_notharvested + currentCohort%l_degrad * &
                  currentCohort%c_area/currentPatch%area
             
          endif
          currentCohort => currentCohort%taller
       enddo !currentCohort

       ! fraction of the logging disturbance rate that is non-harvested
       if (currentPatch%disturbance_rates(dtype_ilog) .gt. nearzero) then
          currentPatch%fract_ldist_not_harvested = dist_rate_ldist_notharvested / &
               currentPatch%disturbance_rates(dtype_ilog)
       endif

       ! Fire Disturbance Rate
       ! Fires can't burn the whole patch, as this causes /0 errors. 
       currentPatch%disturbance_rates(dtype_ifire) = currentPatch%frac_burnt

       if (debug) then
          if (currentPatch%disturbance_rates(dtype_ifire) > 0.98_r8)then
          write(fates_log(),*) 'very high fire areas', &
               currentPatch%disturbance_rates(dtype_ifire),currentPatch%frac_burnt
          endif
       endif



       ! ------------------------------------------------------------------------------------------
       ! Determine which disturbance is dominant, and force mortality diagnostics in the upper 
       ! canopy to be zero for the non-dominant mode.  Note: upper-canopy tree-fall mortality is 
       ! not always disturbance generating, so when tree-fall mort is non-dominant, make sure
       ! to still diagnose and track the non-disturbance rate
       ! ------------------------------------------------------------------------------------------
       
       ! DISTURBANCE IS LOGGING
       if (currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifall) .and. &
             currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifire) ) then 
          
          currentPatch%disturbance_rate = currentPatch%disturbance_rates(dtype_ilog)
          currentPatch%disturbance_mode = dtype_ilog

          ! Update diagnostics
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             if(currentCohort%canopy_layer == 1)then
                currentCohort%cmort = currentCohort%cmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%hmort = currentCohort%hmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%bmort = currentCohort%bmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%dmort = currentCohort%dmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%frmort = currentCohort%frmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%smort = currentCohort%smort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%asmort = currentCohort%asmort*(1.0_r8 - fates_mortality_disturbance_fraction)
             end if
             currentCohort => currentCohort%taller
          enddo !currentCohort
          
       ! DISTURBANCE IS FIRE
       elseif (currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ifall) .and. &
             currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ilog) ) then  

          currentPatch%disturbance_rate = currentPatch%disturbance_rates(dtype_ifire)
          currentPatch%disturbance_mode = dtype_ifire

          ! Update diagnostics, zero non-fire mortality rates
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             if(currentCohort%canopy_layer == 1)then
                currentCohort%cmort = currentCohort%cmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%hmort = currentCohort%hmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%bmort = currentCohort%bmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%dmort = currentCohort%dmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%frmort = currentCohort%frmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%smort = currentCohort%smort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%asmort = currentCohort%asmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%lmort_direct    = 0.0_r8
                currentCohort%lmort_collateral = 0.0_r8
                currentCohort%lmort_infra      = 0.0_r8
                currentCohort%l_degrad         = 0.0_r8
             end if
 
             ! This may be counter-intuitive, but the diagnostic fire-mortality rate
             ! will stay zero in the patch that undergoes fire, this is because
             ! the actual cohorts who experience the fire are only those in the
             ! newly created patch so currentCohort%fmort = 0.0_r8
             ! Don't worry, the cohorts in the newly created patch will reflect burn

             currentCohort => currentCohort%taller
          enddo !currentCohort

       else  ! If fire and logging are not greater than treefall, just set disturbance rate to tree-fall
             ! which is most likely a 0.0

          currentPatch%disturbance_rate = currentPatch%disturbance_rates(dtype_ifall)
          currentPatch%disturbance_mode = dtype_ifall

          ! Update diagnostics, zero non-treefall mortality rates
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             if(currentCohort%canopy_layer == 1)then
                currentCohort%lmort_direct    = 0.0_r8
                currentCohort%lmort_collateral = 0.0_r8
                currentCohort%lmort_infra      = 0.0_r8
                currentCohort%l_degrad         = 0.0_r8
             end if
             currentCohort => currentCohort%taller
          enddo !currentCohort


       endif

       currentPatch => currentPatch%younger

    enddo !patch loop 

  end subroutine disturbance_rates

    ! ============================================================================
  subroutine spawn_patches( currentSite, bc_in)
    !
    ! !DESCRIPTION:
    ! In this subroutine, the following happens
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
    
    use EDParamsMod         , only : ED_val_understorey_death, logging_coll_under_frac
    use EDCohortDynamicsMod , only : zero_cohort, copy_cohort, terminate_cohorts 

    !
    ! !ARGUMENTS:
    type (ed_site_type), intent(inout), target :: currentSite
    type (bc_in_type), intent(in)              :: bc_in
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer :: new_patch
    type (ed_patch_type) , pointer :: new_patch_primary
    type (ed_patch_type) , pointer :: new_patch_secondary
    type (ed_patch_type) , pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    type (ed_cohort_type), pointer :: nc
    type (ed_cohort_type), pointer :: storesmallcohort
    type (ed_cohort_type), pointer :: storebigcohort
    real(r8) :: site_areadis_primary         ! total area disturbed (to primary forest) in m2 per site per day
    real(r8) :: site_areadis_secondary       ! total area disturbed (to secondary forest) in m2 per site per day    
    real(r8) :: patch_site_areadis           ! total area disturbed in m2 per patch per day
    real(r8) :: age                          ! notional age of this patch in years
    integer  :: el                           ! element loop index
    integer  :: tnull                        ! is there a tallest cohort?
    integer  :: snull                        ! is there a shortest cohort?
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
    !---------------------------------------------------------------------

    storesmallcohort => null() ! storage of the smallest cohort for insertion routine
    storebigcohort   => null() ! storage of the largest cohort for insertion routine 

    ! calculate area of disturbed land, in this timestep, by summing contributions from each existing patch. 
    currentPatch => currentSite%youngest_patch

    site_areadis_primary = 0.0_r8
    site_areadis_secondary = 0.0_r8    

    do while(associated(currentPatch))

    
       if(currentPatch%disturbance_rate>1.0_r8) then
          write(fates_log(),*) 'patch disturbance rate > 1 ?',currentPatch%disturbance_rate
          call endrun(msg=errMsg(sourcefile, __LINE__))          
       end if

       ! Check to make sure that the disturbance mode of the patch is set
       if( .not.any(currentPatch%disturbance_mode == [dtype_ilog,dtype_ifall,dtype_ifire])) then
           write(fates_log(),*) 'undefined disturbance mode? : ',currentPatch%disturbance_mode
           call endrun(msg=errMsg(sourcefile, __LINE__))    
       end if


       ! Only create new patches that have non-negligible amount of land
       if((currentPatch%area*currentPatch%disturbance_rate) > nearzero ) then
          
          ! figure out whether the receiver patch for disturbance from this patch will be 
          ! primary or secondary land receiver patch is primary forest only if both the 
          ! donor patch is primary forest and the dominant disturbance type is not logging
          if ( currentPatch%anthro_disturbance_label .eq. primaryforest .and. &
                (currentPatch%disturbance_mode .ne. dtype_ilog) ) then
             
             site_areadis_primary = site_areadis_primary + currentPatch%area * currentPatch%disturbance_rate
          else
             site_areadis_secondary = site_areadis_secondary + currentPatch%area * currentPatch%disturbance_rate          
          endif
          
       end if

       currentPatch => currentPatch%older     
    enddo ! end loop over patches. sum area disturbed for all patches. 

    ! It is possible that no disturbance area was generated
    if ( (site_areadis_primary + site_areadis_secondary) > nearzero) then  
       
       age = 0.0_r8

       ! create two empty patches, to absorb newly disturbed primary and secondary forest area
       ! first create patch to receive primary forest area
       if ( site_areadis_primary .gt. nearzero ) then
          allocate(new_patch_primary)

          call create_patch(currentSite, new_patch_primary, age, &
                site_areadis_primary, primaryforest)
          
          ! Initialize the litter pools to zero, these
          ! pools will be populated by looping over the existing patches
          ! and transfering in mass
          do el=1,num_elements
              call new_patch_primary%litter(el)%InitConditions(init_leaf_fines=0._r8, &
                    init_root_fines=0._r8, &
                    init_ag_cwd=0._r8, &
                    init_bg_cwd=0._r8, &
                    init_seed=0._r8,   &
                    init_seed_germ=0._r8)
          end do
          new_patch_primary%tallest  => null()
          new_patch_primary%shortest => null()

       endif


       ! next create patch to receive secondary forest area
       if ( site_areadis_secondary .gt. nearzero) then
          allocate(new_patch_secondary)
          call create_patch(currentSite, new_patch_secondary, age, &
                site_areadis_secondary, secondaryforest)
          
          ! Initialize the litter pools to zero, these
          ! pools will be populated by looping over the existing patches
          ! and transfering in mass
          do el=1,num_elements
              call new_patch_secondary%litter(el)%InitConditions(init_leaf_fines=0._r8, &
                    init_root_fines=0._r8, &
                    init_ag_cwd=0._r8, &
                    init_bg_cwd=0._r8, &
                    init_seed=0._r8,   &
                    init_seed_germ=0._r8)
          end do
          new_patch_secondary%tallest  => null()
          new_patch_secondary%shortest => null() 

       endif
    
       ! loop round all the patches that contribute surviving indivduals and litter 
       ! pools to the new patch.  We only loop the pre-existing patches, so 
       ! quit the loop if the current patch is either null, or matches the
       ! two new pointers.

       currentPatch => currentSite%oldest_patch
       do while(associated(currentPatch))

          ! This is the amount of patch area that is disturbed, and donated by the donor
          patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate

          
          if ( patch_site_areadis > nearzero ) then

              ! figure out whether the receiver patch for disturbance from this patch 
              ! will be primary or secondary land receiver patch is primary forest 
              ! only if both the donor patch is primary forest and the dominant 
              ! disturbance type is not logging
              if (currentPatch%anthro_disturbance_label .eq. primaryforest .and. &
                    (currentPatch%disturbance_mode .ne. dtype_ilog)) then
                  new_patch => new_patch_primary
              else
                  new_patch => new_patch_secondary
              endif
              
              if(.not.associated(new_patch))then
                  write(fates_log(),*) 'Patch spawning has attempted to point to'
                  write(fates_log(),*) 'an un-allocated patch'
                  call endrun(msg=errMsg(sourcefile, __LINE__)) 
              end if

             ! for the case where the donating patch is secondary forest, if 
             ! the dominant disturbance from this patch is non-anthropogenic,
             ! we need to average in the time-since-anthropogenic-disturbance 
             ! from the donor patch into that of the receiver patch
             if ( currentPatch%anthro_disturbance_label .eq. secondaryforest .and. &
                   (currentPatch%disturbance_mode .ne. dtype_ilog) ) then

                new_patch%age_since_anthro_disturbance = new_patch%age_since_anthro_disturbance + &
                     currentPatch%age_since_anthro_disturbance * (patch_site_areadis / site_areadis_secondary)

             endif
          
             
             ! Transfer the litter existing already in the donor patch to the new patch
             ! This call will only transfer non-burned litter to new patch
             ! and burned litter to atmosphere. Thus it is important to zero burnt_frac_litter when
             ! fire is not the dominant disturbance regime. 

             if(currentPatch%disturbance_mode .ne. dtype_ifire) then
                 currentPatch%burnt_frac_litter(:) = 0._r8
             end if

             call TransLitterNewPatch( currentSite, currentPatch, new_patch, patch_site_areadis)

             ! Transfer in litter fluxes from plants in various contexts of death and destruction

             if(currentPatch%disturbance_mode .eq. dtype_ilog) then
                call logging_litter_fluxes(currentSite, currentPatch, new_patch, patch_site_areadis)
             elseif(currentPatch%disturbance_mode .eq. dtype_ifire) then
                call fire_litter_fluxes(currentSite, currentPatch, new_patch, patch_site_areadis)  
             else
                call mortality_litter_fluxes(currentSite, currentPatch, new_patch, patch_site_areadis)
             endif

             ! --------------------------------------------------------------------------
             ! The newly formed patch from disturbance (new_patch), has now been given 
             ! some litter from dead plants and pre-existing litter from the donor patches.
             !
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
                 call InitPRTBoundaryConditions(nc)
                 
                 call zero_cohort(nc)

                 ! nc is the new cohort that goes in the disturbed patch (new_patch)... currentCohort
                 ! is the curent cohort that stays in the donor patch (currentPatch) 
                 call copy_cohort(currentCohort, nc)

                 !this is the case as the new patch probably doesn't have a closed canopy, and
                 ! even if it does, that will be sorted out in canopy_structure. 
                 nc%canopy_layer = 1 
                 nc%canopy_layer_yesterday = 1._r8 
                 
                 sapw_c   = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
                 struct_c = currentCohort%prt%GetState(struct_organ, all_carbon_elements)
                 leaf_c   = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
                 fnrt_c   = currentCohort%prt%GetState(fnrt_organ, all_carbon_elements)
                 store_c  = currentCohort%prt%GetState(store_organ, all_carbon_elements)
                 total_c  = sapw_c + struct_c + leaf_c + fnrt_c + store_c

                 ! treefall mortality is the dominant disturbance
                 if(currentPatch%disturbance_mode .eq. dtype_ifall) then
                   
                     if(currentCohort%canopy_layer == 1)then

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
                      nc%lmort_direct     = nan
                      nc%lmort_collateral = nan
                      nc%lmort_infra      = nan
                      nc%l_degrad         = nan
                      
                   else
                      ! small trees 
                      if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then
                         
                         
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
                         
                         
                         currentSite%imort_carbonflux = currentSite%imort_carbonflux + &
                              (nc%n * ED_val_understorey_death / hlm_freq_day ) * &
                              total_c * g_per_kg * days_per_sec * years_per_day * ha_per_m2
                         
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
                         nc%dmort            = currentCohort%dmort
                         nc%lmort_direct    = currentCohort%lmort_direct
                         nc%lmort_collateral = currentCohort%lmort_collateral
                         nc%lmort_infra      = currentCohort%lmort_infra
                         
                      endif
                   endif
                   
                   ! Fire is the dominant disturbance 
                elseif (currentPatch%disturbance_mode .eq. dtype_ifire ) then
                   
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
                      
                      currentSite%fmort_carbonflux_canopy = currentSite%fmort_carbonflux_canopy + &
                           (nc%n * currentCohort%fire_mort) * &
                           total_c * g_per_kg * days_per_sec * ha_per_m2
                      
                   else
                      currentSite%fmort_rate_ustory(currentCohort%size_class, currentCohort%pft) = &
                           currentSite%fmort_rate_ustory(currentCohort%size_class, currentCohort%pft) + &
                           nc%n * currentCohort%fire_mort / hlm_freq_day
                      
                      currentSite%fmort_carbonflux_ustory = currentSite%fmort_carbonflux_ustory + &
                           (nc%n * currentCohort%fire_mort) * &
                           total_c * g_per_kg * days_per_sec * ha_per_m2
                   end if
                   
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
                   nc%dmort            = currentCohort%dmort
                   nc%lmort_direct     = currentCohort%lmort_direct
                   nc%lmort_collateral = currentCohort%lmort_collateral
                   nc%lmort_infra      = currentCohort%lmort_infra


                   ! Some of of the leaf mass from living plants has been
                   ! burned off.  Here, we remove that mass, and
                   ! tally it in the flux we sent to the atmosphere
                   
                   if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then
                       leaf_burn_frac = currentCohort%fraction_crown_burned
                   else

                       ! Grasses determine their fraction of leaves burned here

                       leaf_burn_frac = currentPatch%burnt_frac_litter(lg_sf)
                   endif
                   
                   ! Perform a check to make sure that spitfire gave
                   ! us reasonable mortality and burn fraction rates
                   
                   if( (leaf_burn_frac < 0._r8) .or. &
                       (leaf_burn_frac > 1._r8) .or. &
                       (currentCohort%fire_mort < 0._r8) .or. &
                       (currentCohort%fire_mort > 1._r8)) then
                       write(fates_log(),*) 'unexpected fire fractions'
                       write(fates_log(),*) EDPftvarcon_inst%woody(currentCohort%pft)
                       write(fates_log(),*) leaf_burn_frac
                       write(fates_log(),*) currentCohort%fire_mort
                       call endrun(msg=errMsg(sourcefile, __LINE__))
                   end if

                   do el = 1,num_elements

                       leaf_m = nc%prt%GetState(leaf_organ, element_list(el))

                       currentSite%mass_balance(el)%burn_flux_to_atm = &
                             currentSite%mass_balance(el)%burn_flux_to_atm + & 
                             leaf_burn_frac * leaf_m * nc%n
                   end do

                   ! Here the mass is removed from the plant

                   call PRTBurnLosses(nc%prt, leaf_organ, leaf_burn_frac)
                   currentCohort%fraction_crown_burned = 0.0_r8
                   nc%fraction_crown_burned            = 0.0_r8



               ! Logging is the dominant disturbance  
               elseif (currentPatch%disturbance_mode .eq. dtype_ilog ) then
                   
                   ! If this cohort is in the upper canopy. It generated 
                   if(currentCohort%canopy_layer == 1)then
                      
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
                      nc%dmort            = currentCohort%dmort

                      ! since these are the ones that weren't logged, 
                      ! set the logging mortality rates as zero
                      nc%lmort_direct     = 0._r8
                      nc%lmort_collateral = 0._r8
                      nc%lmort_infra      = 0._r8
                      
                   else
                      
                      ! WHat to do with cohorts in the understory of a logging generated
                      ! disturbance patch?
                      
                      if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then
                         
                         
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

                         currentSite%imort_carbonflux = currentSite%imort_carbonflux + &
                              (nc%n * currentPatch%fract_ldist_not_harvested * &
                              logging_coll_under_frac/ hlm_freq_day ) * &
                              total_c * g_per_kg * days_per_sec * years_per_day * ha_per_m2
                         
                         
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
                         nc%dmort            = currentCohort%dmort
                         nc%lmort_direct     = currentCohort%lmort_direct
                         nc%lmort_collateral = currentCohort%lmort_collateral
                         nc%lmort_infra      = currentCohort%lmort_infra
                         
                      endif  ! is/is-not woody
                      
                   endif  ! Select canopy layer
                   
               else
                  write(fates_log(),*) 'unknown disturbance mode?'
                  write(fates_log(),*) 'disturbance_mode: ',currentPatch%disturbance_mode 
                  call endrun(msg=errMsg(sourcefile, __LINE__))          
               end if   ! Select disturbance mode
                
                if (nc%n > 0.0_r8) then   
                   storebigcohort   =>  new_patch%tallest
                   storesmallcohort =>  new_patch%shortest 
                   if(associated(new_patch%tallest))then
                      tnull = 0
                   else
                      tnull = 1
                      new_patch%tallest => nc
                      nc%taller => null()
                   endif
                   
                   if(associated(new_patch%shortest))then
                      snull = 0
                   else
                      snull = 1
                      new_patch%shortest => nc
                      nc%shorter => null()
                   endif
                   nc%patchptr => new_patch
                   call insert_cohort(nc, new_patch%tallest, new_patch%shortest, &
                         tnull, snull, storebigcohort, storesmallcohort)
                   
                   new_patch%tallest  => storebigcohort 
                   new_patch%shortest => storesmallcohort   
                else
                   
                   ! Get rid of the new temporary cohort
                   call DeallocateCohort(nc)
                   deallocate(nc)
                   
                endif
                
                currentCohort => currentCohort%taller      
             enddo ! currentCohort 
             call sort_cohorts(currentPatch)
             
             !update area of donor patch
             currentPatch%area = currentPatch%area - patch_site_areadis

             ! sort out the cohorts, since some of them may be so small as to need removing. 
             ! the first call to terminate cohorts removes sparse number densities,
             ! the second call removes for all other reasons (sparse culling must happen
             ! before fusion)
             call terminate_cohorts(currentSite, currentPatch, 1,16)
             call fuse_cohorts(currentSite,currentPatch, bc_in)
             call terminate_cohorts(currentSite, currentPatch, 2,16)
             call sort_cohorts(currentPatch)

          end if    ! if ( new_patch%area > nearzero ) then 
       
          !zero disturbance rate trackers
          currentPatch%disturbance_rate  = 0._r8
          currentPatch%disturbance_rates = 0._r8
          currentPatch%fract_ldist_not_harvested = 0._r8
          
          currentPatch => currentPatch%younger
          
      enddo ! currentPatch patch loop. 

      !*************************/
      !**  INSERT NEW PATCH(ES) INTO LINKED LIST    
      !**********`***************/
       
      if ( site_areadis_primary .gt. nearzero) then
          currentPatch               => currentSite%youngest_patch
          new_patch_primary%older    => currentPatch
          new_patch_primary%younger  => null()
          currentPatch%younger       => new_patch_primary
          currentSite%youngest_patch => new_patch_primary
      endif
      
      if ( site_areadis_secondary .gt. nearzero) then
          currentPatch               => currentSite%youngest_patch
          new_patch_secondary%older  => currentPatch
          new_patch_secondary%younger=> null()
          currentPatch%younger       => new_patch_secondary
          currentSite%youngest_patch => new_patch_secondary
      endif
 
       
       ! sort out the cohorts, since some of them may be so small as to need removing. 
       ! the first call to terminate cohorts removes sparse number densities,
       ! the second call removes for all other reasons (sparse culling must happen
       ! before fusion)

       if ( site_areadis_primary .gt. nearzero) then
          call terminate_cohorts(currentSite, new_patch_primary, 1,17)
          call fuse_cohorts(currentSite,new_patch_primary, bc_in)
          call terminate_cohorts(currentSite, new_patch_primary, 2,17)
          call sort_cohorts(new_patch_primary)
       endif
       
       if ( site_areadis_secondary .gt. nearzero) then
          call terminate_cohorts(currentSite, new_patch_secondary, 1,18)
          call fuse_cohorts(currentSite,new_patch_secondary, bc_in)
          call terminate_cohorts(currentSite, new_patch_secondary, 2,18)
          call sort_cohorts(new_patch_secondary)
       endif
       
    endif !end new_patch area 
    
    
    call check_patch_area(currentSite)
    call set_patchno(currentSite)
    
    return
  end subroutine spawn_patches

  ! ============================================================================

  subroutine check_patch_area( currentSite )
    !
    ! !DESCRIPTION:
    !  Check to see that total area is not exceeded.  
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(in), target  :: currentSite 
    !
    ! !LOCAL VARIABLES:
    real(r8)                     :: areatot
    type(ed_patch_type), pointer :: currentPatch 
    type(ed_patch_type), pointer :: largestPatch
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
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

       if(debug) then
          write(fates_log(),*) 'Total patch area precision being fixed, adjusting'
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
  subroutine set_patchno( currentSite )
    !
    ! !DESCRIPTION:
    !  Give patches an order number from the oldest to youngest. 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type),intent(in), target :: currentSite 
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: currentPatch 
    integer patchno
    !---------------------------------------------------------------------

    patchno = 1
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       currentPatch%patchno = patchno
       patchno = patchno + 1
       currentPatch => currentPatch%younger
    enddo

  end subroutine set_patchno

  ! ============================================================================

  subroutine TransLitterNewPatch(currentSite,        &
                                 currentPatch,       &
                                 newPatch,           &
                                 patch_site_areadis)

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
    type(ed_site_type)  , intent(in), target  :: currentSite        ! site
    type(ed_patch_type) , intent(in), target  :: currentPatch       ! Donor patch
    type(ed_patch_type) , intent(inout)       :: newPatch           ! New patch
    real(r8)            , intent(in)          :: patch_site_areadis ! Area being donated
                                                                    ! by current patch

    
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
             
          ! Transfer above ground CWD
          donatable_mass     = curr_litt%ag_cwd(c) * patch_site_areadis * &
                               (1._r8 - currentPatch%burnt_frac_litter(c))

          burned_mass        = curr_litt%ag_cwd(c) * patch_site_areadis * &
                               currentPatch%burnt_frac_litter(c)
 
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
          
       do dcmpy=1,ndcmpy

           ! Transfer leaf fines
           donatable_mass           = curr_litt%leaf_fines(dcmpy) * patch_site_areadis * &
                                      (1._r8 - currentPatch%burnt_frac_litter(dl_sf))

           burned_mass              = curr_litt%leaf_fines(dcmpy) * patch_site_areadis * &
                                      currentPatch%burnt_frac_litter(dl_sf)

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

  subroutine fire_litter_fluxes(currentSite, currentPatch, newPatch, patch_site_areadis)
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
    type(ed_patch_type) , intent(inout), target :: currentPatch   ! Donor Patch
    type(ed_patch_type) , intent(inout), target :: newPatch   ! New Patch
    real(r8)            , intent(in)            :: patch_site_areadis ! Area being donated
                                                                      ! by current patch
    !
    ! !LOCAL VARIABLES:

    type(ed_cohort_type), pointer      :: currentCohort
    type(litter_type), pointer         :: new_litt
    type(litter_type), pointer         :: curr_litt
    type(site_massbal_type), pointer   :: site_mass
    type(site_fluxdiags_type), pointer :: flux_diags

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
       flux_diags => currentSite%flux_diags(el)
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
             
             sapw_m   = currentCohort%prt%GetState(sapw_organ, element_id)
             struct_m = currentCohort%prt%GetState(struct_organ, element_id)
             leaf_m   = currentCohort%prt%GetState(leaf_organ, element_id)
             fnrt_m   = currentCohort%prt%GetState(fnrt_organ, element_id)
             store_m  = currentCohort%prt%GetState(store_organ, element_id)
             repro_m  = currentCohort%prt%GetState(repro_organ, element_id)
             
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

             call set_root_fraction(currentSite%rootfrac_scr, pft, currentSite%zi_soil)

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
             flux_diags%leaf_litter_input(pft) = &
                  flux_diags%leaf_litter_input(pft) + &
                  num_dead_trees * (leaf_m+repro_m) * (1.0_r8-currentCohort%fraction_crown_burned)

             flux_diags%root_litter_input(pft) = &
                  flux_diags%root_litter_input(pft) + &
                  (fnrt_m + store_m) * num_dead_trees
             
             ! coarse root biomass per tree
             bcroot = (sapw_m + struct_m) * (1.0_r8 - EDPftvarcon_inst%allom_agb_frac(pft) )
      
             ! below ground coarse woody debris from burned trees
             do c = 1,ncwd
                do sl = 1,currentSite%nlevsoil
                   donatable_mass =  num_dead_trees * SF_val_CWD_frac(c) * &
                         bcroot * currentSite%rootfrac_scr(sl)

                   new_litt%bg_cwd(c,sl) = new_litt%bg_cwd(c,sl) + &
                         donatable_mass * donate_m2
                   curr_litt%bg_cwd(c,sl) = curr_litt%bg_cwd(c,sl) + &
                         donatable_mass * retain_m2

                   ! track diagnostics
                   flux_diags%cwd_bg_input(c) = &
                         flux_diags%cwd_bg_input(c) + &
                         donatable_mass
                enddo
             end do

             ! stem biomass per tree
             bstem  = (sapw_m + struct_m) * EDPftvarcon_inst%allom_agb_frac(pft)

             ! Above ground coarse woody debris from twigs and small branches
             ! a portion of this pool may burn
             do c = 1,ncwd
                 donatable_mass = num_dead_trees * SF_val_CWD_frac(c) * bstem
                 if (c == 1 .or. c == 2) then
                      donatable_mass = donatable_mass * (1.0_r8-currentCohort%fraction_crown_burned)
                      burned_mass = num_dead_trees * SF_val_CWD_frac(c) * bstem * & 
                      currentCohort%fraction_crown_burned
                      site_mass%burn_flux_to_atm = site_mass%burn_flux_to_atm + burned_mass
                endif
                new_litt%ag_cwd(c) = new_litt%ag_cwd(c) + donatable_mass * donate_m2
                curr_litt%ag_cwd(c) = curr_litt%ag_cwd(c) + donatable_mass * retain_m2
                flux_diags%cwd_ag_input(c) = flux_diags%cwd_ag_input(c) + donatable_mass
             enddo   
                  

            currentCohort => currentCohort%taller
        enddo
    end do
    
    return
  end subroutine fire_litter_fluxes

  ! ============================================================================

  subroutine mortality_litter_fluxes(currentSite, currentPatch, newPatch, patch_site_areadis)
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
    type(ed_patch_type) , intent(inout), target :: currentPatch
    type(ed_patch_type) , intent(inout), target :: newPatch
    real(r8)            , intent(in)            :: patch_site_areadis

    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer      :: currentCohort
    type(litter_type), pointer         :: new_litt
    type(litter_type), pointer         :: curr_litt
    type(site_massbal_type), pointer   :: site_mass
    type(site_fluxdiags_type), pointer :: flux_diags

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
       flux_diags => currentSite%flux_diags(el)
       curr_litt  => currentPatch%litter(el)   ! Litter pool of "current" patch
       new_litt   => newPatch%litter(el)       ! Litter pool of "new" patch

       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))       

          pft = currentCohort%pft
   
          sapw_m   = currentCohort%prt%GetState(sapw_organ, element_id)
          struct_m = currentCohort%prt%GetState(struct_organ, element_id)
          leaf_m   = currentCohort%prt%GetState(leaf_organ, element_id)
          fnrt_m   = currentCohort%prt%GetState(fnrt_organ, element_id)
          store_m  = currentCohort%prt%GetState(store_organ, element_id)
          repro_m  = currentCohort%prt%GetState(repro_organ, element_id)

          if(currentCohort%canopy_layer == 1)then

             ! Upper canopy trees. The total dead is based on their disturbance
             ! generating mortality rate.
             
             num_dead = currentCohort%n * min(1.0_r8,currentCohort%dmort * &
                   hlm_freq_day * fates_mortality_disturbance_fraction)
             
          elseif(EDPftvarcon_inst%woody(pft) == itrue) then
             
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
          ag_wood = num_dead * (struct_m + sapw_m) * EDPftvarcon_inst%allom_agb_frac(pft)
          bg_wood = num_dead * (struct_m + sapw_m) * (1.0_r8-EDPftvarcon_inst%allom_agb_frac(pft))
          
          call set_root_fraction(currentSite%rootfrac_scr, pft, currentSite%zi_soil)


          do c=1,ncwd

             ! Transfer wood of dying trees to AG CWD pools
             new_litt%ag_cwd(c) = new_litt%ag_cwd(c) + ag_wood * &
                    SF_val_CWD_frac(c) * donate_m2

             curr_litt%ag_cwd(c) = curr_litt%ag_cwd(c) + ag_wood * &
                   SF_val_CWD_frac(c) * retain_m2
             
             ! Transfer wood of dying trees to BG CWD pools
             do sl = 1,currentSite%nlevsoil
                new_litt%bg_cwd(c,sl) = new_litt%bg_cwd(c,sl) + bg_wood * &
                       currentSite%rootfrac_scr(sl) * SF_val_CWD_frac(c) * &
                       donate_m2

                curr_litt%bg_cwd(c,sl) = curr_litt%bg_cwd(c,sl) + bg_wood * &
                      currentSite%rootfrac_scr(sl) * SF_val_CWD_frac(c) * &
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
             flux_diags%cwd_ag_input(c) = & 
                  flux_diags%cwd_ag_input(c) + SF_val_CWD_frac(c) * ag_wood
             
             flux_diags%cwd_bg_input(c) = &
                  flux_diags%cwd_bg_input(c) + SF_val_CWD_frac(c) * bg_wood
          end do

          flux_diags%leaf_litter_input(pft) = flux_diags%leaf_litter_input(pft) + &
               num_dead*(leaf_m + repro_m)

          flux_diags%root_litter_input(pft) = flux_diags%root_litter_input(pft) + & 
               num_dead * (fnrt_m + store_m*(1.0_r8-EDPftvarcon_inst%allom_frbstor_repro(pft)))
          

          
          currentCohort => currentCohort%taller      
       enddo !currentCohort         

    enddo


    return
  end subroutine mortality_litter_fluxes

  ! ============================================================================

  subroutine create_patch(currentSite, new_patch, age, areap, label)

    !
    ! !DESCRIPTION:
    !  Set default values for creating a new patch
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout), target :: currentSite
    type(ed_patch_type), intent(inout), target :: new_patch
    real(r8), intent(in) :: age                  ! notional age of this patch in years
    real(r8), intent(in) :: areap                ! initial area of this patch in m2. 
    integer, intent(in)  :: label                ! anthropogenic disturbance label

    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------
    integer :: el                                ! element loop index


    allocate(new_patch%tr_soil_dir(hlm_numSWb))
    allocate(new_patch%tr_soil_dif(hlm_numSWb))
    allocate(new_patch%tr_soil_dir_dif(hlm_numSWb))
    allocate(new_patch%fab(hlm_numSWb))
    allocate(new_patch%fabd(hlm_numSWb))
    allocate(new_patch%fabi(hlm_numSWb))
    allocate(new_patch%sabs_dir(hlm_numSWb))
    allocate(new_patch%sabs_dif(hlm_numSWb))


    ! Litter
    ! Allocate, Zero Fluxes, and Initialize to "unset" values

    allocate(new_patch%litter(num_elements))
    do el=1,num_elements
        call new_patch%litter(el)%InitAllocate(numpft,currentSite%nlevsoil,element_list(el))
        call new_patch%litter(el)%ZeroFlux()
        call new_patch%litter(el)%InitConditions(init_leaf_fines = fates_unset_r8, &
              init_root_fines = fates_unset_r8, &
              init_ag_cwd = fates_unset_r8, &
              init_bg_cwd = fates_unset_r8, &
              init_seed = fates_unset_r8,   &
              init_seed_germ = fates_unset_r8)
    end do

    call zero_patch(new_patch) !The nan value in here is not working??

    new_patch%tallest  => null() ! pointer to patch's tallest cohort    
    new_patch%shortest => null() ! pointer to patch's shortest cohort   
    new_patch%older    => null() ! pointer to next older patch   
    new_patch%younger  => null() ! pointer to next shorter patch      

    ! assign known patch attributes 

    new_patch%age                = age   
    new_patch%age_class          = 1
    new_patch%area               = areap 

    ! assign anthropgenic disturbance category and label
    new_patch%anthro_disturbance_label = label
    if (label .eq. secondaryforest) then
       new_patch%age_since_anthro_disturbance = age
    else
       new_patch%age_since_anthro_disturbance = -1._r8   ! replace with fates_unset_r8 when possible
    endif

    ! This new value will be generated when the calculate disturbance
    ! rates routine is called. This does not need to be remembered or in the restart file.
    new_patch%disturbance_mode   = fates_unset_int
 
    new_patch%f_sun              = 0._r8
    new_patch%ed_laisun_z(:,:,:) = 0._r8 
    new_patch%ed_laisha_z(:,:,:) = 0._r8 
    new_patch%ed_parsun_z(:,:,:) = 0._r8 
    new_patch%ed_parsha_z(:,:,:) = 0._r8 
    new_patch%fabi               = 0._r8
    new_patch%fabd               = 0._r8
    new_patch%tr_soil_dir(:)     = 1._r8
    new_patch%tr_soil_dif(:)     = 1._r8
    new_patch%tr_soil_dir_dif(:) = 0._r8
    new_patch%fabd_sun_z(:,:,:)  = 0._r8 
    new_patch%fabd_sha_z(:,:,:)  = 0._r8 
    new_patch%fabi_sun_z(:,:,:)  = 0._r8 
    new_patch%fabi_sha_z(:,:,:)  = 0._r8  
    new_patch%scorch_ht(:)       = 0._r8  
    new_patch%frac_burnt         = 0._r8  
    new_patch%total_tree_area    = 0.0_r8  
    new_patch%NCL_p              = 1

   
    return
  end subroutine create_patch

  ! ============================================================================
  subroutine zero_patch(cp_p)
    !
    ! !DESCRIPTION:
    !  Sets all the variables in the patch to nan or zero 
    ! (this needs to be two seperate routines, one for nan & one for zero
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_patch_type), intent(inout), target :: cp_p
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: currentPatch
    !---------------------------------------------------------------------

    currentPatch  => cp_p  

    currentPatch%tallest  => null()          
    currentPatch%shortest => null()         
    currentPatch%older    => null()               
    currentPatch%younger  => null()           

    currentPatch%patchno  = 999                            

    currentPatch%age                        = nan                          
    currentPatch%age_class                  = 1
    currentPatch%area                       = nan                                           
    currentPatch%canopy_layer_tlai(:)       = nan               
    currentPatch%total_canopy_area          = nan

    currentPatch%tlai_profile(:,:,:)        = nan 
    currentPatch%elai_profile(:,:,:)        = 0._r8 
    currentPatch%tsai_profile(:,:,:)        = nan 
    currentPatch%esai_profile(:,:,:)        = nan       
    currentPatch%canopy_area_profile(:,:,:) = nan       

    currentPatch%fabd_sun_z(:,:,:)          = nan 
    currentPatch%fabd_sha_z(:,:,:)          = nan 
    currentPatch%fabi_sun_z(:,:,:)          = nan 
    currentPatch%fabi_sha_z(:,:,:)          = nan  

    currentPatch%ed_laisun_z(:,:,:)         = nan 
    currentPatch%ed_laisha_z(:,:,:)         = nan 
    currentPatch%ed_parsun_z(:,:,:)         = nan 
    currentPatch%ed_parsha_z(:,:,:)         = nan 
    currentPatch%psn_z(:,:,:)               = 0._r8   

    currentPatch%f_sun(:,:,:)               = nan
    currentPatch%tr_soil_dir(:)             = nan    ! fraction of incoming direct  radiation that is transmitted to the soil as direct
    currentPatch%tr_soil_dif(:)             = nan    ! fraction of incoming diffuse radiation that is transmitted to the soil as diffuse
    currentPatch%tr_soil_dir_dif(:)         = nan    ! fraction of incoming direct  radiation that is transmitted to the soil as diffuse
    currentPatch%fabd(:)                    = nan    ! fraction of incoming direct  radiation that is absorbed by the canopy
    currentPatch%fabi(:)                    = nan    ! fraction of incoming diffuse radiation that is absorbed by the canopy

    currentPatch%canopy_mask(:,:)           = 999    ! is there any of this pft in this layer?
    currentPatch%nrad(:,:)                  = 999    ! number of exposed leaf layers for each canopy layer and pft
    currentPatch%ncan(:,:)                  = 999    ! number of total leaf layers for each canopy layer and pft
    currentPatch%pft_agb_profile(:,:)       = nan    

    ! DISTURBANCE 
    currentPatch%disturbance_rates          = 0._r8 
    currentPatch%disturbance_rate           = 0._r8
    currentPatch%fract_ldist_not_harvested  = 0._r8


    ! FIRE
   currentPatch%litter_moisture(:)         = 0.0_r8 ! litter moisture
    currentPatch%fuel_eff_moist             = 0.0_r8 ! average fuel moisture content of the ground fuel 
    ! (incl. live grasses. omits 1000hr fuels)
    currentPatch%livegrass                  = 0.0_r8 ! total ag grass biomass in patch. 1=c3 grass, 2=c4 grass. gc/m2
    currentPatch%sum_fuel                   = 0.0_r8 ! total ground fuel related to ros (omits 1000hr fuels). gc/m2
    currentPatch%fuel_bulkd                 = 0.0_r8 ! average fuel bulk density of the ground fuel 
    ! (incl. live grasses. omits 1000hr fuels). kgc/m3
    currentPatch%fuel_sav                   = 0.0_r8 ! average surface area to volume ratio of the ground fuel 
    ! (incl. live grasses. omits 1000hr fuels).
    currentPatch%fuel_mef                   = 0.0_r8 ! average moisture of extinction factor of the ground fuel
    ! (incl. live grasses. omits 1000hr fuels).
    currentPatch%ros_front                  = 0.0_r8 ! average rate of forward spread of each fire in the patch. m/min.
    currentPatch%effect_wspeed              = 0.0_r8 ! dailywind modified by fraction of relative grass and tree cover. m/min.
    currentPatch%tau_l                      = 0.0_r8 ! mins p&r(1986)
    currentPatch%fuel_frac(:)               = 0.0_r8 ! fraction of each litter class in the sum_fuel 
    !- for purposes of calculating weighted averages. 
    currentPatch%tfc_ros                    = 0.0_r8 ! used in fi calc
    currentPatch%fi                         = 0._r8  ! average fire intensity of flaming front during day.  
    ! backward ros plays no role. kj/m/s or kw/m.
    currentPatch%fire                       = 999    ! sr decide_fire.1=fire hot enough to proceed. 0=stop everything- no fires today
    currentPatch%fd                         = 0.0_r8 ! fire duration (mins)
    currentPatch%ros_back                   = 0.0_r8 ! backward ros (m/min)
    currentPatch%scorch_ht(:)               = 0.0_r8 ! scorch height of flames on a given PFT
    currentPatch%frac_burnt                 = 0.0_r8 ! fraction burnt daily  
    currentPatch%burnt_frac_litter(:)       = 0.0_r8 
    currentPatch%btran_ft(:)                = 0.0_r8

    currentPatch%canopy_layer_tlai(:)       = 0.0_r8

    currentPatch%fab(:)                     = 0.0_r8
    currentPatch%sabs_dir(:)                = 0.0_r8
    currentPatch%sabs_dif(:)                = 0.0_r8
    currentPatch%zstar                      = 0.0_r8
    currentPatch%c_stomata                  = 0.0_r8 ! This is calculated immediately before use
    currentPatch%c_lblayer                  = 0.0_r8

    currentPatch%solar_zenith_flag          = .false.
    currentPatch%solar_zenith_angle         = nan

    currentPatch%gnd_alb_dir(:)             = nan
    currentPatch%gnd_alb_dif(:)             = nan

  end subroutine zero_patch

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
    type(ed_patch_type), pointer :: currentPatch,tpp,tmpptr
    integer  :: ft,z        !counters for pft and height class
    real(r8) :: norm        !normalized difference between biomass profiles
    real(r8) :: profiletol  !tolerance of patch fusion routine. Starts off high and is reduced if there are too many patches.
    integer  :: nopatches(n_anthro_disturbance_categories)   !number of patches presently in gridcell
    integer  :: iterate     !switch of patch reduction iteration scheme. 1 to keep going, 0 to stop
    integer  :: fuse_flag   !do patches get fused (1) or not (0).
    integer  :: i_disttype  !iterator over anthropogenic disturbance categories
    !
    !---------------------------------------------------------------------

    currentSite => csite 

    profiletol = ED_val_patch_fusion_tol

    nopatches(1:n_anthro_disturbance_categories) = 0
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       nopatches(currentPatch%anthro_disturbance_label) = &
            nopatches(currentPatch%anthro_disturbance_label) + 1
       currentPatch => currentPatch%older
    enddo

    !---------------------------------------------------------------------!
    ! iterate over anthropogenic disturbance categories
    !---------------------------------------------------------------------!    

    do i_disttype = 1, n_anthro_disturbance_categories

       !---------------------------------------------------------------------!
       !  We only really care about fusing patches if nopatches > 1          !
       !---------------------------------------------------------------------!

       iterate = 1

       !---------------------------------------------------------------------!
       !  Keep doing this until nopatches <= maxPatchesPerSite               !
       !---------------------------------------------------------------------!

       do while(iterate == 1)
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
          do while(associated(currentPatch))      
             tpp => currentSite%youngest_patch
             do while(associated(tpp))

                if(.not.associated(currentPatch))then
                   write(fates_log(),*) 'ED: issue with currentPatch'
                endif

                if(associated(tpp).and.associated(currentPatch))then
                   !--------------------------------------------------------------------!
                   ! only fuse patches whose anthropogenic disturbance category matches !
                   ! that of the outer loop that we are in                              !
                   !--------------------------------------------------------------------!
                   if ( tpp%anthro_disturbance_label .eq. i_disttype .and. &
                        currentPatch%anthro_disturbance_label .eq. i_disttype) then

                      !--------------------------------------------------------------------------------------------
                      ! The default is to fuse the patches, unless some criteria is met which keeps them separated.
                      ! there are multiple criteria which all need to be met to keep them distinct:
                      ! (a) one of them is younger than the max age at which we force fusion;
                      ! (b) there is more than a threshold (tiny) amount of biomass in at least one of the patches;
                      ! (c) for at least one pft x size class, where there is biomass in that class in at least one patch,
                      ! and the normalized difference between the patches exceeds a threshold.
                      !--------------------------------------------------------------------------------------------

                      fuse_flag = 1
                      if(currentPatch%patchno /= tpp%patchno) then   !these should be the same patch

                         !-----------------------------------------------------------------------------------
                         ! check to see if both patches are older than the age at which we force them to fuse
                         !-----------------------------------------------------------------------------------

                         if ( tpp%age .le. max_age_of_second_oldest_patch .or. &
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

                            if(sum(currentPatch%pft_agb_profile(:,:)) > force_patchfuse_min_biomass .or. &
                                 sum(tpp%pft_agb_profile(:,:)) > force_patchfuse_min_biomass ) then

                               !---------------------------------------------------------------------!
                               ! Calculate the difference criteria for each pft and dbh class        !
                               !---------------------------------------------------------------------!   

                               do ft = 1,numpft        ! loop over pfts
                                  do z = 1,n_dbh_bins      ! loop over hgt bins 

                                     !----------------------------------
                                     ! is there biomass in this category?
                                     !----------------------------------

                                     if(currentPatch%pft_agb_profile(ft,z)  > 0.0_r8 .or.  &
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

                                        endif ! profile tol           
                                     endif ! biomass(ft,z) .gt. 0
                                  enddo !ht bins
                               enddo ! PFT
                            endif ! sum(biomass(:,:) .gt. force_patchfuse_min_biomass 
                         endif ! maxage

                         !-------------------------------------------------------------------------!
                         ! Call the patch fusion routine if there is not a meaningful difference   !
                         ! any of the pft x height categories                                      !
                         ! or both are older than forced fusion age                                !
                         !-------------------------------------------------------------------------!

                         if(fuse_flag  ==  1)then
                            
                            !-----------------------!
                            ! fuse the two patches  !
                            !-----------------------!
                            
                            tmpptr => currentPatch%older       
                            call fuse_2_patches(csite, currentPatch, tpp)
                            call fuse_cohorts(csite,tpp, bc_in)
                            call sort_cohorts(tpp)
                            currentPatch => tmpptr

                            !------------------------------------------------------------------------!
                            ! since we've just fused two patches, but are still in the midst of      !
                            ! a patch x patch loop, reset the patch fusion tolerance to the starting !
                            ! value so that any subsequent fusions in this loop are done with that   !
                            ! value. otherwise we can end up in a situation where we've loosened the !
                            ! fusion tolerance to get nopatches <= maxPatchesPerSite, but then,      !
                            ! having accomplished that, we continue through all the patch x patch    !
                            ! combinations and then all the patches get fused, ending up with        !
                            ! nopatches << maxPatchesPerSite and losing all heterogeneity.           !
                            !------------------------------------------------------------------------!

                            profiletol = ED_val_patch_fusion_tol
                            
                         else
                            ! write(fates_log(),*) 'patches not fused'
                         endif
                      endif  !are both patches the same anthropogenic disturbance category as the disturbance type loop iterator?
                   endif  !are both patches associated?        
                endif    !are these different patches?   
                tpp => tpp%older
             enddo !tpp loop

             if(associated(currentPatch))then 
                currentPatch => currentPatch%older 
             else
                currentPatch => null()
             endif !associated currentPatch

          enddo ! currentPatch loop

          !---------------------------------------------------------------------!
          ! Is the number of patches larger than the maximum?                   !
          !---------------------------------------------------------------------!   
          nopatches(i_disttype) = 0
          currentPatch => currentSite%youngest_patch
          do while(associated(currentPatch))
             if (currentPatch%anthro_disturbance_label .eq. i_disttype) then
                nopatches(i_disttype) = nopatches(i_disttype) +1
             endif
             currentPatch => currentPatch%older
          enddo

          if(nopatches(i_disttype) > maxPatchesPerSite_by_disttype(i_disttype))then
             iterate = 1
             profiletol = profiletol * patch_fusion_tolerance_relaxation_increment

             !---------------------------------------------------------------------!
             ! Making profile tolerance larger means that more fusion will happen  !
             !---------------------------------------------------------------------!        
          else
             iterate = 0
          endif

       enddo !do while nopatches>maxPatchesPerSite

    end do  ! i_disttype loop
 
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
    type (ed_patch_type) , pointer :: dp                ! Donor Patch
    type (ed_patch_type) , target, intent(inout) :: rp  ! Recipient Patch

    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type), pointer :: currentCohort ! Current Cohort
    type (ed_cohort_type), pointer :: nextc         ! Remembers next cohort in list 
    type (ed_cohort_type), pointer :: storesmallcohort
    type (ed_cohort_type), pointer :: storebigcohort  
    integer                        :: c,p          !counters for pft and litter size class. 
    integer                        :: tnull,snull  ! are the tallest and shortest cohorts associated?
    integer                        :: el           ! loop counting index for elements
    type(ed_patch_type), pointer   :: youngerp     ! pointer to the patch younger than donor
    type(ed_patch_type), pointer   :: olderp       ! pointer to the patch older than donor
    real(r8)                       :: inv_sum_area ! Inverse of the sum of the two patches areas
    !-----------------------------------------------------------------------------------------------

    ! Generate a litany of area weighted averages

    inv_sum_area = 1.0_r8/(dp%area + rp%area)
    
    rp%age = (dp%age * dp%area + rp%age * rp%area) * inv_sum_area

    rp%age_class = get_age_class_index(rp%age)
    
    do el = 1,num_elements
       call rp%litter(el)%FuseLitter(rp%area,dp%area,dp%litter(el))
    end do

    
    rp%fuel_eff_moist       = (dp%fuel_eff_moist*dp%area + rp%fuel_eff_moist*rp%area) * inv_sum_area
    rp%livegrass            = (dp%livegrass*dp%area + rp%livegrass*rp%area) * inv_sum_area
    rp%sum_fuel             = (dp%sum_fuel*dp%area + rp%sum_fuel*rp%area) * inv_sum_area
    rp%fuel_bulkd           = (dp%fuel_bulkd*dp%area + rp%fuel_bulkd*rp%area) * inv_sum_area
    rp%fuel_sav             = (dp%fuel_sav*dp%area + rp%fuel_sav*rp%area) * inv_sum_area
    rp%fuel_mef             = (dp%fuel_mef*dp%area + rp%fuel_mef*rp%area) * inv_sum_area
    rp%ros_front            = (dp%ros_front*dp%area + rp%ros_front*rp%area) * inv_sum_area
    rp%effect_wspeed        = (dp%effect_wspeed*dp%area + rp%effect_wspeed*rp%area) * inv_sum_area
    rp%tau_l                = (dp%tau_l*dp%area + rp%tau_l*rp%area) * inv_sum_area
    rp%fuel_frac(:)         = (dp%fuel_frac(:)*dp%area + rp%fuel_frac(:)*rp%area) * inv_sum_area
    rp%tfc_ros              = (dp%tfc_ros*dp%area + rp%tfc_ros*rp%area) * inv_sum_area
    rp%fi                   = (dp%fi*dp%area + rp%fi*rp%area) * inv_sum_area
    rp%fd                   = (dp%fd*dp%area + rp%fd*rp%area) * inv_sum_area
    rp%ros_back             = (dp%ros_back*dp%area + rp%ros_back*rp%area) * inv_sum_area
    rp%scorch_ht(:)         = (dp%scorch_ht(:)*dp%area + rp%scorch_ht(:)*rp%area) * inv_sum_area
    rp%frac_burnt           = (dp%frac_burnt*dp%area + rp%frac_burnt*rp%area) * inv_sum_area
    rp%burnt_frac_litter(:) = (dp%burnt_frac_litter(:)*dp%area + rp%burnt_frac_litter(:)*rp%area) * inv_sum_area
    rp%btran_ft(:)          = (dp%btran_ft(:)*dp%area + rp%btran_ft(:)*rp%area) * inv_sum_area
    rp%zstar                = (dp%zstar*dp%area + rp%zstar*rp%area) * inv_sum_area
    rp%c_stomata            = (dp%c_stomata*dp%area + rp%c_stomata*rp%area) * inv_sum_area
    rp%c_lblayer            = (dp%c_lblayer*dp%area + rp%c_lblayer*rp%area) * inv_sum_area
    
    rp%area = rp%area + dp%area !THIS MUST COME AT THE END!

    !insert donor cohorts into recipient patch
    if(associated(dp%shortest))then

       currentCohort => dp%shortest
       if(associated(currentCohort)) then
          nextc => currentCohort%taller
       endif

       do while(associated(dp%shortest))

          storebigcohort   => rp%tallest
          storesmallcohort => rp%shortest

          if(associated(rp%tallest))then
             tnull = 0
          else
             tnull = 1
             rp%tallest => currentCohort
          endif

          if(associated(rp%shortest))then
             snull = 0
          else
             snull = 1
             rp%shortest => currentCohort
          endif

          call insert_cohort(currentCohort, rp%tallest, rp%shortest, tnull, snull, storebigcohort, storesmallcohort)

          rp%tallest  => storebigcohort 
          rp%shortest => storesmallcohort    

          currentCohort%patchptr => rp

          currentCohort => nextc

          dp%shortest => currentCohort

          if(associated(currentCohort)) then
             nextc => currentCohort%taller
          endif

       enddo !cohort
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
    call dealloc_patch(dp)
    deallocate(dp)


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


  end subroutine fuse_2_patches

  ! ============================================================================

  subroutine terminate_patches(currentSite)
    !
    ! !DESCRIPTION:
    !  Terminate Patches if they  are too small                          
    !
    !
    ! !ARGUMENTS:
    type(ed_site_type), target, intent(inout) :: currentSite
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: currentPatch
    type(ed_patch_type), pointer :: olderPatch
    type(ed_patch_type), pointer :: youngerPatch
    integer, parameter           :: max_cycles = 10  ! After 10 loops through
                                                     ! You should had fused
    integer                      :: count_cycles

    real(r8) areatot ! variable for checking whether the total patch area is wrong. 
    !---------------------------------------------------------------------
 
    count_cycles = 0

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch)) 
       
       if(currentPatch%area <= min_patch_area)then
          
          ! Even if the patch area is small, avoid fusing it into its neighbor
          ! if it is the youngest of all patches. We do this in attempts to maintain
          ! a discrete patch for very young patches
          ! However, if the patch to be fused is excessivlely small, then fuse
          ! at all costs.  If it is not fused, it will make

          if ( .not.associated(currentPatch,currentSite%youngest_patch) .or. &
               currentPatch%area <= min_patch_area_forced ) then
             
             if(associated(currentPatch%older) )then
                
                if(debug) &
                     write(fates_log(),*) 'fusing to older patch because this one is too small',&
                     currentPatch%area, &
                     currentPatch%older%area
                
                ! We set a pointer to this patch, because
                ! it will be returned by the subroutine as de-referenced
                
                olderPatch => currentPatch%older
                call fuse_2_patches(currentSite, olderPatch, currentPatch)
                
                ! The fusion process has updated the "older" pointer on currentPatch
                ! for us.
                
                ! This logic checks to make sure that the younger patch is not the youngest
                ! patch. As mentioned earlier, we try not to fuse it.
                
             elseif( associated(currentPatch%younger) ) then
                
                if(debug) &
                      write(fates_log(),*) 'fusing to younger patch because oldest one is too small', &
                      currentPatch%area

                youngerPatch => currentPatch%younger
                call fuse_2_patches(currentSite, youngerPatch, currentPatch)
                
                ! The fusion process has updated the "younger" pointer on currentPatch
                
             endif
          endif
       endif
       
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
          call endrun(msg=errMsg(sourcefile, __LINE__))
          
          ! Note to user. If you DO decide to remove the end-run above this line
          ! Make sure that you keep the pointer below this line, or you will get
          ! an infinite loop.
          currentPatch => currentPatch%older
          count_cycles = 0
       end if

    enddo
    
    !check area is not exceeded
    call check_patch_area( currentSite )

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
      type(ed_patch_type), pointer              :: currentPatch
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


  ! =====================================================================================

  subroutine dealloc_patch(cpatch)

    ! This Subroutine is intended to de-allocate the allocatable memory that is pointed
    ! to via the patch structure.  This subroutine DOES NOT deallocate the patch
    ! structure itself.

    type(ed_patch_type), target :: cpatch

    type(ed_cohort_type), pointer :: ccohort  ! current
    type(ed_cohort_type), pointer :: ncohort  ! next
    integer                       :: el       ! loop counter for elements
    
    ! First Deallocate the cohort space
    ! -----------------------------------------------------------------------------------
    ccohort => cpatch%shortest
    do while(associated(ccohort))
       
       ncohort => ccohort%taller

       call DeallocateCohort(ccohort)
       deallocate(ccohort)
       ccohort => ncohort

    end do

    ! Deallocate all litter objects
    do el=1,num_elements
       call cpatch%litter(el)%DeallocateLitt()
    end do
    deallocate(cpatch%litter)

    ! Secondly, and lastly, deallocate the allocatable vector spaces in the patch
    if(allocated(cpatch%tr_soil_dir))then
       deallocate(cpatch%tr_soil_dir)
       deallocate(cpatch%tr_soil_dif)
       deallocate(cpatch%tr_soil_dir_dif)
       deallocate(cpatch%fab)
       deallocate(cpatch%fabd)
       deallocate(cpatch%fabi)
       deallocate(cpatch%sabs_dir)
       deallocate(cpatch%sabs_dif)
      
    end if

    return
  end subroutine dealloc_patch

  ! ============================================================================
  subroutine patch_pft_size_profile(cp_pnt)
    !
    ! !DESCRIPTION:
    !  Binned patch size profiles generated for patch fusion routine        
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_patch_type), target, intent(inout) :: cp_pnt
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type) , pointer  :: currentPatch
    type(ed_cohort_type), pointer  :: currentCohort
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
                  currentCohort%prt%GetState(struct_organ, all_carbon_elements) * &
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
    type (ed_patch_type), pointer :: currentPatch
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

 end module EDPatchDynamicsMod
