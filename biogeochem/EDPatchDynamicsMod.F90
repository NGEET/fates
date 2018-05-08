module EDPatchDynamicsMod

  ! ============================================================================
  ! Controls formation, creation, fusing and termination of patch level processes. 
  ! ============================================================================
  use FatesGlobals         , only : fates_log 
  use FatesInterfaceMod    , only : hlm_freq_day
  use EDPftvarcon          , only : EDPftvarcon_inst
  use EDCohortDynamicsMod  , only : fuse_cohorts, sort_cohorts, insert_cohort
  use EDtypesMod           , only : ncwd, n_dbh_bins, area, patchfusion_dbhbin_loweredges
  use EDtypesMod           , only : force_patchfuse_min_biomass
  use EDTypesMod           , only : maxPatchesPerSite
  use EDTypesMod           , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod           , only : min_patch_area
  use EDTypesMod           , only : nclmax
  use EDTypesMod           , only : maxpft
  use EDTypesMod           , only : dtype_ifall
  use EDTypesMod           , only : dtype_ilog
  use EDTypesMod           , only : dtype_ifire
  use FatesInterfaceMod    , only : hlm_use_planthydro
  use FatesInterfaceMod    , only : hlm_numlevgrnd
  use FatesInterfaceMod    , only : hlm_numlevsoil
  use FatesInterfaceMod    , only : hlm_numSWb
  use FatesInterfaceMod    , only : bc_in_type
  use FatesInterfaceMod    , only : hlm_days_per_year
  use FatesInterfaceMod    , only : numpft
  use FatesGlobals         , only : endrun => fates_endrun
  use FatesConstantsMod    , only : r8 => fates_r8
  use FatesConstantsMod    , only : itrue
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use FatesPlantHydraulicsMod, only : DeallocateHydrCohort
  use EDLoggingMortalityMod, only : logging_litter_fluxes 
  use EDLoggingMortalityMod, only : logging_time
  use EDParamsMod          , only : fates_mortality_disturbance_fraction
  use FatesAllometryMod    , only : carea_allom
  use FatesConstantsMod    , only : g_per_kg
  use FatesConstantsMod    , only : ha_per_m2
  use FatesConstantsMod    , only : days_per_sec
  use FatesConstantsMod    , only : years_per_day


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
  public :: set_root_fraction
  private:: fuse_2_patches

  character(len=*), parameter, private :: sourcefile = &
        __FILE__

  logical, parameter :: debug = .false.

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

    real(r8) :: lmort_direct
    real(r8) :: lmort_collateral
    real(r8) :: lmort_infra

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

          call mortality_rates(currentCohort,bc_in,cmort,hmort,bmort,frmort)
          currentCohort%dmort  = cmort+hmort+bmort+frmort
          call carea_allom(currentCohort%dbh,currentCohort%n,site_in%spread,currentCohort%pft, &
               currentCohort%c_area)

          ! Initialize diagnostic mortality rates
          currentCohort%cmort = cmort
          currentCohort%bmort = bmort
          currentCohort%hmort = hmort
          currentCohort%frmort = frmort
          currentCohort%fmort = 0.0_r8 ! Fire mortality is initialized as zero, but may be changed

          call LoggingMortality_frac(currentCohort%pft, currentCohort%dbh, &
                lmort_direct,lmort_collateral,lmort_infra )
         
          currentCohort%lmort_direct    = lmort_direct
          currentCohort%lmort_collateral = lmort_collateral
          currentCohort%lmort_infra      = lmort_infra

          
          currentCohort => currentCohort%taller
       end do
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
       
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))   

          if(currentCohort%canopy_layer == 1)then

             ! Treefall Disturbance Rate
             currentPatch%disturbance_rates(dtype_ifall) = currentPatch%disturbance_rates(dtype_ifall) + &
                  fates_mortality_disturbance_fraction * &
                  min(1.0_r8,currentCohort%dmort)*hlm_freq_day*currentCohort%c_area/currentPatch%area

             ! Logging Disturbance Rate
             currentPatch%disturbance_rates(dtype_ilog) = currentPatch%disturbance_rates(dtype_ilog) + &
                   min(1.0_r8, currentCohort%lmort_direct +                         & 
                               currentCohort%lmort_collateral +                      &
                               currentCohort%lmort_infra ) *                         &
                               currentCohort%c_area/currentPatch%area
             
          endif
          currentCohort => currentCohort%taller
       enddo !currentCohort

       ! Fire Disturbance Rate
       ! Fudge - fires can't burn the whole patch, as this causes /0 errors.
       ! This is accumulating the daily fires over the whole 30 day patch generation phase.  
       currentPatch%disturbance_rates(dtype_ifire) = &
             min(0.99_r8,currentPatch%disturbance_rates(dtype_ifire) + currentPatch%frac_burnt)

       if (currentPatch%disturbance_rates(dtype_ifire) > 0.98_r8)then
          write(fates_log(),*) 'very high fire areas', &
                currentPatch%disturbance_rates(dtype_ifire),currentPatch%frac_burnt
       endif



       ! ------------------------------------------------------------------------------------------
       ! Determine which disturbance is dominant, and force mortality diagnostics in the upper 
       ! canopy to be zero for the non-dominant mode.  Note: upper-canopy tree-fall mortality is 
       ! not always disturbance generating, so when tree-fall mort is non-dominant, make sure
       ! to still diagnose and track the non-disturbance rate
       ! ------------------------------------------------------------------------------------------
       
       
       if (currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifall) .and. &
             currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifire) ) then 
          
          currentPatch%disturbance_rate = currentPatch%disturbance_rates(dtype_ilog)

          ! Update diagnostics
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             if(currentCohort%canopy_layer == 1)then
                currentCohort%fmort = 0.0_r8
                currentCohort%cmort = currentCohort%cmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%hmort = currentCohort%hmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%bmort = currentCohort%bmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%dmort = currentCohort%dmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%frmort = currentCohort%frmort*(1.0_r8 - fates_mortality_disturbance_fraction)
             end if
             currentCohort => currentCohort%taller
          enddo !currentCohort
          
          
       elseif (currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ifall) .and. &
             currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ilog) ) then  ! DISTURBANCE IS FIRE

          currentPatch%disturbance_rate = currentPatch%disturbance_rates(dtype_ifire)

          ! Update diagnostics, zero non-fire mortality rates
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             if(currentCohort%canopy_layer == 1)then
                currentCohort%cmort = currentCohort%cmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%hmort = currentCohort%hmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%bmort = currentCohort%bmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%dmort = currentCohort%dmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%frmort = currentCohort%frmort*(1.0_r8 - fates_mortality_disturbance_fraction)
                currentCohort%lmort_direct    = 0.0_r8
                currentCohort%lmort_collateral = 0.0_r8
                currentCohort%lmort_infra      = 0.0_r8
             end if
 
             ! This may be counter-intuitive, but the diagnostic fire-mortality rate
             ! will stay zero in the patch that undergoes fire, this is because
             ! the actual cohorts who experience the fire are only those in the
             ! newly created patch so currentCohort%fmort = 0.0_r8
             ! Don't worry, the cohorts in the newly created patch will reflect burn

             currentCohort => currentCohort%taller
          enddo !currentCohort

       else  ! If fire and loggin are not greater than treefall, just set disturbance rate to tree-fall
             ! which is most likely a 0.0

          currentPatch%disturbance_rate = currentPatch%disturbance_rates(dtype_ifall)
          
          ! Update diagnostics, zero non-treefall mortality rates
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             if(currentCohort%canopy_layer == 1)then
                currentCohort%lmort_direct    = 0.0_r8
                currentCohort%lmort_collateral = 0.0_r8
                currentCohort%lmort_infra      = 0.0_r8
                currentCohort%fmort            = 0.0_r8
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
    type (ed_patch_type) , pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    type (ed_cohort_type), pointer :: nc
    type (ed_cohort_type), pointer :: storesmallcohort
    type (ed_cohort_type), pointer :: storebigcohort  
    real(r8) :: site_areadis                 ! total area disturbed in m2 per site per day
    real(r8) :: patch_site_areadis           ! total area disturbed in m2 per patch per day
    real(r8) :: age                          ! notional age of this patch in years
    integer  :: tnull                        ! is there a tallest cohort?
    integer  :: snull                        ! is there a shortest cohort?
    real(r8) :: root_litter_local(maxpft)    ! initial value of root litter. KgC/m2
    real(r8) :: leaf_litter_local(maxpft)    ! initial value of leaf litter. KgC/m2
    real(r8) :: cwd_ag_local(ncwd)           ! initial value of above ground coarse woody debris. KgC/m2
    real(r8) :: cwd_bg_local(ncwd)           ! initial value of below ground coarse woody debris. KgC/m2
    !---------------------------------------------------------------------

    storesmallcohort => null() ! storage of the smallest cohort for insertion routine
    storebigcohort   => null() ! storage of the largest cohort for insertion routine 

    ! calculate area of disturbed land, in this timestep, by summing contributions from each existing patch. 
    currentPatch => currentSite%youngest_patch

    ! zero site-level fire fluxes
    currentSite%cwd_ag_burned       = 0.0_r8
    currentSite%leaf_litter_burned  = 0.0_r8
    currentSite%total_burn_flux_to_atm = 0.0_r8    

    site_areadis = 0.0_r8
    do while(associated(currentPatch))

       !FIX(RF,032414) Does using the max(fire,mort) actually make sense here?
       site_areadis = site_areadis + currentPatch%area * min(1.0_r8,currentPatch%disturbance_rate) 
       currentPatch => currentPatch%older     

    enddo ! end loop over patches. sum area disturbed for all patches. 

    if (site_areadis > 0.0_r8) then  
       cwd_ag_local = 0.0_r8
       cwd_bg_local = 0.0_r8
       leaf_litter_local = 0.0_r8
       root_litter_local = 0.0_r8
       age = 0.0_r8

       allocate(new_patch)
       call create_patch(currentSite, new_patch, age, site_areadis, &
            cwd_ag_local, cwd_bg_local, leaf_litter_local, &
            root_litter_local)

       new_patch%tallest  => null()
       new_patch%shortest => null()

       currentPatch => currentSite%oldest_patch
       ! loop round all the patches that contribute surviving indivduals and litter pools to the new patch.     
       do while(associated(currentPatch))   

          ! This is the amount of patch area that is disturbed, and donated by the donor
          patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate

          call average_patch_properties(currentPatch, new_patch, patch_site_areadis)
          
          if (currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifall) .and. &
                currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifire) ) then 
             
             call logging_litter_fluxes(currentSite, currentPatch, new_patch, patch_site_areadis)
             
          elseif (currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ifall) .and. &
                currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ilog) ) then
             
             call fire_litter_fluxes(currentSite, currentPatch, new_patch, patch_site_areadis)  
             
          else
             
             call mortality_litter_fluxes(currentSite, currentPatch, new_patch, patch_site_areadis)
             
          endif

          !INSERT SURVIVORS FROM DISTURBANCE INTO NEW PATCH 
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))       

             allocate(nc)             
             if(hlm_use_planthydro.eq.itrue) call InitHydrCohort(nc)
             call zero_cohort(nc)

             ! nc is the new cohort that goes in the disturbed patch (new_patch)... currentCohort
             ! is the curent cohort that stays in the donor patch (currentPatch) 
             call copy_cohort(currentCohort, nc)

             !this is the case as the new patch probably doesn't have a closed canopy, and
             ! even if it does, that will be sorted out in canopy_structure. 
             nc%canopy_layer = 1 
             nc%canopy_layer_yesterday = 1._r8 

             ! treefall mortality is the dominant disturbance
             if(currentPatch%disturbance_rates(dtype_ifall) > currentPatch%disturbance_rates(dtype_ifire) .and. &
                    currentPatch%disturbance_rates(dtype_ifall) > currentPatch%disturbance_rates(dtype_ilog))then 

                if(currentCohort%canopy_layer == 1)then

                   ! In the donor patch we are left with fewer trees because the area has decreased
                   ! the plant density for large trees does not actually decrease in the donor patch
                   ! because this is the part of the original patch where no trees have actually fallen
                   ! The diagnostic cmort,bmort,hmort, and frmort  rates have already been saved         

                   currentCohort%n = currentCohort%n * (1.0_r8 - fates_mortality_disturbance_fraction * &
                        min(1.0_r8,currentCohort%dmort * hlm_freq_day))

                   nc%n = 0.0_r8      ! kill all of the trees who caused the disturbance.  
       
                   nc%cmort = nan     ! The mortality diagnostics are set to nan because the cohort should dissappear
                   nc%hmort = nan
                   nc%bmort = nan
                   nc%fmort = nan
                   nc%frmort = nan
                   nc%lmort_direct     = nan
                   nc%lmort_collateral = nan
                   nc%lmort_infra      = nan

                else
                   ! small trees 
                   if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then


                      ! Survivorship of undestory woody plants.  Two step process.
                      ! Step 1:  Reduce current number of plants to reflect the change in area.
                      !          The number density per square are doesn't change, but since the patch is smaller
                      !          and cohort counts are absolute, reduce this number.
                      nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                      ! because the mortality rate due to impact for the cohorts which had been in the understory and are now in the newly-
                      ! disturbed patch is very high, passing the imort directly to history results in large numerical errors, on account
                      ! of the sharply reduced number densities.  so instead pass this info via a site-level diagnostic variable before reducing 
                      ! the number density.
                      currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) = &
                           currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) + &
                           nc%n * ED_val_understorey_death / hlm_freq_day
                      currentSite%imort_carbonflux = currentSite%imort_carbonflux + &
                           (nc%n * ED_val_understorey_death / hlm_freq_day ) * &
                           currentCohort%b_total() * g_per_kg * days_per_sec * years_per_day * ha_per_m2
                      
                      ! Step 2:  Apply survivor ship function based on the understory death fraction
                      ! remaining of understory plants of those that are knocked over by the overstorey trees dying...  
                      nc%n = nc%n * (1.0_r8 - ED_val_understorey_death)
                      
                      ! since the donor patch split and sent a fraction of its members
                      ! to the new patch and a fraction to be preserved in itself,
                      ! when reporting diagnostic rates, we must carry over the mortality rates from
                      ! the donor that were applied before the patch split.  Remember this is only
                      ! for diagnostics.  But think of it this way, the rates are weighted by 
                      ! number density in EDCLMLink, and the number density of this new patch is donated
                      ! so with the number density must come the effective mortality rates.

                      nc%fmort            = 0.0_r8               ! Should had also been zero in the donor
                      nc%cmort            = currentCohort%cmort
                      nc%hmort            = currentCohort%hmort
                      nc%bmort            = currentCohort%bmort
                      nc%frmort           = currentCohort%frmort
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

                      nc%fmort            = 0.0_r8
                      nc%cmort            = currentCohort%cmort
                      nc%hmort            = currentCohort%hmort
                      nc%bmort            = currentCohort%bmort
                      nc%frmort           = currentCohort%frmort
                      nc%dmort            = currentCohort%dmort
                      nc%lmort_direct    = currentCohort%lmort_direct
                      nc%lmort_collateral = currentCohort%lmort_collateral
                      nc%lmort_infra      = currentCohort%lmort_infra
                      
                   endif
                endif

             ! Fire is the dominant disturbance 
             elseif (currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ifall) .and. &
                     currentPatch%disturbance_rates(dtype_ifire) > currentPatch%disturbance_rates(dtype_ilog)) then !fire

                ! Number of members in the new patch, before we impose fire survivorship
                nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                ! loss of individuals from source patch due to area shrinking
                currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area) 

                ! loss of individual from fire in new patch.
                nc%n = nc%n * (1.0_r8 - currentCohort%fire_mort) 

                nc%fmort            = currentCohort%fire_mort/hlm_freq_day
                
                nc%cmort            = currentCohort%cmort
                nc%hmort            = currentCohort%hmort
                nc%bmort            = currentCohort%bmort
                nc%frmort           = currentCohort%frmort
                nc%dmort            = currentCohort%dmort
                nc%lmort_direct     = currentCohort%lmort_direct
                nc%lmort_collateral = currentCohort%lmort_collateral
                nc%lmort_infra      = currentCohort%lmort_infra
                
             ! Logging is the dominant disturbance  
             elseif (currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifall) .and. &
                     currentPatch%disturbance_rates(dtype_ilog) > currentPatch%disturbance_rates(dtype_ifire)) then  ! Logging 

                ! If this cohort is in the upper canopy. It generated 
                if(currentCohort%canopy_layer == 1)then
                   
                   ! Trees generating this disturbance are not there by definition
                   nc%n            = 0.0_r8 

                   ! Reduce counts in the existing/donor patch according to the logging rate
                   currentCohort%n = currentCohort%n * (1.0_r8 - min(1.0_r8,(currentCohort%lmort_direct +    &
                                                                             currentCohort%lmort_collateral + &
                                                                             currentCohort%lmort_infra)))

                   ! The mortality diagnostics are set to nan because the cohort should dissappear
                   nc%cmort            = nan
                   nc%hmort            = nan
                   nc%bmort            = nan
                   nc%frmort           = nan
                   nc%fmort            = nan
                   nc%lmort_direct     = nan
                   nc%lmort_collateral = nan
                   nc%lmort_infra      = nan

                else

                   ! WHat to do with cohorts in the understory of a logging generated
                   ! disturbance patch?

                   if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then


                      ! Survivorship of undestory woody plants.  Two step process.
                      ! Step 1:  Reduce current number of plants to reflect the change in area.
                      !          The number density per square are doesn't change, but since the patch is smaller
                      !          and cohort counts are absolute, reduce this number.
                      nc%n = currentCohort%n * patch_site_areadis/currentPatch%area

                      ! because the mortality rate due to impact for the cohorts which had been in the understory and are now in the newly-
                      ! disturbed patch is very high, passing the imort directly to history results in large numerical errors, on account
                      ! of the sharply reduced number densities.  so instead pass this info via a site-level diagnostic variable before reducing 
                      ! the number density.
                      currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) = &
                           currentSite%imort_rate(currentCohort%size_class, currentCohort%pft) + &
                           nc%n * logging_coll_under_frac / hlm_freq_day
                      currentSite%imort_carbonflux = currentSite%imort_carbonflux + &
                           (nc%n * logging_coll_under_frac/ hlm_freq_day ) * &
                           currentCohort%b_total() * g_per_kg * days_per_sec * years_per_day * ha_per_m2

                      
                      ! Step 2:  Apply survivor ship function based on the understory death fraction
                     
                      ! remaining of understory plants of those that are knocked over by the overstorey trees dying...  
                      ! LOGGING SURVIVORSHIP OF UNDERSTORY PLANTS IS SET AS A NEW PARAMETER in the fatesparameter files 
                      nc%n = nc%n * (1.0_r8 - logging_coll_under_frac)

                      ! Step 3: Reduce the number count of cohorts in the original/donor/non-disturbed patch 
                      !         to reflect the area change
                      currentCohort%n = currentCohort%n * (1._r8 -  patch_site_areadis/currentPatch%area)


                      nc%fmort            = 0.0_r8
                      nc%cmort            = currentCohort%cmort
                      nc%hmort            = currentCohort%hmort
                      nc%bmort            = currentCohort%bmort
                      nc%frmort           = currentCohort%frmort
                      nc%dmort            = currentCohort%dmort
                      nc%lmort_direct     = currentCohort%lmort_direct
                      nc%lmort_collateral = currentCohort%lmort_collateral
                      nc%lmort_infra      = currentCohort%lmort_infra

                   else
                      
                      ! grass is not killed by mortality disturbance events. Just move it into the new patch area. 
                      ! Just split the grass into the existing and new patch structures
                      nc%n = currentCohort%n * patch_site_areadis/currentPatch%area
                      
                      ! Those remaining in the existing
                      currentCohort%n = currentCohort%n * (1._r8 - patch_site_areadis/currentPatch%area)

                      ! No grass impact mortality imposed on the newly created patch
                      nc%fmort            = 0.0_r8
                      nc%cmort            = currentCohort%cmort
                      nc%hmort            = currentCohort%hmort
                      nc%bmort            = currentCohort%bmort
                      nc%frmort           = currentCohort%frmort
                      nc%dmort            = currentCohort%dmort
                      nc%lmort_direct    = currentCohort%lmort_direct
                      nc%lmort_collateral = currentCohort%lmort_collateral
                      nc%lmort_infra      = currentCohort%lmort_infra
                      
                   endif  ! is/is-not woody
                   
                endif  ! Select canopy layer

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
                call insert_cohort(nc, new_patch%tallest, new_patch%shortest, tnull, snull, storebigcohort, storesmallcohort)

                new_patch%tallest  => storebigcohort 
                new_patch%shortest => storesmallcohort   
             else
                if(hlm_use_planthydro.eq.itrue) call DeallocateHydrCohort(nc)
                deallocate(nc) !get rid of the new memory.
             endif

             currentCohort => currentCohort%taller      
          enddo ! currentCohort 
          call sort_cohorts(currentPatch)

          !zero disturbance accumulators
          currentPatch%disturbance_rate  = 0._r8
          currentPatch%disturbance_rates = 0._r8

          !update area of donor patch
          currentPatch%area = currentPatch%area - patch_site_areadis

          ! sort out the cohorts, since some of them may be so small as to need removing. 
          ! the first call to terminate cohorts removes sparse number densities,
          ! the second call removes for all other reasons (sparse culling must happen
          ! before fusion)
          call terminate_cohorts(currentSite, currentPatch, 1)
          call fuse_cohorts(currentPatch, bc_in)
          call terminate_cohorts(currentSite, currentPatch, 2)
          call sort_cohorts(currentPatch)

          currentPatch => currentPatch%younger

       enddo ! currentPatch patch loop. 

          
       !*************************/
       !**  INSERT NEW PATCH INTO LINKED LIST    
       !**********`***************/        
       currentPatch               => currentSite%youngest_patch
       new_patch%older            => currentPatch
       new_patch%younger          => NULL()
       currentPatch%younger       => new_patch
       currentSite%youngest_patch => new_patch

       ! sort out the cohorts, since some of them may be so small as to need removing. 
       ! the first call to terminate cohorts removes sparse number densities,
       ! the second call removes for all other reasons (sparse culling must happen
       ! before fusion)
       call terminate_cohorts(currentSite, new_patch, 1)
       call fuse_cohorts(new_patch, bc_in)
       call terminate_cohorts(currentSite, new_patch, 2)
       call sort_cohorts(new_patch)

    endif !end new_patch area 

    call check_patch_area(currentSite)
    call set_patchno(currentSite)

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
    real(r8) :: areatot
    type(ed_patch_type), pointer :: currentPatch 
    !---------------------------------------------------------------------

    areatot = 0._r8
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       areatot = areatot + currentPatch%area
       currentPatch => currentPatch%younger
       if (( areatot - area ) > 0._r8 ) then 
          write(fates_log(),*) 'trimming patch area - is too big' , areatot-area
          currentSite%oldest_patch%area = currentSite%oldest_patch%area - (areatot - area)
       endif
    enddo

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
  subroutine average_patch_properties( currentPatch, newPatch, patch_site_areadis )
    !
    ! !DESCRIPTION:
    ! Average together the state properties of all of the donor patches that
    ! make up the new patch. 
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_patch_type) , intent(in), target  :: currentPatch
    type(ed_patch_type) , intent(inout)       :: newPatch 
    real(r8)            , intent(out)         :: patch_site_areadis   ! amount of land disturbed in this patch. m2
    !
    ! !LOCAL VARIABLES:
    integer  :: c,p ! counters for PFT and litter size class. 
    !---------------------------------------------------------------------

    patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate ! how much land is disturbed in this donor patch? 
 
    do c = 1,ncwd !move litter pool en mass into the new patch. 
       newPatch%cwd_ag(c) = newPatch%cwd_ag(c) + currentPatch%cwd_ag(c) * patch_site_areadis/newPatch%area
       newPatch%cwd_bg(c) = newPatch%cwd_bg(c) + currentPatch%cwd_bg(c) * patch_site_areadis/newPatch%area
    enddo

    do p = 1,numpft !move litter pool en mass into the new patch
       newPatch%root_litter(p) = newPatch%root_litter(p) + &
             currentPatch%root_litter(p) * patch_site_areadis/newPatch%area
       newPatch%leaf_litter(p) = newPatch%leaf_litter(p) + &
             currentPatch%leaf_litter(p) * patch_site_areadis/newPatch%area

       ! The fragmentation/decomposition flux from donor patches has already occured in existing patches.  However
       ! some of their area has been carved out for this new patches which is receiving donations.
       ! Lets maintain conservation on that pre-existing mass flux in these newly disturbed patches
       
       newPatch%root_litter_out(p) = newPatch%root_litter_out(p) + &
             currentPatch%root_litter_out(p) * patch_site_areadis/newPatch%area
       newPatch%leaf_litter_out(p) = newPatch%leaf_litter_out(p) + &
             currentPatch%leaf_litter_out(p) * patch_site_areadis/newPatch%area

    enddo

  end subroutine average_patch_properties

  ! ============================================================================
  subroutine fire_litter_fluxes(currentSite, cp_target, new_patch_target, patch_site_areadis)
    !
    ! !DESCRIPTION:
    !  CWD pool burned by a fire. 
    !  Carbon going from burned trees into CWD pool
    !  Burn parts of trees that don't die in fire
    !  Burn live grasses and kill them. 
    !
    ! !USES:
    use SFParamsMod,          only : SF_VAL_CWD_FRAC
    use EDtypesMod          , only : dl_sf
    !
    ! !ARGUMENTS:
    type(ed_site_type)  , intent(inout), target :: currentSite
    type(ed_patch_type) , intent(inout), target :: cp_target
    type(ed_patch_type) , intent(inout), target :: new_patch_target
    real(r8)            , intent(inout)         :: patch_site_areadis
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type) , pointer :: currentPatch
    type(ed_patch_type) , pointer :: new_patch
    type(ed_cohort_type), pointer :: currentCohort
    real(r8) :: bcroot               ! amount of below ground coarse root per cohort  kgC. (goes into CWD_BG)
    real(r8) :: bstem                ! amount of above ground stem biomass per cohort  kgC.(goes into CWG_AG)
    real(r8) :: dead_tree_density    ! no trees killed by fire per m2
    reaL(r8) :: burned_litter        ! amount of each litter pool burned by fire.  kgC/m2/day
    real(r8) :: burned_leaves        ! amount of tissue consumed by fire for grass. KgC/individual/day
    integer  :: c, p 
    !---------------------------------------------------------------------

    !check that total area is not exceeded. 
    currentPatch => cp_target
    new_patch => new_patch_target

    if ( currentPatch%fire  ==  1 ) then !only do this if there was a fire in this actual patch. 
       patch_site_areadis = currentPatch%area * currentPatch%disturbance_rate ! how much land is disturbed in this donor patch? 

       !************************************/ 
       !PART 1)  Burn the fractions of existing litter in the new patch that were consumed by the fire. 
       !************************************/ 
       do c = 1,ncwd
          burned_litter = new_patch%cwd_ag(c) * patch_site_areadis/new_patch%area * &
                currentPatch%burnt_frac_litter(c+1) !kG/m2/day
          new_patch%cwd_ag(c) = new_patch%cwd_ag(c) - burned_litter
          currentSite%flux_out = currentSite%flux_out + burned_litter * new_patch%area !kG/site/day
          currentSite%total_burn_flux_to_atm = currentSite%total_burn_flux_to_atm + &
                burned_litter * new_patch%area !kG/site/day
       enddo

       do p = 1,numpft
          burned_litter = new_patch%leaf_litter(p) * patch_site_areadis/new_patch%area * &
                currentPatch%burnt_frac_litter(dl_sf)
          new_patch%leaf_litter(p) = new_patch%leaf_litter(p) - burned_litter
          currentSite%flux_out = currentSite%flux_out + burned_litter * new_patch%area !kG/site/day
          currentSite%total_burn_flux_to_atm = currentSite%total_burn_flux_to_atm + &
                burned_litter * new_patch%area !kG/site/day
      enddo

       !************************************/     
       !PART 2) Put unburned parts of plants that died in the fire into the litter pool of new and old patches 
       ! This happens BEFORE the plant numbers have been updated. So we are working with the 
       ! pre-fire population of plants, which is the right way round. 
       !************************************/ 
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))
          p = currentCohort%pft
          if(EDPftvarcon_inst%woody(p) == 1)then !DEAD (FROM FIRE) TREES
             !************************************/ 
             ! Number of trees that died because of the fire, per m2 of ground. 
             ! Divide their litter into the four litter streams, and spread evenly across ground surface. 
             !************************************/  
             ! stem biomass per tree
             bstem  = (currentCohort%bsw + currentCohort%bdead) * EDPftvarcon_inst%allom_agb_frac(p)
             ! coarse root biomass per tree
             bcroot = (currentCohort%bsw + currentCohort%bdead) * (1.0_r8 - EDPftvarcon_inst%allom_agb_frac(p) )
             ! density of dead trees per m2. 
             dead_tree_density  = (currentCohort%fire_mort * currentCohort%n*patch_site_areadis/currentPatch%area) / AREA  

             ! Unburned parts of dead tree pool. 
             ! Unburned leaves and roots    
             
             new_patch%leaf_litter(p) = new_patch%leaf_litter(p) + dead_tree_density * (currentCohort%bl) &
             * (1.0_r8-currentCohort%cfa)
             new_patch%root_litter(p) = new_patch%root_litter(p) + dead_tree_density * (currentCohort%br+currentCohort%bstore)
             currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + dead_tree_density * &
                  (currentCohort%bl) * (1.0_r8-currentCohort%cfa)
             currentPatch%root_litter(p) = currentPatch%root_litter(p) + dead_tree_density * &
                  (currentCohort%br+currentCohort%bstore)

             ! track as diagnostic fluxes
             currentSite%leaf_litter_diagnostic_input_carbonflux(p) = currentSite%leaf_litter_diagnostic_input_carbonflux(p) + &
                  (currentCohort%bl) * (1.0_r8-currentCohort%cfa) * currentCohort%fire_mort * currentCohort%n * &
                  hlm_days_per_year / AREA
             currentSite%root_litter_diagnostic_input_carbonflux(p) = currentSite%root_litter_diagnostic_input_carbonflux(p) + &
                  (currentCohort%br+currentCohort%bstore) * (1.0_r8-currentCohort%cfa) * currentCohort%fire_mort * &
                  currentCohort%n * hlm_days_per_year / AREA
      
             ! below ground coarse woody debris from burned trees
             do c = 1,ncwd
                new_patch%cwd_bg(c) = new_patch%cwd_bg(c) + dead_tree_density * SF_val_CWD_frac(c) * bcroot
                currentPatch%cwd_bg(c) = currentPatch%cwd_bg(c) + dead_tree_density * SF_val_CWD_frac(c) * bcroot

                ! track as diagnostic fluxes
                currentSite%CWD_BG_diagnostic_input_carbonflux(c) = currentSite%CWD_BG_diagnostic_input_carbonflux(c) + &
                     SF_val_CWD_frac(c) * bcroot * currentCohort%fire_mort * currentCohort%n * &
                     hlm_days_per_year / AREA
             enddo

             ! above ground coarse woody debris from unburned twigs and small branches
             do c = 1,2
                new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem &
                * (1.0_r8-currentCohort%cfa)
                currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * &
                     bstem * (1.0_r8-currentCohort%cfa)

                ! track as diagnostic fluxes
                currentSite%CWD_AG_diagnostic_input_carbonflux(c) = currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
                     SF_val_CWD_frac(c) * bstem * (1.0_r8-currentCohort%cfa) * currentCohort%fire_mort * currentCohort%n * &
                     hlm_days_per_year / AREA
             enddo
             
             ! above ground coarse woody debris from large branches and stems: these do not burn in crown fires. 
             do c = 3,4
                new_patch%cwd_ag(c) = new_patch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem
                currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + dead_tree_density * SF_val_CWD_frac(c) * bstem

                ! track as diagnostic fluxes
                currentSite%CWD_AG_diagnostic_input_carbonflux(c) = currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
                     SF_val_CWD_frac(c) * bstem * currentCohort%fire_mort * currentCohort%n * &
                     hlm_days_per_year / AREA
             enddo
             
             ! Burned parts of dead tree pool.  
             ! Burned twigs and small branches. 
             do c = 1,2

                currentSite%cwd_ag_burned(c) = currentSite%cwd_ag_burned(c) + dead_tree_density * &
                     SF_val_CWD_frac(c) * bstem * currentCohort%cfa
                currentSite%flux_out  = currentSite%flux_out + dead_tree_density * &
                     AREA * SF_val_CWD_frac(c) * bstem * currentCohort%cfa
                currentSite%total_burn_flux_to_atm  = currentSite%total_burn_flux_to_atm + dead_tree_density * &
                     AREA * SF_val_CWD_frac(c) * bstem * currentCohort%cfa

             enddo
             
             !burned leaves. 
             do p = 1,numpft                  

                currentSite%leaf_litter_burned(p) = currentSite%leaf_litter_burned(p) + &
                     dead_tree_density * currentCohort%bl * currentCohort%cfa
                currentSite%flux_out  = currentSite%flux_out + &
                     dead_tree_density * AREA * currentCohort%bl * currentCohort%cfa
                currentSite%total_burn_flux_to_atm  = currentSite%total_burn_flux_to_atm + &
                     dead_tree_density * AREA * currentCohort%bl * currentCohort%cfa

             enddo

         endif

          currentCohort => currentCohort%taller

       enddo  ! currentCohort

       !************************************/     
       ! PART 3) Burn parts of trees that did *not* die in the fire.
       ! PART 4) Burn parts of grass that are consumed by the fire. 
       ! grasses are not killed directly by fire. They die by losing all of their leaves and starving. 
       !************************************/ 
       currentCohort => new_patch%shortest
       do while(associated(currentCohort))

          call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread,currentCohort%pft,currentCohort%c_area)
          if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then
             burned_leaves = min(currentCohort%bl, (currentCohort%bl+currentCohort%bsw) * currentCohort%cfa)
          else
             burned_leaves = min(currentCohort%bl, (currentCohort%bl+currentCohort%bsw) * currentPatch%burnt_frac_litter(6))
          endif
          if (burned_leaves > 0.0_r8) then

             currentCohort%bl     = currentCohort%bl - burned_leaves

             !KgC/gridcell/day
             currentSite%flux_out = currentSite%flux_out + burned_leaves * currentCohort%n * &
                  patch_site_areadis/currentPatch%area * AREA 
             currentSite%total_burn_flux_to_atm = currentSite%total_burn_flux_to_atm+ burned_leaves * currentCohort%n * &
                  patch_site_areadis/currentPatch%area * AREA 

          endif
          currentCohort%cfa = 0.0_r8        

          currentCohort => currentCohort%taller

       enddo

    endif !currentPatch%fire. 

  end subroutine fire_litter_fluxes

  ! ============================================================================
  subroutine mortality_litter_fluxes(currentSite, cp_target, new_patch_target, patch_site_areadis)
    !
    ! !DESCRIPTION:
    !  Carbon going from ongoing mortality into CWD pools. 
    !
    ! !USES:
    use EDParamsMod,  only : ED_val_understorey_death
    use SFParamsMod,  only : SF_val_cwd_frac
    !
    ! !ARGUMENTS:
    type(ed_site_type)  , intent(inout), target :: currentSite 
    type(ed_patch_type) , intent(inout), target :: cp_target 
    type(ed_patch_type) , intent(inout), target :: new_patch_target
    real(r8)            , intent(in)            :: patch_site_areadis
    !
    ! !LOCAL VARIABLES:
    real(r8) :: cwd_litter_density
    real(r8) :: litter_area ! area over which to distribute this litter. 
    type(ed_cohort_type), pointer :: currentCohort
    type(ed_patch_type) , pointer :: currentPatch 
    type(ed_patch_type) , pointer :: new_patch 
    real(r8) :: understorey_dead  !Number of individual dead from the canopy layer /day
    real(r8) :: canopy_dead       !Number of individual dead from the understorey layer /day
    real(r8) :: np_mult           !Fraction of the new patch which came from the current patch (and so needs the same litter) 
    integer :: p,c
    real(r8) :: canopy_mortality_woody_litter(maxpft)    ! flux of wood litter in to litter pool: KgC/m2/day
    real(r8) :: canopy_mortality_leaf_litter(maxpft)     ! flux in to  leaf litter from tree death: KgC/m2/day
    real(r8) :: canopy_mortality_root_litter(maxpft)     ! flux in to froot litter  from tree death: KgC/m2/day
    !---------------------------------------------------------------------

    currentPatch => cp_target
    new_patch => new_patch_target
    canopy_mortality_woody_litter(:) = 0.0_r8 ! mortality generated litter. KgC/m2/day
    canopy_mortality_leaf_litter(:)  = 0.0_r8
    canopy_mortality_root_litter(:)  = 0.0_r8

    currentCohort => currentPatch%shortest
    do while(associated(currentCohort))       
       p = currentCohort%pft

          if(currentCohort%canopy_layer == 1)then         
             !currentCohort%dmort = mortality_rates(currentCohort) 
             !the disturbance calculations are done with the previous n, c_area and d_mort. So it's probably &
             !not right to recalcualte dmort here.
             canopy_dead = currentCohort%n * min(1.0_r8,currentCohort%dmort * hlm_freq_day * fates_mortality_disturbance_fraction)

             canopy_mortality_woody_litter(p)= canopy_mortality_woody_litter(p) + &
                  canopy_dead*(currentCohort%bdead+currentCohort%bsw)
             canopy_mortality_leaf_litter(p) = canopy_mortality_leaf_litter(p) + &
                  canopy_dead*(currentCohort%bl)
             canopy_mortality_root_litter(p) = canopy_mortality_root_litter(p) + &
                  canopy_dead*(currentCohort%br+currentCohort%bstore)

         else 
             if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then

                understorey_dead = ED_val_understorey_death * currentCohort%n * (patch_site_areadis/currentPatch%area)  !kgC/site/day
                canopy_mortality_woody_litter(p) = canopy_mortality_woody_litter(p)  + &
                     understorey_dead*(currentCohort%bdead+currentCohort%bsw)  
                canopy_mortality_leaf_litter(p)= canopy_mortality_leaf_litter(p)+ &
                     understorey_dead* currentCohort%bl 
                canopy_mortality_root_litter(p)= canopy_mortality_root_litter(p)+ &
                      understorey_dead*(currentCohort%br+currentCohort%bstore)

             ! FIX(SPM,040114) - clarify this comment
             ! grass is not killed by canopy mortality disturbance events.
             ! Just move it into the new patch area. 
             else 
                ! no-op
             endif
          endif
       

       currentCohort => currentCohort%taller      

    enddo !currentCohort         

    !************************************/
    !Evenly distribute the litter from the trees that died across the new and old patches
    !************************************/
    !************************************/
    !Evenly distribute the litter from the trees that died across the new and old patches
    !'litter' fluxes here are in KgC
    !************************************/
    litter_area = currentPatch%area 
    np_mult =  patch_site_areadis/new_patch%area
    ! This litter is distributed between the current and new patches, &
    ! not to any other patches. This is really the eventually area of the current patch &
    ! (currentPatch%area-patch_site_areadis) +patch_site_areadis...
    ! For the new patch, only some fraction of its land area (patch_areadis/np%area) is derived from the current patch
    ! so we need to multiply by patch_areadis/np%area

    do p = 1,numpft
    
       do c = 1,ncwd
          
          cwd_litter_density = SF_val_CWD_frac(c) * canopy_mortality_woody_litter(p) / litter_area
          
          new_patch%cwd_ag(c)    = new_patch%cwd_ag(c)    + EDPftvarcon_inst%allom_agb_frac(p) * cwd_litter_density * np_mult
          currentPatch%cwd_ag(c) = currentPatch%cwd_ag(c) + EDPftvarcon_inst%allom_agb_frac(p) * cwd_litter_density
          new_patch%cwd_bg(c)    = new_patch%cwd_bg(c)    + (1._r8-EDPftvarcon_inst%allom_agb_frac(p)) * cwd_litter_density * np_mult 
          currentPatch%cwd_bg(c) = currentPatch%cwd_bg(c) + (1._r8-EDPftvarcon_inst%allom_agb_frac(p)) * cwd_litter_density 
          
          ! track as diagnostic fluxes
          currentSite%CWD_AG_diagnostic_input_carbonflux(c) = currentSite%CWD_AG_diagnostic_input_carbonflux(c) + &
                SF_val_CWD_frac(c) * canopy_mortality_woody_litter(p) * hlm_days_per_year * EDPftvarcon_inst%allom_agb_frac(p)/ AREA 
          currentSite%CWD_BG_diagnostic_input_carbonflux(c) = currentSite%CWD_BG_diagnostic_input_carbonflux(c) + &
                SF_val_CWD_frac(c) * canopy_mortality_woody_litter(p) * hlm_days_per_year * (1.0_r8 - EDPftvarcon_inst%allom_agb_frac(p)) / AREA
       enddo

       new_patch%leaf_litter(p) = new_patch%leaf_litter(p) + canopy_mortality_leaf_litter(p) / litter_area * np_mult
       new_patch%root_litter(p) = new_patch%root_litter(p) + canopy_mortality_root_litter(p) / litter_area * np_mult 

       currentPatch%leaf_litter(p) = currentPatch%leaf_litter(p) + canopy_mortality_leaf_litter(p) / litter_area
       currentPatch%root_litter(p) = currentPatch%root_litter(p) + canopy_mortality_root_litter(p) / litter_area

       ! track as diagnostic fluxes
       currentSite%leaf_litter_diagnostic_input_carbonflux(p) = currentSite%leaf_litter_diagnostic_input_carbonflux(p) + &
            canopy_mortality_leaf_litter(p) * hlm_days_per_year / AREA

       currentSite%root_litter_diagnostic_input_carbonflux(p) = currentSite%root_litter_diagnostic_input_carbonflux(p) + &
            canopy_mortality_root_litter(p) * hlm_days_per_year / AREA
    enddo

  end subroutine mortality_litter_fluxes

  ! ============================================================================
  subroutine create_patch(currentSite, new_patch, age, areap,cwd_ag_local,cwd_bg_local, &
       leaf_litter_local,root_litter_local)
    !
    ! !DESCRIPTION:
    !  Set default values for creating a new patch
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout), target :: currentSite
    type(ed_patch_type), intent(inout), target :: new_patch
    real(r8), intent(in) :: age                 ! notional age of this patch in years
    real(r8), intent(in) :: areap               ! initial area of this patch in m2. 
    real(r8), intent(in) :: cwd_ag_local(:)     ! initial value of above ground coarse woody debris. KgC/m2
    real(r8), intent(in) :: cwd_bg_local(:)     ! initial value of below ground coarse woody debris. KgC/m2
    real(r8), intent(in) :: root_litter_local(:)! initial value of root litter. KgC/m2
    real(r8), intent(in) :: leaf_litter_local(:)! initial value of leaf litter. KgC/m2
    !
    ! !LOCAL VARIABLES:
    !---------------------------------------------------------------------

    allocate(new_patch%tr_soil_dir(hlm_numSWb))
    allocate(new_patch%tr_soil_dif(hlm_numSWb))
    allocate(new_patch%tr_soil_dir_dif(hlm_numSWb))
    allocate(new_patch%fab(hlm_numSWb))
    allocate(new_patch%fabd(hlm_numSWb))
    allocate(new_patch%fabi(hlm_numSWb))
    allocate(new_patch%sabs_dir(hlm_numSWb))
    allocate(new_patch%sabs_dif(hlm_numSWb))
    allocate(new_patch%rootfr_ft(numpft,hlm_numlevgrnd))
    allocate(new_patch%rootr_ft(numpft,hlm_numlevgrnd)) 
    
    call zero_patch(new_patch) !The nan value in here is not working??

    new_patch%tallest  => null() ! pointer to patch's tallest cohort    
    new_patch%shortest => null() ! pointer to patch's shortest cohort   
    new_patch%older    => null() ! pointer to next older patch   
    new_patch%younger  => null() ! pointer to next shorter patch      

    ! assign known patch attributes 

    new_patch%age                = age   
    new_patch%age_class          = 1
    new_patch%area               = areap 
    new_patch%cwd_ag             = cwd_ag_local
    new_patch%cwd_bg             = cwd_bg_local
    new_patch%leaf_litter        = leaf_litter_local
    new_patch%root_litter        = root_litter_local
 
    !zeroing things because of the surfacealbedo problem... shouldnt really be necesary
    new_patch%cwd_ag_in(:)       = 0._r8
    new_patch%cwd_bg_in(:)       = 0._r8

    new_patch%cwd_ag_out(:)      = 0._r8
    new_patch%cwd_bg_out(:)      = 0._r8

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
    new_patch%frac_burnt         = 0._r8  
    new_patch%total_tree_area    = 0.0_r8  
    new_patch%NCL_p              = 1
 
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
    currentPatch%canopy_layer_tai(:)        = nan               
    currentPatch%total_canopy_area          = nan
    currentPatch%bare_frac_area             = nan                             

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

    ! LITTER
    currentPatch%cwd_ag(:)                  = 0.0_r8 ! above ground coarse woody debris gc/m2. 
    currentPatch%cwd_bg(:)                  = 0.0_r8 ! below ground coarse woody debris
    currentPatch%root_litter(:)             = 0.0_r8 ! In new disturbed patches, loops over donors to increment total, needs zero here
    currentPatch%leaf_litter(:)             = 0.0_r8 ! In new disturbed patches, loops over donors to increment total, needs zero here

    ! Cold-start initialized patches should have no litter flux in/out as they have not undergone any time.
    ! Litter fluxes in/out also need to be initialized to zero for newly disturbed patches, as they
    ! will incorporate the fluxes from donors over a loop, and need an initialization

    currentPatch%leaf_litter_in(:)  = 0.0_r8 ! As a newly created patch with no age, there is no flux in
    currentPatch%leaf_litter_out(:) = 0.0_r8 ! As a newly created patch with no age, no frag or decomp has happened yet
    currentPatch%root_litter_in(:)  = 0.0_r8 ! As a newly created patch with no age, there is no flux in
    currentPatch%root_litter_out(:) = 0.0_r8 ! As a newly created patch with no age, no frag or decomp has happened yet

    ! FIRE
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
    currentPatch%ab                         = 0.0_r8 ! area burnt daily m2
    currentPatch%nf                         = 0.0_r8 ! number of fires initiated daily 
    currentPatch%sh                         = 0.0_r8 ! average scorch height for the patch(m)
    currentPatch%frac_burnt                 = 0.0_r8 ! fraction burnt in each timestep. 
    currentPatch%burnt_frac_litter(:)       = 0.0_r8 
    currentPatch%btran_ft(:)                = 0.0_r8

    currentPatch%canopy_layer_tai(:)        = 0.0_r8

    currentPatch%seeds_in(:)                = 0.0_r8
    currentPatch%seed_decay(:)              = 0.0_r8
    currentPatch%seed_germination(:)        = 0.0_r8

    currentPatch%fab(:)                     = 0.0_r8
    currentPatch%sabs_dir(:)                = 0.0_r8
    currentPatch%sabs_dif(:)                = 0.0_r8
    currentPatch%zstar                      = 0.0_r8
    currentPatch%c_stomata                  = 0.0_r8 ! This is calculated immediately before use
    currentPatch%c_lblayer                  = 0.0_r8

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
    integer  :: nopatches   !number of patches presently in gridcell
    integer  :: iterate     !switch of patch reduction iteration scheme. 1 to keep going, 0 to stop
    integer  :: fuse_flag   !do patches get fused (1) or not (0). 
    !
    !---------------------------------------------------------------------

    currentSite => csite 

    profiletol = ED_val_patch_fusion_tol

    nopatches = 0
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       nopatches = nopatches +1
       currentPatch => currentPatch%older
    enddo
    !---------------------------------------------------------------------!
    !  We only really care about fusing patches if nopatches > 1           !
    !---------------------------------------------------------------------!
    iterate = 1

    !---------------------------------------------------------------------!
    !  Keep doing this until nopatches >= maxPatchesPerSite                         !
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

       !---------------------------------------------------------------------!
       ! Loop round current & target (currentPatch,tpp) patches to assess combinations !
       !---------------------------------------------------------------------!   
       currentPatch => currentSite%youngest_patch
       do while(associated(currentPatch))      
          tpp => currentSite%youngest_patch
          do while(associated(tpp))

             if(.not.associated(currentPatch))then
                write(fates_log(),*) 'ED: issue with currentPatch'
             endif

             if(associated(tpp).and.associated(currentPatch))then

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

                      
                      !---------------------------------------------------------------------------------------------------------
                      ! the next bit of logic forces fusion of two patches which both have tiny biomass densities. without this,
                      ! fates gives a bunch of really young patches which all have almost no biomass and so don't need to be 
                      ! distinguished from each other. but if force_patchfuse_min_biomass is too big, it takes too long for the 
                      ! youngest patch to build up enough biomass to be its own distinct entity, which leads to large oscillations 
                      ! in the patch dynamics and dependent variables.
                      !---------------------------------------------------------------------------------------------------------
                      
                      if(sum(currentPatch%pft_agb_profile(:,:)) > force_patchfuse_min_biomass .or. &
                           sum(tpp%pft_agb_profile(:,:)) > force_patchfuse_min_biomass ) then

                         !---------------------------------------------------------------------!
                         ! Calculate the difference criteria for each pft and dbh class        !
                         !---------------------------------------------------------------------!   

                         do ft = 1,numpft        ! loop over pfts
                            do z = 1,n_dbh_bins      ! loop over hgt bins 

                               !----------------------------------
                               !is there biomass in this category?
                               !----------------------------------

                               if(currentPatch%pft_agb_profile(ft,z)  > 0.0_r8 .or.  &
                                    tpp%pft_agb_profile(ft,z) > 0.0_r8)then 

                                  !-------------------------------------------------------------------------------------
                                  ! what is the relative difference in biomass i nthis category between the two patches?
                                  !-------------------------------------------------------------------------------------

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
                      tmpptr => currentPatch%older       
                      call fuse_2_patches(csite, currentPatch, tpp)
                      call fuse_cohorts(tpp, bc_in)
                      call sort_cohorts(tpp)
                      currentPatch => tmpptr
                   else
                     ! write(fates_log(),*) 'patches not fused'
                   endif
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
       nopatches = 0
       currentPatch => currentSite%youngest_patch
       do while(associated(currentPatch))
          nopatches = nopatches +1
          currentPatch => currentPatch%older
       enddo

       if(nopatches > maxPatchesPerSite)then
          iterate = 1
          profiletol = profiletol * patch_fusion_tolerance_relaxation_increment

          !---------------------------------------------------------------------!
          ! Making profile tolerance larger means that more fusion will happen  !
          !---------------------------------------------------------------------!        
       else
          iterate = 0
       endif

    enddo !do while nopatches>maxPatchesPerSite
 
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
    type(ed_patch_type), pointer   :: youngerp     ! pointer to the patch younger than donor
    type(ed_patch_type), pointer   :: olderp       ! pointer to the patch older than donor
    real(r8)                       :: inv_sum_area ! Inverse of the sum of the two patches areas
    !-----------------------------------------------------------------------------------------------

    ! Generate a litany of area weighted averages

    inv_sum_area = 1.0_r8/(dp%area + rp%area)
    
    rp%age = (dp%age * dp%area + rp%age * rp%area) * inv_sum_area

    rp%age_class = get_age_class_index(rp%age)


    do c = 1,ncwd
       rp%cwd_ag(c) = (dp%cwd_ag(c)*dp%area + rp%cwd_ag(c)*rp%area) * inv_sum_area
       rp%cwd_bg(c) = (dp%cwd_bg(c)*dp%area + rp%cwd_bg(c)*rp%area) * inv_sum_area
    enddo
    
    do p = 1,numpft 
       rp%seeds_in(p)         = (rp%seeds_in(p)*rp%area + dp%seeds_in(p)*dp%area) * inv_sum_area
       rp%seed_decay(p)       = (rp%seed_decay(p)*rp%area + dp%seed_decay(p)*dp%area) * inv_sum_area
       rp%seed_germination(p) = (rp%seed_germination(p)*rp%area + dp%seed_germination(p)*dp%area) * inv_sum_area

       rp%leaf_litter(p)      = (dp%leaf_litter(p)*dp%area + rp%leaf_litter(p)*rp%area) * inv_sum_area
       rp%root_litter(p)      = (dp%root_litter(p)*dp%area + rp%root_litter(p)*rp%area) * inv_sum_area

       rp%root_litter_out(p)  = (dp%root_litter_out(p)*dp%area + rp%root_litter_out(p)*rp%area) * inv_sum_area
       rp%leaf_litter_out(p)  = (dp%leaf_litter_out(p)*dp%area + rp%leaf_litter_out(p)*rp%area) * inv_sum_area

       rp%root_litter_in(p)   = (dp%root_litter_in(p)*dp%area + rp%root_litter_in(p)*rp%area) * inv_sum_area
       rp%leaf_litter_in(p)   = (dp%leaf_litter_in(p)*dp%area + rp%leaf_litter_in(p)*rp%area) * inv_sum_area

       rp%dleaf_litter_dt(p)  = (dp%dleaf_litter_dt(p)*dp%area + rp%dleaf_litter_dt(p)*rp%area) * inv_sum_area
       rp%droot_litter_dt(p)  = (dp%droot_litter_dt(p)*dp%area + rp%droot_litter_dt(p)*rp%area) * inv_sum_area
    enddo
    
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
    rp%ab                   = (dp%ab*dp%area + rp%ab*rp%area) * inv_sum_area
    rp%nf                   = (dp%nf*dp%area + rp%nf*rp%area) * inv_sum_area
    rp%sh                   = (dp%sh*dp%area + rp%sh*rp%area) * inv_sum_area
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

    real(r8) areatot ! variable for checking whether the total patch area is wrong. 
    !---------------------------------------------------------------------
 
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch)) 
       
       if(currentPatch%area <= min_patch_area)then
          
          if ( .not.associated(currentPatch,currentSite%youngest_patch) ) then
             
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
                
             elseif( .not. associated(currentPatch%younger,currentSite%youngest_patch) ) then
                
                if(debug) &
                      write(fates_log(),*) 'fusing to younger patch because oldest one is too small', &
                      currentPatch%area
                
                youngerPatch => currentPatch%younger
                call fuse_2_patches(currentSite, youngerPatch, currentPatch)
                
                ! The fusion process has updated the "younger" pointer on currentPatch
                
             endif
          endif
       endif
          
       currentPatch => currentPatch%older
       
    enddo
    
    !check area is not exceeded
    areatot = 0._r8
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))
       areatot = areatot + currentPatch%area
       currentPatch => currentPatch%younger
       if((areatot-area) > 0.0000001_r8)then
          write(fates_log(),*) 'ED: areatot too large. end terminate', areatot
       endif
    enddo

  end subroutine terminate_patches

  ! =====================================================================================

  subroutine dealloc_patch(cpatch)

    ! This Subroutine is intended to de-allocate the allocatable memory that is pointed
    ! to via the patch structure.  This subroutine DOES NOT deallocate the patch
    ! structure itself.

    type(ed_patch_type), target :: cpatch
    type(ed_cohort_type), pointer :: ccohort  ! current
    type(ed_cohort_type), pointer :: ncohort  ! next
    
    ! First Deallocate the cohort space
    ! -----------------------------------------------------------------------------------
    ccohort => cpatch%shortest
    do while(associated(ccohort))
       
       ncohort => ccohort%taller
       if(hlm_use_planthydro.eq.itrue) call DeallocateHydrCohort(ccohort)
       deallocate(ccohort)
       ccohort => ncohort

    end do

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
       deallocate(cpatch%rootfr_ft)
       deallocate(cpatch%rootr_ft)
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
                  currentCohort%bdead*currentCohort%n/currentPatch%area

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

   ! ====================================================================================

  subroutine set_root_fraction( cpatch , zi )
    !
    ! !DESCRIPTION:
    !  Calculates the fractions of the root biomass in each layer for each pft. 
    !
    ! !USES:

    !
    ! !ARGUMENTS
    type(ed_patch_type),intent(inout), target :: cpatch
    real(r8),intent(in)  :: zi(0:hlm_numlevsoil)
    !
    ! !LOCAL VARIABLES:
    integer :: lev,p,c,ft
    !----------------------------------------------------------------------
    
    do ft = 1,numpft
       do lev = 1, hlm_numlevgrnd
          cpatch%rootfr_ft(ft,lev) = 0._r8
       enddo

       do lev = 1, hlm_numlevsoil-1
          cpatch%rootfr_ft(ft,lev) = .5_r8*( &
                 exp(-EDPftvarcon_inst%roota_par(ft) * zi(lev-1))  &
               + exp(-EDPftvarcon_inst%rootb_par(ft) * zi(lev-1))  &
               - exp(-EDPftvarcon_inst%roota_par(ft) * zi(lev))    &
               - exp(-EDPftvarcon_inst%rootb_par(ft) * zi(lev)))
       end do
    end do

  end subroutine set_root_fraction

 end module EDPatchDynamicsMod
