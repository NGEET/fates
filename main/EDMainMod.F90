module EDMainMod

  ! ===========================================================================
  ! Main ED module.    
  ! ============================================================================

  use shr_kind_mod             , only : r8 => shr_kind_r8
  
  use FatesGlobals             , only : fates_log
  use FatesInterfaceMod        , only : hlm_freq_day
  use FatesInterfaceMod        , only : hlm_day_of_year
  use FatesInterfaceMod        , only : hlm_days_per_year
  use FatesInterfaceMod        , only : hlm_current_year
  use FatesInterfaceMod        , only : hlm_current_month
  use FatesInterfaceMod        , only : hlm_current_day 
  use FatesInterfaceMod        , only : hlm_use_planthydro 
  use FatesInterfaceMod        , only : hlm_reference_date
  use FatesInterfaceMod        , only : hlm_use_ed_prescribed_phys
  use FatesInterfaceMod        , only : hlm_use_ed_st3 
  use FatesInterfaceMod        , only : bc_in_type
  use FatesInterfaceMod        , only : hlm_masterproc
  use FatesInterfaceMod        , only : numpft
  use EDCohortDynamicsMod      , only : terminate_cohorts
  use EDCohortDynamicsMod      , only : fuse_cohorts
  use EDCohortDynamicsMod      , only : sort_cohorts
  use EDCohortDynamicsMod      , only : count_cohorts
  use EDCohortDynamicsMod      , only : EvaluateAndCorrectDBH
  use EDPatchDynamicsMod       , only : disturbance_rates
  use EDPatchDynamicsMod       , only : fuse_patches
  use EDPatchDynamicsMod       , only : spawn_patches
  use EDPatchDynamicsMod       , only : terminate_patches
  use EDPhysiologyMod          , only : non_canopy_derivs
  use EDPhysiologyMod          , only : phenology
  use EDPhysiologyMod          , only : recruitment
  use EDPhysiologyMod          , only : trim_canopy
  use EDPhysiologyMod          , only : ZeroAllocationRates
  use EDCohortDynamicsMod      , only : UpdateCohortBioPhysRates
  use SFMainMod                , only : fire_model 
  use FatesSizeAgeTypeIndicesMod, only : get_age_class_index
  use EDtypesMod               , only : ncwd
  use EDtypesMod               , only : ed_site_type
  use EDtypesMod               , only : ed_patch_type
  use EDtypesMod               , only : ed_cohort_type
  use EDTypesMod               , only : AREA
  use EDTypesMod               , only : phen_dstat_moiston
  use EDTypesMod               , only : phen_dstat_timeon
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesConstantsMod        , only : primaryforest, secondaryforest
  use FatesPlantHydraulicsMod  , only : do_growthrecruiteffects
  use FatesPlantHydraulicsMod  , only : updateSizeDepTreeHydProps
  use FatesPlantHydraulicsMod  , only : updateSizeDepTreeHydStates
  use FatesPlantHydraulicsMod  , only : initTreeHydStates
  use FatesPlantHydraulicsMod  , only : updateSizeDepRhizHydProps 
  use FatesAllometryMod        , only : h_allom,tree_sai,tree_lai
  use FatesPlantHydraulicsMod , only : updateSizeDepRhizHydStates
  use EDLoggingMortalityMod    , only : IsItLoggingTime
  use FatesGlobals             , only : endrun => fates_endrun
  use ChecksBalancesMod        , only : SiteCarbonStock
  use EDMortalityFunctionsMod  , only : Mortality_Derivative

  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : all_carbon_elements
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ

  use PRTLossFluxesMod,       only : PRTMaintTurnover
  use PRTLossFluxesMod,       only : PRTReproRelease

  use EDPftvarcon,            only : EDPftvarcon_inst

  ! CIME Globals
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_infnan_mod      , only : nan => shr_infnan_nan, assignment(=)

  implicit none
  private

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: ed_ecosystem_dynamics
  public  :: ed_update_site
  !
  ! !PRIVATE MEMBER FUNCTIONS:
  
  private :: ed_integrate_state_variables
  private :: ed_total_balance_check
  private :: bypass_dynamics
  
  logical :: debug  = .false.
  
  character(len=*), parameter, private :: sourcefile = &
         __FILE__
  
  !
  ! 10/30/09: Created by Rosie Fisher
  !-----------------------------------------------------------------------

contains

  !-------------------------------------------------------------------------------!
  subroutine ed_ecosystem_dynamics(currentSite, bc_in)
    !
    ! !DESCRIPTION:
    !  Core of ed model, calling all subsequent vegetation dynamics routines         
    !
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), target  :: currentSite
    type(bc_in_type)        , intent(in)             :: bc_in
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type), pointer :: currentPatch
    !-----------------------------------------------------------------------

    if ( hlm_masterproc==itrue ) write(fates_log(),'(A,I4,A,I2.2,A,I2.2)') 'FATES Dynamics: ',&
          hlm_current_year,'-',hlm_current_month,'-',hlm_current_day

    currentSite%flux_in = 0.0_r8


    ! Call a routine that simply identifies if logging should occur
    ! This is limited to a global event until more structured event handling is enabled
    call IsItLoggingTime(hlm_masterproc,currentSite)

    !**************************************************************************
    ! Fire, growth, biogeochemistry. 
    !**************************************************************************
    
    !FIX(SPM,032414) take this out.  On startup these values are all zero and on restart it
    !zeros out values read in the restart file

    ! Zero turnover rates and growth diagnostics
    call ZeroAllocationRates(currentSite)

   
    call ed_total_balance_check(currentSite, 0)

    ! We do not allow phenology while in ST3 mode either, it is hypothetically
    ! possible to allow this, but we have not plugged in the litter fluxes
    ! of flushing or turning over leaves for non-dynamics runs
    if (hlm_use_ed_st3.eq.ifalse) then
       call phenology(currentSite, bc_in )
    end if


    if (hlm_use_ed_st3.eq.ifalse) then   ! Bypass if ST3
       call fire_model(currentSite, bc_in) 

       ! Calculate disturbance and mortality based on previous timestep vegetation.
       ! disturbance_rates calls logging mortality and other mortalities, Yi Xu
       call disturbance_rates(currentSite, bc_in)
    end if

    if (hlm_use_ed_st3.eq.ifalse) then
       ! Integrate state variables from annual rates to daily timestep
       call ed_integrate_state_variables(currentSite, bc_in ) 
    else
       ! ed_intergrate_state_variables is where the new cohort flag
       ! is set. This flag designates wether a cohort has
       ! experienced a day, and therefore has been populated with non-nonsense
       ! values.  If we aren't entering that sequence, we need to set the flag
       ! Make sure cohorts are marked as non-recruits

       call bypass_dynamics(currentSite)
       
    end if

    !******************************************************************************
    ! Reproduction, Recruitment and Cohort Dynamics : controls cohort organisation 
    !******************************************************************************

    if(hlm_use_ed_st3.eq.ifalse) then 
       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))                 
          
          ! adds small cohort of each PFT
          call recruitment(currentSite, currentPatch, bc_in)
          
          currentPatch => currentPatch%younger
       enddo
    end if
    
       
    call ed_total_balance_check(currentSite,1)

    if( hlm_use_ed_st3.eq.ifalse ) then 
       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))
          
          ! puts cohorts in right order
          call sort_cohorts(currentPatch)            

          ! kills cohorts that are too few
          call terminate_cohorts(currentSite, currentPatch, 1)

          ! fuses similar cohorts
          call fuse_cohorts(currentSite,currentPatch, bc_in )
          
          ! kills cohorts for various other reasons
          call terminate_cohorts(currentSite, currentPatch, 2)
          
          
          currentPatch => currentPatch%younger
       enddo
    end if
       
    call ed_total_balance_check(currentSite,2)

    !*********************************************************************************
    ! Patch dynamics sub-routines: fusion, new patch creation (spwaning), termination.
    !*********************************************************************************

    ! make new patches from disturbed land
    if ( hlm_use_ed_st3.eq.ifalse ) then
       call spawn_patches(currentSite, bc_in)
    end if
   
    call ed_total_balance_check(currentSite,3)

    ! fuse on the spawned patches.
    if ( hlm_use_ed_st3.eq.ifalse ) then
       call fuse_patches(currentSite, bc_in )        
       
       ! If using BC FATES hydraulics, update the rhizosphere geometry
       ! based on the new cohort-patch structure
       ! 'rhizosphere geometry' (column-level root biomass + rootfr --> root length 
       ! density --> node radii and volumes)
       if( (hlm_use_planthydro.eq.itrue) .and. do_growthrecruiteffects) then
          call updateSizeDepRhizHydProps(currentSite, bc_in)
          call updateSizeDepRhizHydStates(currentSite, bc_in)
          !       if(nshell > 1) then  (THIS BEING CHECKED INSIDE OF the update)
          !          call updateSizeDepRhizHydStates(currentSite, c, soilstate_inst, &
          !                waterstate_inst)
          !       end if
       end if
    end if

    call ed_total_balance_check(currentSite,4)

    ! kill patches that are too small
    if ( hlm_use_ed_st3.eq.ifalse ) then
       call terminate_patches(currentSite)   
    end if
   
    call ed_total_balance_check(currentSite,5)

  end subroutine ed_ecosystem_dynamics

  !-------------------------------------------------------------------------------!
  subroutine ed_integrate_state_variables(currentSite, bc_in )
    !
    ! !DESCRIPTION:
    ! FIX(SPM,032414) refactor so everything goes through interface
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type)     , intent(inout) :: currentSite
    type(bc_in_type)        , intent(in)   :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort

    integer  :: c                     ! Counter for litter size class 
    integer  :: ft                    ! Counter for PFT
    real(r8) :: small_no              ! to circumvent numerical errors that cause negative values of things that can't be negative
    real(r8) :: cohort_biomass_store  ! remembers the biomass in the cohort for balance checking
    real(r8) :: dbh_old               ! dbh of plant before daily PRT [cm]
    real(r8) :: hite_old              ! height of plant before daily PRT [m]
    logical  :: is_drought            ! logical for if the plant (site) is in a drought state
    real(r8) :: leaf_c
    real(r8) :: delta_dbh             ! correction for dbh
    real(r8) :: delta_hite            ! correction for hite

    !-----------------------------------------------------------------------

    small_no = 0.0000000000_r8  ! Obviously, this is arbitrary.  RF - changed to zero

    currentSite%dseed_dt(:)       = 0._r8
    currentSite%seed_rain_flux(:) = 0._r8  

    
    
    currentPatch => currentSite%youngest_patch

    do while(associated(currentPatch))

       currentPatch%age = currentPatch%age + hlm_freq_day
       ! FIX(SPM,032414) valgrind 'Conditional jump or move depends on uninitialised value'
       if( currentPatch%age  <  0._r8 )then
          write(fates_log(),*) 'negative patch age?',currentPatch%age, &
               currentPatch%patchno,currentPatch%area
       endif

       ! add age increment to secondary forest patches as well
       if (currentPatch%anthro_disturbance_label .eq. secondaryforest) then
          currentPatch%age_since_anthro_disturbance = &
               currentPatch%age_since_anthro_disturbance + hlm_freq_day
       endif

       ! check to see if the patch has moved to the next age class
       currentPatch%age_class = get_age_class_index(currentPatch%age)

       ! Update Canopy Biomass Pools
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 


          ft = currentCohort%pft

          ! Calculate the mortality derivatives
          call Mortality_Derivative( currentSite, currentCohort, bc_in )

          ! -----------------------------------------------------------------------------
          ! Apply Plant Allocation and Reactive Transport
          ! -----------------------------------------------------------------------------

          

          ! -----------------------------------------------------------------------------
          !  Identify the net carbon gain for this dynamics interval
          !    Set the available carbon pool, identify allocation portions, and 
          !    decrement the available carbon pool to zero.
          ! -----------------------------------------------------------------------------
          !
          ! convert from kgC/indiv/day into kgC/indiv/year
          ! <x>_acc_hold is remembered until the next dynamics step (used for I/O)
          ! <x>_acc will be reset soon and will be accumulated on the next leaf 
          !         photosynthesis step
          ! -----------------------------------------------------------------------------
          
          if (hlm_use_ed_prescribed_phys .eq. itrue) then
             if (currentCohort%canopy_layer .eq. 1) then
                currentCohort%npp_acc_hold = EDPftvarcon_inst%prescribed_npp_canopy(ft) &
                     * currentCohort%c_area / currentCohort%n
                ! add these for balance checking purposes
                currentCohort%npp_acc = currentCohort%npp_acc_hold / hlm_days_per_year 
             else
                currentCohort%npp_acc_hold = EDPftvarcon_inst%prescribed_npp_understory(ft) &
                     * currentCohort%c_area / currentCohort%n
                ! add these for balance checking purposes
                currentCohort%npp_acc = currentCohort%npp_acc_hold / hlm_days_per_year
             endif
          else
             currentCohort%npp_acc_hold  = currentCohort%npp_acc  * real(hlm_days_per_year,r8)
             currentCohort%gpp_acc_hold  = currentCohort%gpp_acc  * real(hlm_days_per_year,r8)
             currentCohort%resp_acc_hold = currentCohort%resp_acc * real(hlm_days_per_year,r8)
          endif
          
          currentSite%flux_in = currentSite%flux_in + currentCohort%npp_acc * currentCohort%n

          leaf_c = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
          currentCohort%treelai = tree_lai(leaf_c, currentCohort%pft, currentCohort%c_area, currentCohort%n, &
               currentCohort%canopy_layer, currentPatch%canopy_layer_tlai, &
               currentCohort%vcmax25top)
          currentCohort%treesai = tree_sai(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_trim, &
               currentCohort%c_area, currentCohort%n, currentCohort%canopy_layer, &
               currentPatch%canopy_layer_tlai, currentCohort%treelai,currentCohort%vcmax25top,6 )
          

          ! Conduct Maintenance Turnover (parteh)
          call currentCohort%prt%CheckMassConservation(ft,3)
          if(any(currentSite%dstatus == [phen_dstat_moiston,phen_dstat_timeon])) then
             is_drought = .false.
          else
             is_drought = .true.
          end if
          call PRTMaintTurnover(currentCohort%prt,ft,is_drought)
          call currentCohort%prt%CheckMassConservation(ft,4)

          leaf_c = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
          currentCohort%treelai = tree_lai(leaf_c, currentCohort%pft, currentCohort%c_area, currentCohort%n, &
               currentCohort%canopy_layer, currentPatch%canopy_layer_tlai, &
               currentCohort%vcmax25top)
          currentCohort%treesai = tree_sai(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_trim, &
               currentCohort%c_area, currentCohort%n, currentCohort%canopy_layer, &
               currentPatch%canopy_layer_tlai, currentCohort%treelai,currentCohort%vcmax25top,7 )
          

          ! If the current diameter of a plant is somehow less than what is consistent
          ! with what is allometrically consistent with the stuctural biomass, then
          ! correct the dbh to match.

          call EvaluateAndCorrectDBH(currentCohort,delta_dbh,delta_hite)

          hite_old = currentCohort%hite
          dbh_old  = currentCohort%dbh

          ! Conduct Growth (parteh)
          call currentCohort%prt%DailyPRT()
          call currentCohort%prt%CheckMassConservation(ft,5)

          ! Update the leaf biophysical rates based on proportion of leaf
          ! mass in the different leaf age classes. Following growth
          ! and turnover, these proportions won't change again. This
          ! routine is also called following fusion
          call UpdateCohortBioPhysRates(currentCohort)

          leaf_c = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
          currentCohort%treelai = tree_lai(leaf_c, currentCohort%pft, currentCohort%c_area, currentCohort%n, &
               currentCohort%canopy_layer, currentPatch%canopy_layer_tlai, &
               currentCohort%vcmax25top)
          currentCohort%treesai = tree_sai(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_trim, &
               currentCohort%c_area, currentCohort%n, currentCohort%canopy_layer, &
               currentPatch%canopy_layer_tlai, currentCohort%treelai,currentCohort%vcmax25top,3 ) 
          
          ! Transfer all reproductive tissues into seed production
          call PRTReproRelease(currentCohort%prt,repro_organ,carbon12_element, &
                               1.0_r8, currentCohort%seed_prod)
          currentCohort%seed_prod =  currentCohort%seed_prod / hlm_freq_day

          ! This cohort has grown, it is no longer "new"
          currentCohort%isnew = .false.
          
          ! Update the plant height (if it has grown)
          call h_allom(currentCohort%dbh,ft,currentCohort%hite)
          
          currentCohort%dhdt      = (currentCohort%hite-hite_old)/hlm_freq_day
          currentCohort%ddbhdt    = (currentCohort%dbh-dbh_old)/hlm_freq_day

          ! Carbon assimilate has been spent at this point
          ! and can now be safely zero'd

          currentCohort%npp_acc  = 0.0_r8
          currentCohort%gpp_acc  = 0.0_r8
          currentCohort%resp_acc = 0.0_r8
          
          ! BOC...update tree 'hydraulic geometry' 
          ! (size --> heights of elements --> hydraulic path lengths --> 
          ! maximum node-to-node conductances)
          if( (hlm_use_planthydro.eq.itrue) .and. do_growthrecruiteffects) then
             call updateSizeDepTreeHydProps(currentSite,currentCohort, bc_in)
             call updateSizeDepTreeHydStates(currentSite,currentCohort)
          end if
  
          currentCohort => currentCohort%taller

       enddo
      
       call non_canopy_derivs( currentSite, currentPatch, bc_in)

       !update state variables simultaneously according to derivatives for this time period. 

       ! first update the litter variables that are tracked at the patch level
       do c = 1,ncwd
          currentPatch%cwd_ag(c) =  currentPatch%cwd_ag(c) + currentPatch%dcwd_ag_dt(c)* hlm_freq_day
          currentPatch%cwd_bg(c) =  currentPatch%cwd_bg(c) + currentPatch%dcwd_bg_dt(c)* hlm_freq_day
       enddo

       do ft = 1,numpft
          currentPatch%leaf_litter(ft) = currentPatch%leaf_litter(ft) + &
               currentPatch%dleaf_litter_dt(ft)* hlm_freq_day
          currentPatch%root_litter(ft) = currentPatch%root_litter(ft) + &
               currentPatch%droot_litter_dt(ft)* hlm_freq_day
       enddo

       do c = 1,ncwd
          if(currentPatch%cwd_ag(c)<small_no)then
            write(fates_log(),*) 'negative CWD_AG', currentPatch%cwd_ag(c),CurrentSite%lat,currentSite%lon
            currentPatch%cwd_ag(c) = small_no
          endif
          if(currentPatch%cwd_bg(c)<small_no)then
            write(fates_log(),*) 'negative CWD_BG', currentPatch%cwd_bg(c),CurrentSite%lat,CurrentSite%lon
            currentPatch%cwd_bg(c) = small_no
          endif
       enddo

       do ft = 1,numpft
          if(currentPatch%leaf_litter(ft)<small_no)then
            write(fates_log(),*) 'negative leaf litter numerical error', &
                  currentPatch%leaf_litter(ft),CurrentSite%lat,CurrentSite%lon,&
            currentPatch%dleaf_litter_dt(ft),currentPatch%leaf_litter_in(ft), &
            currentPatch%leaf_litter_out(ft),currentpatch%age
            currentPatch%leaf_litter(ft) = small_no
          endif
          if(currentPatch%root_litter(ft)<small_no)then
               write(fates_log(),*) 'negative root litter numerical error', currentPatch%root_litter(ft), &
               currentPatch%droot_litter_dt(ft)* hlm_freq_day, &
               CurrentSite%lat,CurrentSite%lon
            currentPatch%root_litter(ft) = small_no
          endif
       enddo

     
       ! update cohort number. This needs to happen after the CWD_input and seed_input calculations as they 
       ! assume the pre-mortality currentCohort%n. 
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 
         currentCohort%n = max(small_no,currentCohort%n + currentCohort%dndt * hlm_freq_day )  
         currentCohort => currentCohort%taller
       enddo

       currentPatch => currentPatch%older

    enddo

    ! at the site level, update the seed bank mass
    do ft = 1,numpft
       currentSite%seed_bank(ft) = currentSite%seed_bank(ft) + currentSite%dseed_dt(ft)*hlm_freq_day
    enddo

    ! Check for negative values. Write out warning to show carbon balance. 
    do ft = 1,numpft
       if(currentSite%seed_bank(ft)<small_no)then
          write(fates_log(),*) 'negative seedbank', currentSite%seed_bank(ft)
          currentSite%seed_bank(ft) = small_no
       endif
    enddo

  end subroutine ed_integrate_state_variables

  !-------------------------------------------------------------------------------!
  subroutine ed_update_site( currentSite, bc_in )
    !
    ! !DESCRIPTION:
    ! Calls routines to consolidate the ED growth process.
    ! Canopy Structure to assign canopy layers to cohorts
    ! Canopy Spread to figure out the size of tree crowns
    ! Trim_canopy to figure out the target leaf biomass. 
    ! Extra recruitment to fill empty patches.  
    !
    ! !USES:
    use EDCanopyStructureMod , only : canopy_spread, canopy_structure
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout), target :: currentSite
    type(bc_in_type)        , intent(in)             :: bc_in
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer :: currentPatch   
    integer :: cohort_number ! To print out the number of cohorts.  
    integer :: g             ! Counter for sites
    !-----------------------------------------------------------------------

    call canopy_spread(currentSite)

    call ed_total_balance_check(currentSite,6)

    call canopy_structure(currentSite, bc_in)

    call ed_total_balance_check(currentSite,7)

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))

       ! Is termination really needed here? canopy_structure just called it several times! (rgk)
       call terminate_cohorts(currentSite, currentPatch, 1) 
       call terminate_cohorts(currentSite, currentPatch, 2) 

       ! FIX(SPM,040314) why is this needed for BFB restarts? Look into this at some point
       cohort_number = count_cohorts(currentPatch)  
       if ( debug ) then
          write(fates_log(),*) 'tempCount ',cohort_number
       endif

       ! Note (RF)
       ! This breaks the balance check, but if we leave it out, then 
       ! the first new patch that isn't fused has no cohorts at the end of the spawn process
       ! and so there are radiation errors instead. 
       ! Fixing this would likely require a re-work of how seed germination works which would be tricky. 
       if(currentPatch%countcohorts < 1)then
          !write(fates_log(),*) 'ED: calling recruitment for no cohorts',currentSite%clmgcell,currentPatch%patchno
          !call recruitment(1, currentSite, currentPatch)
          ! write(fates_log(),*) 'patch empty',currentPatch%area,currentPatch%age
       endif

       currentPatch => currentPatch%younger    

    enddo

    ! FIX(RF,032414). This needs to be monthly, not annual
    ! If this is the second to last day of the year, then perform trimming
    if( hlm_day_of_year == hlm_days_per_year-1) then

       write(fates_log(),*) 'calling trim canopy' 
       call trim_canopy(currentSite)  
    endif

  end subroutine ed_update_site

  !-------------------------------------------------------------------------------!
  subroutine ed_total_balance_check (currentSite, call_index )
    !
    ! !DESCRIPTION:
    ! This routine looks at the carbon in and out of the ED model and compares it to 
    ! the change in total carbon stocks. 
    ! Fluxes in are NPP. Fluxes out are decay of CWD and litter into SOM pools.  
    ! ed_allsites_inst%flux_out and ed_allsites_inst%flux_in are set where they occur 
    ! in the code. 
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout) :: currentSite
    integer            , intent(in)    :: call_index
    !
    ! !LOCAL VARIABLES:
    real(r8) :: biomass_stock   ! total biomass   in KgC/site
    real(r8) :: litter_stock    ! total litter    in KgC/site
    real(r8) :: seed_stock      ! total seed mass in KgC/site
    real(r8) :: total_stock     ! total ED carbon in KgC/site
    real(r8) :: change_in_stock ! Change since last time we set ed_allsites_inst%old_stock in this routine.  KgC/site
    real(r8) :: error           ! How much carbon did we gain or lose (should be zero!) 
    real(r8) :: error_frac      ! Error as a fraction of total biomass
    real(r8) :: net_flux        ! Difference between recorded fluxes in and out. KgC/site
    real(r8) :: leaf_c
    real(r8) :: fnrt_c
    real(r8) :: sapw_c
    real(r8) :: store_c
    real(r8) :: struct_c
    real(r8) :: repro_c

    ! nb. There is no time associated with these variables 
    ! because this routine can be called between any two 
    ! arbitrary points in code, even if no time has passed. 
    ! Also, the carbon pools are per site/gridcell, so that 
    ! we can account for the changing areas of patches. 

    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort
    !-----------------------------------------------------------------------

    change_in_stock = 0.0_r8


    call SiteCarbonStock(currentSite,total_stock,biomass_stock,litter_stock,seed_stock)


    change_in_stock = total_stock - currentSite%old_stock  

    net_flux        = currentSite%flux_in - currentSite%flux_out
    error           = abs(net_flux - change_in_stock)   

    if(change_in_stock>0.0)then
       error_frac      = error/abs(total_stock)
    else
       error_frac      = 0.0_r8
    end if

    ! -----------------------------------------------------------------------------------
    ! Terms:
    ! %flux_in:  accumulates npp over all cohorts,  
    !               currentSite%flux_in = currentSite%flux_in + &
    !               currentCohort%npp_acc * currentCohort%n
    ! %flux_out: coarse woody debris going into fragmentation pools:
    !               currentSite%flux_out + sum(currentPatch%leaf_litter_out) * &
    !               currentPatch%area *hlm_freq_day!kgC/site/day
    !            burn fractions:  
    !               currentSite%flux_out = currentSite%flux_out + &
    !               burned_litter * new_patch%area !kG/site/day
    ! -----------------------------------------------------------------------------------
    
    if ( error_frac > 10e-6_r8 ) then
       write(fates_log(),*) 'carbon balance error detected'
       write(fates_log(),*) 'error fraction relative to biomass stock:',error_frac
       write(fates_log(),*) 'call index: ',call_index
       write(fates_log(),*) 'flux in (npp):  ',currentSite%flux_in
       write(fates_log(),*) 'flux out (fragmentation/harvest): ',currentSite%flux_out
       write(fates_log(),*) 'net: ',net_flux
       write(fates_log(),*) 'dstock: ',change_in_stock
       write(fates_log(),*) 'error=net_flux-dstock:', error
       write(fates_log(),*) 'biomass', biomass_stock
       write(fates_log(),*) 'litter',litter_stock
       write(fates_log(),*) 'seeds',seed_stock
       write(fates_log(),*) 'previous total',currentSite%old_stock  

       write(fates_log(),*) 'lat lon',currentSite%lat,currentSite%lon
       
       ! If this is the first day of simulation, carbon balance reports but does not end the run
       if( int(hlm_current_year*10000 + hlm_current_month*100 + hlm_current_day).ne.hlm_reference_date ) then

          change_in_stock = 0.0_r8
          biomass_stock   = 0.0_r8
          litter_stock    = 0.0_r8
          
          seed_stock   =  sum(currentSite%seed_bank)*AREA
          currentPatch => currentSite%oldest_patch
          do while(associated(currentPatch))
             write(fates_log(),*) '---------------------------------------'
             write(fates_log(),*) currentPatch%area , sum(currentPatch%cwd_ag), sum(currentPatch%cwd_bg)
             write(fates_log(),*) sum(currentPatch%leaf_litter),sum(currentPatch%root_litter)
             write(fates_log(),*)'---'
             currentCohort => currentPatch%tallest
             do while(associated(currentCohort))
                write(fates_log(),*) 'pft: ',currentCohort%pft
                write(fates_log(),*) 'dbh: ',currentCohort%dbh
                leaf_c   = currentCohort%prt%GetState(leaf_organ,all_carbon_elements)
                struct_c = currentCohort%prt%GetState(struct_organ,all_carbon_elements)
                store_c  = currentCohort%prt%GetState(store_organ,all_carbon_elements)
                fnrt_c   = currentCohort%prt%GetState(fnrt_organ,all_carbon_elements)
                repro_c  = currentCohort%prt%GetState(repro_organ,all_carbon_elements)
                sapw_c   = currentCohort%prt%GetState(sapw_organ,all_carbon_elements)

                write(fates_log(),*) 'lc: ',leaf_c,' dc: ',struct_c,' stc: ',store_c
                write(fates_log(),*) 'fc: ',fnrt_c,' rc: ',repro_c,' sac: ',sapw_c
                write(fates_log(),*) 'N plant: ',currentCohort%n
                currentCohort => currentCohort%shorter
             enddo !end cohort loop                                                                                                   
             currentPatch => currentPatch%younger
          enddo !end patch loop                                                                                                      
          write(fates_log(),*) 'aborting on date:',hlm_current_year,hlm_current_month,hlm_current_day
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if

    endif

    currentSite%flux_in   = 0.0_r8
    currentSite%flux_out  = 0.0_r8  
    currentSite%old_stock = total_stock

 end subroutine ed_total_balance_check
 
 ! =====================================================================================
 
 subroutine bypass_dynamics(currentSite)

    ! ----------------------------------------------------------------------------------
    ! If dynamics are bypassed, various fluxes, rates and flags need to be set
    ! to trivial values.
    ! WARNING: Turning off things like dynamics is experimental. The setting of
    ! variables to trivial values may not be complete, use at your own risk.
    ! ----------------------------------------------------------------------------------

    ! Arguments
    type(ed_site_type)      , intent(inout), target  :: currentSite
    
    ! Locals
    type(ed_patch_type), pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort
    
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 

          currentCohort%isnew=.false.

          currentCohort%npp_acc_hold  = currentCohort%npp_acc  * real(hlm_days_per_year,r8)
          currentCohort%gpp_acc_hold  = currentCohort%gpp_acc  * real(hlm_days_per_year,r8)
          currentCohort%resp_acc_hold = currentCohort%resp_acc * real(hlm_days_per_year,r8)

          currentCohort%npp_acc  = 0.0_r8
          currentCohort%gpp_acc  = 0.0_r8
          currentCohort%resp_acc = 0.0_r8

          ! No need to set the "net_art" terms to zero
          ! they are zeroed at the beginning of the daily step
          ! If DailyPRT, maintenance, and phenology are not called
          ! then these should stay zero.

          currentCohort%bmort = 0.0_r8
          currentCohort%hmort = 0.0_r8
          currentCohort%cmort = 0.0_r8
          currentCohort%frmort = 0.0_r8

          currentCohort%dndt      = 0.0_r8
          currentCohort%dhdt      = 0.0_r8
          currentCohort%ddbhdt    = 0.0_r8

          currentCohort => currentCohort%taller
       enddo
       currentPatch => currentPatch%older
    enddo
    
 end subroutine bypass_dynamics


end module EDMainMod
