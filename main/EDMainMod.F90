module EDMainMod

  ! ===========================================================================
  ! Main ED module.
  ! ============================================================================

  use shr_kind_mod             , only : r8 => shr_kind_r8

  use FatesGlobals             , only : fates_log

  use FatesInterfaceTypesMod        , only : hlm_freq_day
  use FatesInterfaceTypesMod        , only : hlm_day_of_year
  use FatesInterfaceTypesMod        , only : hlm_days_per_year
  use FatesInterfaceTypesMod        , only : hlm_current_year
  use FatesInterfaceTypesMod        , only : hlm_current_month
  use FatesInterfaceTypesMod        , only : hlm_current_day
  use FatesInterfaceTypesMod        , only : hlm_use_planthydro
  use FatesInterfaceTypesMod        , only : hlm_parteh_mode
  use FatesInterfaceTypesMod        , only : hlm_use_cohort_age_tracking
  use FatesInterfaceTypesMod        , only : hlm_reference_date
  use FatesInterfaceTypesMod        , only : hlm_use_ed_prescribed_phys
  use FatesInterfaceTypesMod        , only : hlm_use_tree_damage
  use FatesInterfaceTypesMod        , only : hlm_use_ed_st3
  use FatesInterfaceTypesMod        , only : hlm_use_sp
  use FatesInterfaceTypesMod        , only : bc_in_type
  use FatesInterfaceTypesMod        , only : bc_out_type
  use FatesInterfaceTypesMod        , only : hlm_masterproc
  use FatesInterfaceTypesMod        , only : numpft
  use FatesInterfaceTypesMod        , only : hlm_use_nocomp
  use PRTGenericMod            , only : prt_carbon_allom_hyp
  use PRTGenericMod            , only : prt_cnp_flex_allom_hyp
  use PRTGenericMod            , only : nitrogen_element
  use PRTGenericMod            , only : phosphorus_element
  use EDCohortDynamicsMod      , only : terminate_cohorts
  use EDCohortDynamicsMod      , only : fuse_cohorts
  use EDCohortDynamicsMod      , only : EvaluateAndCorrectDBH
  use EDCohortDynamicsMod      , only : DamageRecovery
  use EDPatchDynamicsMod       , only : disturbance_rates
  use EDPatchDynamicsMod       , only : fuse_patches
  use EDPatchDynamicsMod       , only : spawn_patches
  use EDPatchDynamicsMod       , only : terminate_patches
  use EDPhysiologyMod          , only : phenology
  use EDPhysiologyMod          , only : satellite_phenology
  use EDPhysiologyMod          , only : recruitment
  use EDPhysiologyMod          , only : trim_canopy
  use EDPhysiologyMod          , only : SeedUpdate
  use EDPhysiologyMod          , only : ZeroAllocationRates
  use EDPhysiologyMod          , only : ZeroLitterFluxes
  use EDPhysiologyMod          , only : PreDisturbanceLitterFluxes
  use EDPhysiologyMod          , only : PreDisturbanceIntegrateLitter
  use EDPhysiologyMod          , only : UpdateRecruitL2FR
  use EDPhysiologyMod          , only : UpdateRecruitStoich
  use EDPhysiologyMod          , only : SetRecruitL2FR
  use EDPhysiologyMod          , only : GenerateDamageAndLitterFluxes
  use FatesSoilBGCFluxMod      , only : FluxIntoLitterPools
  use FatesSoilBGCFluxMod      , only : EffluxIntoLitterPools
  use FatesSoilBGCFluxMod      , only : PrepNutrientAquisitionBCs
  use FatesSoilBGCFluxMod      , only : PrepCH4BCs
  use SFMainMod                , only : DailyFireModel
  use FatesSizeAgeTypeIndicesMod, only : get_age_class_index
  use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index
  use FatesLitterMod           , only : litter_type
  use FatesLitterMod           , only : ncwd
  use EDtypesMod               , only : ed_site_type
  use EDTypesMod               , only : set_patchno
  use FatesPatchMod            , only : fates_patch_type
  use FatesCohortMod           , only : fates_cohort_type
  use EDTypesMod               , only : AREA
  use EDTypesMod               , only : site_massbal_type
  use PRTGenericMod            , only : num_elements
  use PRTGenericMod            , only : element_list
  use PRTGenericMod            , only : element_pos
  use EDTypesMod               , only : phen_dstat_moiston
  use EDTypesMod               , only : phen_dstat_timeon
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesConstantsMod        , only : primaryland, secondaryland
  use FatesConstantsMod        , only : n_landuse_cats  
  use FatesConstantsMod        , only : nearzero
  use FatesConstantsMod        , only : m2_per_ha
  use FatesConstantsMod        , only : sec_per_day
  use FatesConstantsMod        , only : nocomp_bareground
  use FatesPlantHydraulicsMod  , only : do_growthrecruiteffects
  use FatesPlantHydraulicsMod  , only : UpdateSizeDepPlantHydProps
  use FatesPlantHydraulicsMod  , only : UpdateSizeDepPlantHydStates
  use FatesPlantHydraulicsMod  , only : InitPlantHydStates
  use FatesPlantHydraulicsMod  , only : UpdateSizeDepRhizHydProps
  use FatesPlantHydraulicsMod  , only : AccumulateMortalityWaterStorage
  use FatesAllometryMod        , only : h_allom
  use EDLoggingMortalityMod    , only : IsItLoggingTime
  use EDLoggingMortalityMod    , only : get_harvestable_carbon
  use DamageMainMod            , only : IsItDamageTime
  use FatesGlobals             , only : endrun => fates_endrun
  use ChecksBalancesMod        , only : SiteMassStock
  use ChecksBalancesMod        , only : CheckIntegratedMassPools
  use EDMortalityFunctionsMod  , only : Mortality_Derivative
  use EDTypesMod               , only : AREA_INV
  use PRTGenericMod,          only : carbon12_element
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTLossFluxesMod,       only : PRTMaintTurnover
  use PRTParametersMod      , only : prt_params
  use EDPftvarcon,            only : EDPftvarcon_inst
  use FatesHistoryInterfaceMod, only : fates_hist
  use FatesLandUseChangeMod,  only: FatesGrazing

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
  private :: TotalBalanceCheck
  private :: bypass_dynamics

  logical :: debug  = .false.

  integer, parameter :: final_check_id = -1

  character(len=*), parameter, private :: sourcefile = &
         __FILE__

  !
  ! 10/30/09: Created by Rosie Fisher
  !-----------------------------------------------------------------------

contains

  !-------------------------------------------------------------------------------!
  subroutine ed_ecosystem_dynamics(currentSite, bc_in, bc_out)
    !
    ! !DESCRIPTION:
    !  Core of ed model, calling all subsequent vegetation dynamics routines
    !
    ! !ARGUMENTS:
    type(ed_site_type)      , intent(inout), target  :: currentSite
    type(bc_in_type)        , intent(in)             :: bc_in
    type(bc_out_type)       , intent(inout)          :: bc_out
    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type), pointer :: currentPatch
    integer :: el                ! Loop counter for variables 
    integer :: do_patch_dynamics ! for some modes, we turn off patch dynamics

    !-----------------------------------------------------------------------

    if (debug .and.( hlm_masterproc==itrue)) write(fates_log(),'(A,I4,A,I2.2,A,I2.2)') 'FATES Dynamics: ',&
          hlm_current_year,'-',hlm_current_month,'-',hlm_current_day

    ! Consider moving this towards the end, because some of these
    ! are being integrated over the short time-step

    do el = 1,num_elements
       call currentSite%mass_balance(el)%ZeroMassBalFlux()
    end do
    call currentSite%flux_diags%ZeroFluxDiags()

    
    ! Call a routine that simply identifies if logging should occur
    ! This is limited to a global event until more structured event handling is enabled
    call IsItLoggingTime(hlm_masterproc,currentSite)

    ! Call a routine that identifies if damage should occur
    call IsItDamageTime(hlm_masterproc)
 
    !**************************************************************************
    ! Fire, growth, biogeochemistry.
    !**************************************************************************

    !FIX(SPM,032414) take this out.  On startup these values are all zero and on restart it
    !zeros out values read in the restart file

    ! Zero turnover rates and growth diagnostics
    call ZeroAllocationRates(currentSite)

    ! Zero fluxes in and out of litter pools
    call ZeroLitterFluxes(currentSite)

    ! Zero mass balance
    call TotalBalanceCheck(currentSite, 0)

    ! We do not allow phenology while in ST3 mode either, it is hypothetically
    ! possible to allow this, but we have not plugged in the litter fluxes
    ! of flushing or turning over leaves for non-dynamics runs
    if (hlm_use_ed_st3.eq.ifalse)then
      if(hlm_use_sp.eq.ifalse) then
        call phenology(currentSite, bc_in )
      else
        call satellite_phenology(currentSite, bc_in )
      end if ! SP phenology
    end if


    if (hlm_use_ed_st3.eq.ifalse.and.hlm_use_sp.eq.ifalse) then   ! Bypass if ST3
       
       ! Check that the site doesn't consist solely of a single bareground patch.
       ! If so, skip the fire model.  Since the bareground patch should be the
       ! oldest patch per set_patchno, we check that the youngest patch isn't zero.
       ! If there are multiple patches on the site, the bareground patch is avoided
       ! at the level of the fire_model subroutines.

       if (currentSite%youngest_patch%nocomp_pft_label .ne. nocomp_bareground)then
          call DailyFireModel(currentSite, bc_in)
       end if

       ! Calculate disturbance and mortality based on previous timestep vegetation.
       ! disturbance_rates calls logging mortality and other mortalities, Yi Xu
       call disturbance_rates(currentSite, bc_in)

       ! Integrate state variables from annual rates to daily timestep
       call ed_integrate_state_variables(currentSite, bc_in, bc_out )

       ! at this point in the call sequence, if flag to transition_landuse_from_off_to_on was set, unset it as it is no longer needed
       if(currentSite%transition_landuse_from_off_to_on) then
          currentSite%transition_landuse_from_off_to_on = .false.
       endif
       
    else
       ! ed_intergrate_state_variables is where the new cohort flag
       ! is set. This flag designates wether a cohort has
       ! experienced a day, and therefore has been populated with non-nonsense
       ! values.  If we aren't entering that sequence, we need to set the flag
       ! Make sure cohorts are marked as non-recruits

       call bypass_dynamics(currentSite,bc_out)

    end if

    !******************************************************************************
    ! Reproduction, Recruitment and Cohort Dynamics : controls cohort organization
    !******************************************************************************

    if(hlm_use_ed_st3.eq.ifalse.and.hlm_use_sp.eq.ifalse) then
       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))

          ! adds small cohort of each PFT
          call recruitment(currentSite, currentPatch, bc_in)
          !YL --------------
          ! call recruitment(currentSite, currentPatch, bc_in, bc_out)

          currentPatch => currentPatch%younger
       enddo

       call TotalBalanceCheck(currentSite,1)

       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))

          ! puts cohorts in right order
          call currentPatch%SortCohorts()

          ! kills cohorts that are too few
          call terminate_cohorts(currentSite, currentPatch, 1, 10, bc_in  )

          ! fuses similar cohorts
          call fuse_cohorts(currentSite,currentPatch, bc_in )

          ! kills cohorts for various other reasons
          call terminate_cohorts(currentSite, currentPatch, 2, 10, bc_in )


          currentPatch => currentPatch%younger
       enddo
         
    end if

    call TotalBalanceCheck(currentSite,2)

    !*********************************************************************************
    ! Patch dynamics sub-routines: fusion, new patch creation (spwaning), termination.
    !*********************************************************************************

    ! turn off patch dynamics if SP or ST3 modes in use
    do_patch_dynamics = itrue
    if(hlm_use_ed_st3.eq.itrue .or. &
       hlm_use_sp.eq.itrue)then
       do_patch_dynamics = ifalse
    end if

    ! make new patches from disturbed land
    if (do_patch_dynamics.eq.itrue ) then

       call spawn_patches(currentSite, bc_in)

       call TotalBalanceCheck(currentSite,3)

       ! fuse on the spawned patches.
       call fuse_patches(currentSite, bc_in )

       ! If using BC FATES hydraulics, update the rhizosphere geometry
       ! based on the new cohort-patch structure
       ! 'rhizosphere geometry' (column-level root biomass + rootfr --> root length
       ! density --> node radii and volumes)
       if( (hlm_use_planthydro.eq.itrue) .and. do_growthrecruiteffects) then
          call UpdateSizeDepRhizHydProps(currentSite, bc_in)
          !! call UpdateSizeDepRhizHydStates(currentSite, bc_in) ! keeping if re-implemented (RGK 12-2021)
       end if

       ! SP has changes in leaf carbon but we don't expect them to be in balance.
       call TotalBalanceCheck(currentSite,4)

       ! kill patches that are too small
       call terminate_patches(currentSite, bc_in)
    end if

    ! Final instantaneous mass balance check
    call TotalBalanceCheck(currentSite,5)

    
  end subroutine ed_ecosystem_dynamics

  !-------------------------------------------------------------------------------!
  subroutine ed_integrate_state_variables(currentSite, bc_in, bc_out )
    !

    ! !DESCRIPTION:
    ! FIX(SPM,032414) refactor so everything goes through interface
    !
    ! !USES:
    use FatesInterfaceTypesMod, only : hlm_num_lu_harvest_cats
    use PRTGenericMod        , only : leaf_organ
    use PRTGenericMod        , only : repro_organ
    use PRTGenericMod        , only : sapw_organ
    use PRTGenericMod        , only : struct_organ
    use PRTGenericMod        , only : store_organ
    use PRTGenericMod        , only : fnrt_organ
    use FatesInterfaceTypesMod, only : hlm_use_cohort_age_tracking
    use FatesConstantsMod, only : itrue
    use FatesConstantsMod     , only : nearzero
    use EDCanopyStructureMod  , only : canopy_structure


    ! !ARGUMENTS:

    type(ed_site_type)     , intent(inout) :: currentSite
    type(bc_in_type)        , intent(in)   :: bc_in
    type(bc_out_type)       , intent(inout)  :: bc_out

    !
    ! !LOCAL VARIABLES:
    type(site_massbal_type), pointer :: site_cmass
    type(fates_patch_type)  , pointer :: currentPatch
    type(fates_cohort_type) , pointer :: currentCohort
    type(fates_cohort_type) , pointer :: nc
    type(fates_cohort_type) , pointer :: storesmallcohort
    type(fates_cohort_type) , pointer :: storebigcohort
    
    integer :: snull
    integer :: tnull 

    integer  :: c                     ! Counter for litter size class
    integer  :: ft                    ! Counter for PFT
    integer  :: io_si                 ! global site index for history writing
    integer  :: iscpf                 ! index for the size-class x pft multiplexed bins
    integer  :: el                    ! Counter for element type (c,n,p,etc)
    real(r8) :: cohort_biomass_store  ! remembers the biomass in the cohort for balance checking
    real(r8) :: dbh_old               ! dbh of plant before daily PRT [cm]
    real(r8) :: height_old            ! height of plant before daily PRT [m]
    logical  :: is_drought            ! logical for if the plant (site) is in a drought state
    real(r8) :: delta_dbh             ! correction for dbh
    real(r8) :: delta_height          ! correction for height
    real(r8) :: mean_temp

    logical  :: newly_recovered       ! If the current loop is dealing with a newly created cohort, which
                                      ! was created because it is a clone of the previous cohort in
                                      ! a lowered damage state. This cohort should bypass several calculations
                                      ! because it inherited them (such as daily carbon balance)
    real(r8) :: target_leaf_c
    real(r8) :: current_fates_landuse_state_vector(n_landuse_cats)

    real(r8) :: harvestable_forest_c(hlm_num_lu_harvest_cats)
    integer  :: harvest_tag(hlm_num_lu_harvest_cats)

    real(r8) :: n_old
    real(r8) :: n_recover
    real(r8) :: sapw_c
    real(r8) :: leaf_c
    real(r8) :: fnrt_c
    real(r8) :: struct_c
    real(r8) :: repro_c
    real(r8) :: total_c
    real(r8) :: store_c

    real(r8) :: cc_leaf_c
    real(r8) :: cc_fnrt_c
    real(r8) :: cc_struct_c
    real(r8) :: cc_repro_c
    real(r8) :: cc_store_c
    real(r8) :: cc_sapw_c
    
    real(r8) :: sapw_c0
    real(r8) :: leaf_c0
    real(r8) :: fnrt_c0
    real(r8) :: struct_c0
    real(r8) :: repro_c0
    real(r8) :: store_c0
    real(r8) :: total_c0
    real(r8) :: nc_carbon
    real(r8) :: cc_carbon
    
    integer,parameter :: leaf_c_id = 1
    
    !-----------------------------------------------------------------------

    current_fates_landuse_state_vector = currentSite%get_current_landuse_statevector()

    ! Clear site GPP and AR passing to HLM
    bc_out%gpp_site = 0._r8
    bc_out%ar_site = 0._r8

    ! Patch level biomass are required for C-based harvest
    call get_harvestable_carbon(currentSite, bc_in%site_area, bc_in%hlm_harvest_catnames, harvestable_forest_c)

    ! Set a pointer to this sites carbon12 mass balance
    site_cmass => currentSite%mass_balance(element_pos(carbon12_element))

    ! This call updates the assessment of the total stoichiometry
    ! for a new recruit, based on its PFT and the L2FR of
    ! a new recruit.  This is called here, because it is
    ! prior to the growth sequence, where reproductive
    ! tissues are allocated
    call UpdateRecruitStoich(currentSite)

    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))

       currentPatch%age = currentPatch%age + hlm_freq_day
       ! FIX(SPM,032414) valgrind 'Conditional jump or move depends on uninitialised value'
       if( currentPatch%age  <  0._r8 )then
          write(fates_log(),*) 'negative patch age?',currentPatch%age, &
               currentPatch%patchno,currentPatch%area
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif

       ! add age increment to secondary forest patches as well
       if (currentPatch%land_use_label .ne. primaryland) then
          currentPatch%age_since_anthro_disturbance = &
               currentPatch%age_since_anthro_disturbance + hlm_freq_day
       endif

       ! check to see if the patch has moved to the next age class
       currentPatch%age_class = get_age_class_index(currentPatch%age)


       ! Within this loop, we may be creating new cohorts, which
       ! are copies of pre-existing cohorts with reduced damage classes.
       ! If that is true, we want to bypass some of the things in
       ! this loop (such as calculation of npp, etc) because they
       ! are derived from the donor and have been modified accordingly
       newly_recovered = .false.
       
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))

          ft = currentCohort%pft
          ! Some cohorts are created and inserted to the list while
          ! the loop is going. These are pointed to the "taller" position
          ! of current, and then inherit properties of their donor (current)
          ! we don't need to repeat things before allocation for these
          ! newly_recovered cohorts
          
          if_not_newlyrecovered: if(.not.newly_recovered) then


             ! Calculate the mortality derivatives
             mean_temp = currentPatch%tveg24%GetMean()
             call Mortality_Derivative(currentSite, currentCohort, bc_in,      &
               currentPatch%btran_ft, mean_temp,                               &
               currentPatch%land_use_label,                                    &
               currentPatch%age_since_anthro_disturbance, current_fates_landuse_state_vector,   &
               harvestable_forest_c, harvest_tag)


             ! -----------------------------------------------------------------------------
             ! Apply Plant Allocation and Reactive Transport
             ! -----------------------------------------------------------------------------
             ! -----------------------------------------------------------------------------
             !  Identify the net carbon gain for this dynamics interval
             !    Set the available carbon pool, identify allocation portions, and
             !    decrement the available carbon pool to zero.
             ! -----------------------------------------------------------------------------


             if (hlm_use_ed_prescribed_phys .eq. itrue) then
                if (currentCohort%canopy_layer .eq. 1) then
                   currentCohort%npp_acc = EDPftvarcon_inst%prescribed_npp_canopy(ft) &
                        * currentCohort%c_area / currentCohort%n / hlm_days_per_year
                else
                   currentCohort%npp_acc = EDPftvarcon_inst%prescribed_npp_understory(ft) &
                        * currentCohort%c_area / currentCohort%n / hlm_days_per_year
                endif

                ! We don't explicitly define a respiration rate for prescribe phys
                ! but we do need to pass mass balance. So we say it is zero respiration
                currentCohort%gpp_acc  = currentCohort%npp_acc
                currentCohort%resp_m_acc = 0._r8

             end if

             ! -----------------------------------------------------------------------------
             ! Save NPP/GPP/R in these "hold" style variables. These variables
             ! persist after this routine is complete, and used in I/O diagnostics.
             ! Whereas the _acc style variables are zero'd because they are key
             ! accumulation state variables.
             !
             ! convert from kgC/indiv/day into kgC/indiv/year
             ! <x>_acc_hold is remembered until the next dynamics step (used for I/O)
             ! <x>_acc will be reset soon and will be accumulated on the next leaf
             !         photosynthesis step
             ! -----------------------------------------------------------------------------

             currentCohort%gpp_acc_hold  = currentCohort%gpp_acc  * real(hlm_days_per_year,r8)
             currentCohort%resp_m_acc_hold = currentCohort%resp_m_acc * real(hlm_days_per_year,r8)

             ! at this point we have the info we need to calculate growth respiration
             ! as a "tax" on the difference between daily GPP and daily maintenance respiration
             if (hlm_use_ed_prescribed_phys .eq. ifalse) then

                currentCohort%resp_g_acc_hold = prt_params%grperc(ft) * &
                     max(0._r8,(currentCohort%gpp_acc - currentCohort%resp_m_acc)) * real(hlm_days_per_year,r8)

             else
                ! set growth respiration to zero in prescribed physiology mode,
                ! that way the npp_acc vars will be set to the nominal gpp values set above.
                currentCohort%resp_g_acc_hold = 0._r8
             endif

             ! calculate the npp as the difference between gpp and autotrophic respiration
             ! (NPP is also updated if there is any excess respiration from nutrient limitations)
             currentCohort%npp_acc       = currentCohort%gpp_acc - &
                  (currentCohort%resp_m_acc + currentCohort%resp_g_acc_hold/real(hlm_days_per_year,r8))
             currentCohort%npp_acc_hold  = currentCohort%gpp_acc_hold - &
                  (currentCohort%resp_m_acc_hold + currentCohort%resp_g_acc_hold)
             

             ! allow herbivores to graze
             call FatesGrazing(currentCohort%prt, ft, currentPatch%land_use_label, currentCohort%height)

             ! Conduct Maintenance Turnover (parteh)
             if(debug) call currentCohort%prt%CheckMassConservation(ft,3)
             if(any(currentSite%dstatus(ft) == [phen_dstat_moiston,phen_dstat_timeon])) then
                is_drought = .false.
             else
                is_drought = .true.
             end if

             call PRTMaintTurnover(currentCohort%prt,ft, currentCohort%canopy_layer,is_drought)
             
             ! -----------------------------------------------------------------------------------
             ! Call the routine that advances leaves in age.
             ! This will move a portion of the leaf mass in each
             ! age bin, to the next bin. This will not handle movement
             ! of mass from the oldest bin into the litter pool, that is something else.
             ! -----------------------------------------------------------------------------------
             call currentCohort%prt%AgeLeaves(ft,currentCohort%canopy_layer, sec_per_day)

             ! Plants can acquire N from 3 sources (excluding re-absorption),
             ! the source doesn't affect how its allocated (yet), so they
             ! are combined into daily_n_gain, which is the value used in the following
             ! allocation scheme
             
             currentCohort%daily_n_gain = currentCohort%daily_nh4_uptake + &
                  currentCohort%daily_no3_uptake + currentCohort%sym_nfix_daily

             currentCohort%resp_excess_hold = 0._r8
             
          end if if_not_newlyrecovered

          ! If the current diameter of a plant is somehow less than what is consistent
          ! with what is allometrically consistent with the stuctural biomass, then
          ! correct the dbh to match.
          call EvaluateAndCorrectDBH(currentCohort,delta_dbh,delta_height)
          
          ! We want to save these values for the newly recovered cohort as well
          height_old = currentCohort%height
          dbh_old  = currentCohort%dbh

          ! -----------------------------------------------------------------------------
          ! Growth and Allocation (PARTEH)
          ! -----------------------------------------------------------------------------
          

          ! We split the allocation into phases (currently for all hypotheses)
          ! In phase 1, allocation, we address prioritized allocation that should
          ! only happen once per day, this is only allocation that does not grow stature.
          ! In phase 2, allocation , we address allocation that can be performed
          ! as many times as necessary. This is allocation that does not contain stature
          ! growth.  This is separate from phase 1, because some recovering plants
          ! will have new allocation targets that need to be updated after they change status.
          ! In Phase 3, we assume that the plant has reached its targets, and any
          ! left-over resources are used to grow the stature of the plant
          
          if(.not.newly_recovered)then
             call currentCohort%prt%DailyPRT(phase=1)
          end if

          call currentCohort%prt%DailyPRT(phase=2)
          
          if((.not.newly_recovered) .and. (hlm_use_tree_damage .eq. itrue) ) then
             ! The loop order is shortest to tallest
             ! The recovered cohort (ie one with larger targets)
             ! is newly created in DamageRecovery(), and
             ! is inserted into the next position, following the 
             ! original and current (unrecovered) cohort.
             ! we pass it back here in case the pointer is
             ! needed for diagnostics
             call DamageRecovery(currentSite,currentPatch,currentCohort,newly_recovered)

          else
             newly_recovered = .false.
          end if

          call currentCohort%prt%DailyPRT(phase=3)

          ! If nutrients are limiting growth, and carbon continues
          ! to accumulate beyond the plant's storage capacity, then
          ! it will burn carbon as what we call "excess respiration"
          ! We must subtract this term from NPP. We do not need to subtract it from
          ! currentCohort%npp_acc, it has already been removed from this in
          ! the daily growth code (PARTEH).
          
          currentCohort%npp_acc_hold  = currentCohort%npp_acc_hold - currentCohort%resp_excess_hold

          ! Passing gpp_acc_hold to HLM 
          bc_out%gpp_site = bc_out%gpp_site + currentCohort%gpp_acc_hold * &
               AREA_INV * currentCohort%n / hlm_days_per_year / sec_per_day
          bc_out%ar_site = bc_out%ar_site + (currentCohort%resp_m_acc_hold + &
               currentCohort%resp_g_acc_hold + currentCohort%resp_excess_hold) * & 
               AREA_INV * currentCohort%n / hlm_days_per_year / sec_per_day
          
          
          ! Update the mass balance tracking for the daily nutrient uptake flux
          ! Then zero out the daily uptakes, they have been used

          ! -----------------------------------------------------------------------------


          
          call EffluxIntoLitterPools(currentSite, currentPatch, currentCohort, bc_in )

          if(element_pos(nitrogen_element)>0) then
             ! Mass balance for N uptake
             currentSite%mass_balance(element_pos(nitrogen_element))%net_root_uptake = &
                  currentSite%mass_balance(element_pos(nitrogen_element))%net_root_uptake + &
                  (currentCohort%daily_n_gain-currentCohort%daily_n_efflux)*currentCohort%n
          end if
          if(element_pos(phosphorus_element)>0) then
             ! Mass balance for P uptake
             currentSite%mass_balance(element_pos(phosphorus_element))%net_root_uptake = &
                  currentSite%mass_balance(element_pos(phosphorus_element))%net_root_uptake + &
                  (currentCohort%daily_p_gain-currentCohort%daily_p_efflux)*currentCohort%n
          end if
          
          ! mass balance for C efflux (if any)
          currentSite%mass_balance(element_pos(carbon12_element))%net_root_uptake = &
               currentSite%mass_balance(element_pos(carbon12_element))%net_root_uptake - &
               currentCohort%daily_c_efflux*currentCohort%n

          ! Save NPP diagnostic for flux accounting [kg/m2/day]

          currentSite%flux_diags%npp = currentSite%flux_diags%npp + &
               currentCohort%npp_acc_hold/hlm_days_per_year * currentCohort%n * area_inv
          
          ! And simultaneously add the input fluxes to mass balance accounting
          site_cmass%gpp_acc   = site_cmass%gpp_acc + &
                currentCohort%gpp_acc * currentCohort%n

          site_cmass%aresp_acc = site_cmass%aresp_acc + &
               currentCohort%resp_m_acc*currentCohort%n + & 
               (currentCohort%resp_g_acc_hold+currentCohort%resp_excess_hold) * &
               currentCohort%n/real( hlm_days_per_year,r8)

          call currentCohort%prt%CheckMassConservation(ft,5)

          ! Update the leaf biophysical rates based on proportion of leaf
          ! mass in the different leaf age classes. Following growth
          ! and turnover, these proportions won't change again. This
          ! routine is also called following fusion
          call currentCohort%UpdateCohortBioPhysRates()

          ! This cohort has grown, it is no longer "new"
          currentCohort%isnew = .false.

          ! Update the plant height (if it has grown)
          call h_allom(currentCohort%dbh,ft,currentCohort%height)

          currentCohort%dhdt      = (currentCohort%height-height_old)/hlm_freq_day
          currentCohort%ddbhdt    = (currentCohort%dbh-dbh_old)/hlm_freq_day

          ! Carbon assimilate has been spent at this point
          ! and can now be safely zero'd

          currentCohort%npp_acc  = 0.0_r8
          currentCohort%gpp_acc  = 0.0_r8
          currentCohort%resp_m_acc = 0.0_r8

          ! BOC...update tree 'hydraulic geometry'
          ! (size --> heights of elements --> hydraulic path lengths -->
          ! maximum node-to-node conductances)
          if( (hlm_use_planthydro.eq.itrue) .and. do_growthrecruiteffects) then
             call UpdateSizeDepPlantHydProps(currentSite,currentCohort, bc_in)
             call UpdateSizeDepPlantHydStates(currentSite,currentCohort)
          end if

          ! if we are in age-dependent mortality mode
          if (hlm_use_cohort_age_tracking .eq. itrue) then
             ! update cohort age
             currentCohort%coage = currentCohort%coage + hlm_freq_day
             if(currentCohort%coage < 0.0_r8)then
                write(fates_log(),*) 'negative cohort age?',currentCohort%coage
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if

             ! update cohort age class and age x pft class
             call coagetype_class_index(currentCohort%coage, currentCohort%pft, &
                  currentCohort%coage_class,currentCohort%coage_by_pft_class)
          end if
        
          currentCohort => currentCohort%taller
       end do

       currentPatch => currentPatch%younger
   end do
          
   ! We keep a record of the L2FRs of plants
   ! that are near the recruit size, for different
   ! pfts and canopy layer. We use this mean to
   ! set the L2FRs of newly recruited plants
   
   call UpdateRecruitL2FR(currentSite)

   ! Update history diagnostics related to Nutrients (if any)
   ! -----------------------------------------------------------------------------
   select case(hlm_parteh_mode)
   case (prt_cnp_flex_allom_hyp)
      call fates_hist%update_history_nutrflux(currentSite)
   end select
   
   ! When plants die, the water goes with them.  This effects

   ! the water balance.

    if( hlm_use_planthydro == itrue ) then
       currentPatch => currentSite%youngest_patch
       do while(associated(currentPatch))
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             call AccumulateMortalityWaterStorage(currentSite,currentCohort,&
                  -1.0_r8 * currentCohort%dndt * hlm_freq_day)
             currentCohort => currentCohort%taller
          end do
          currentPatch => currentPatch%older
      end do
    end if


    ! With growth and mortality rates now calculated we can determine the seed rain
    ! fluxes. However, because this is potentially a cross-patch mixing model
    ! we will calculate this as a group
    
    call SeedUpdate(currentSite)

    ! Calculate all other litter fluxes
    ! -----------------------------------------------------------------------------------

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       
       call GenerateDamageAndLitterFluxes( currentSite, currentPatch, bc_in)

       call PreDisturbanceLitterFluxes( currentSite, currentPatch, bc_in)

       call PreDisturbanceIntegrateLitter(currentPatch )

       currentPatch => currentPatch%older
    enddo


    ! RGK: This call is unecessary for CLM coupling. I believe we
    ! can remove it completely if/when this call is added in ELM to 
    ! subroutine UpdateLitterFluxes(this,bounds_clump) in elmfates_interfaceMod.F90

    call FluxIntoLitterPools(currentsite, bc_in, bc_out)


    ! Update cohort number.
    ! This needs to happen after the CWD_input and seed_input calculations as they
    ! assume the pre-mortality currentCohort%n.

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))
          currentCohort%n = max(0._r8,currentCohort%n + currentCohort%dndt * hlm_freq_day )
          currentCohort%sym_nfix_daily = 0._r8
          currentCohort => currentCohort%taller
       enddo
       currentPatch => currentPatch%older
   enddo


   return
  end subroutine ed_integrate_state_variables

  !-------------------------------------------------------------------------------!
  subroutine ed_update_site( currentSite, bc_in, bc_out, is_restarting )
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
    type(bc_in_type)   , intent(in)       :: bc_in
    type(bc_out_type)  , intent(inout)    :: bc_out
    logical,intent(in)                    :: is_restarting ! is this called during restart read?
    !
    ! !LOCAL VARIABLES:
    type (fates_patch_type) , pointer :: currentPatch
    !-----------------------------------------------------------------------

    ! check patch order (set second argument to true)
    if (debug) then
       call set_patchno(currentSite,.true.,1)
    end if
    
    if(hlm_use_sp.eq.ifalse .and. (.not.is_restarting))then
      call canopy_spread(currentSite)
    end if

    call TotalBalanceCheck(currentSite,6)

    if(hlm_use_sp.eq.ifalse .and. (.not.is_restarting) )then
       call canopy_structure(currentSite, bc_in)
    endif

    call TotalBalanceCheck(currentSite,final_check_id)

    ! Update recruit L2FRs based on new canopy position
    call SetRecruitL2FR(currentSite)
    
    currentSite%area_by_age(:) = 0._r8
    
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))

       if(.not.is_restarting)then
          call terminate_cohorts(currentSite, currentPatch, 1, 11, bc_in)
          call terminate_cohorts(currentSite, currentPatch, 2, 11, bc_in)
       end if

       ! This cohort count is used in the photosynthesis loop
       call currentPatch%CountCohorts()
       
       ! Update the total area of by patch age class array 
       currentSite%area_by_age(currentPatch%age_class) = &
            currentSite%area_by_age(currentPatch%age_class) + currentPatch%area
       
       currentPatch => currentPatch%younger
       
    enddo

    ! Check to see if the time integrated fluxes match the state
    ! Dont call this if we are restarting, it will double count the flux
    if(.not.is_restarting)then
       call CheckIntegratedMassPools(currentSite)
    end if
    
    ! The HLMs need to know about nutrient demand, and/or
    ! root mass and affinities
    call PrepNutrientAquisitionBCs(currentSite,bc_in,bc_out)

    ! The HLM methane module needs information about
    ! rooting mass, distributions, respiration rates and NPP
    call PrepCH4BCs(currentSite,bc_in,bc_out)


    ! FIX(RF,032414). This needs to be monthly, not annual
    ! If this is the second to last day of the year, then perform trimming
    if( hlm_day_of_year == hlm_days_per_year-1 .and. (.not.is_restarting)) then
       if(hlm_use_sp.eq.ifalse)then
          call trim_canopy(currentSite)
     endif
    endif

  end subroutine ed_update_site

  !-------------------------------------------------------------------------------!

  subroutine TotalBalanceCheck (currentSite, call_index )

    !
    ! !DESCRIPTION:
    ! This routine looks at the mass flux in and out of the FATES and compares it to
    ! the change in total stocks (states).
    ! Fluxes in are NPP. Fluxes out are decay of CWD and litter into SOM pools.
    !
    ! !ARGUMENTS:
    type(ed_site_type) , intent(inout) :: currentSite
    integer            , intent(in)    :: call_index
    !
    ! !LOCAL VARIABLES:
    type(site_massbal_type),pointer :: site_mass
    real(r8) :: biomass_stock   ! total biomass   in Kg/site
    real(r8) :: litter_stock    ! total litter    in Kg/site
    real(r8) :: seed_stock      ! total seed mass in Kg/site
    real(r8) :: total_stock     ! total ED carbon in Kg/site
    real(r8) :: change_in_stock ! Change since last time we set ed_allsites_inst%old_stock in this routine.  KgC/site
    real(r8) :: error           ! How much carbon did we gain or lose (should be zero!)
    real(r8) :: error_frac      ! Error as a fraction of total biomass
    real(r8) :: net_flux        ! Difference between recorded fluxes in and out. KgC/site
    real(r8) :: flux_in         ! mass flux into fates control volume
    real(r8) :: flux_out        ! mass flux out of fates control volume
    real(r8) :: leaf_m          ! Mass in leaf tissues kg
    real(r8) :: fnrt_m          ! "" fine root
    real(r8) :: sapw_m          ! "" sapwood
    real(r8) :: store_m         ! "" storage
    real(r8) :: struct_m        ! "" structure
    real(r8) :: repro_m         ! "" reproduction

    integer  :: el              ! loop counter for element types

    ! nb. There is no time associated with these variables
    ! because this routine can be called between any two
    ! arbitrary points in code, even if no time has passed.
    ! Also, the carbon pools are per site/gridcell, so that
    ! we can account for the changing areas of patches.

    type(fates_patch_type)  , pointer :: currentPatch
    type(fates_cohort_type) , pointer :: currentCohort
    type(litter_type), pointer     :: litt
    logical, parameter :: print_cohorts = .true.   ! Set to true if you want
                                                    ! to print cohort data
                                                    ! upon fail (lots of text)
    !-----------------------------------------------------------------------

  if(hlm_use_sp.eq.ifalse)then

    change_in_stock = 0.0_r8


    ! Loop through the number of elements in the system

    do el = 1, num_elements

       site_mass => currentSite%mass_balance(el)

       call SiteMassStock(currentSite,el,total_stock,biomass_stock,litter_stock,seed_stock)

       change_in_stock = total_stock - site_mass%old_stock

       flux_in  = site_mass%seed_in + &
                  site_mass%net_root_uptake + &
                  site_mass%gpp_acc + &
                  site_mass%flux_generic_in + &
                  site_mass%patch_resize_err

       flux_out = sum(site_mass%wood_product_harvest(:)) + &
                  sum(site_mass%wood_product_landusechange(:)) + &
                  site_mass%burn_flux_to_atm + &
                  site_mass%seed_out + &
                  site_mass%flux_generic_out + &
                  site_mass%frag_out + &
                  site_mass%aresp_acc + &
                  site_mass%herbivory_flux_out

       net_flux        = flux_in - flux_out
       error           = abs(net_flux - change_in_stock)


       if(change_in_stock>0.0)then
          error_frac      = error/abs(total_stock)
       else
          error_frac      = 0.0_r8
       end if

       if ( error_frac > 10e-6_r8 .or. (error /= error) ) then
          write(fates_log(),*) 'mass balance error detected'
          write(fates_log(),*) 'element type (see PRTGenericMod.F90): ',element_list(el)
          write(fates_log(),*) 'error fraction relative to biomass stock: ',error_frac
          write(fates_log(),*) 'absolut error (flux in - change): ',net_flux - change_in_stock
          write(fates_log(),*) 'call index: ',call_index
          write(fates_log(),*) 'Element index (PARTEH global):',element_list(el)
          write(fates_log(),*) 'net: ',net_flux
          write(fates_log(),*) 'dstock: ',change_in_stock
          write(fates_log(),*) 'seed_in: ',site_mass%seed_in
          write(fates_log(),*) 'net_root_uptake: ',site_mass%net_root_uptake
          write(fates_log(),*) 'gpp_acc: ',site_mass%gpp_acc
          write(fates_log(),*) 'flux_generic_in: ',site_mass%flux_generic_in
          write(fates_log(),*) 'wood_product_harvest: ',site_mass%wood_product_harvest(:)
          write(fates_log(),*) 'wood_product_landusechange: ',site_mass%wood_product_landusechange(:)
          write(fates_log(),*) 'error from patch resizing: ',site_mass%patch_resize_err
          write(fates_log(),*) 'burn_flux_to_atm: ',site_mass%burn_flux_to_atm
          write(fates_log(),*) 'seed_out: ',site_mass%seed_out
          write(fates_log(),*) 'flux_generic_out: ',site_mass%flux_generic_out
          write(fates_log(),*) 'frag_out: ',site_mass%frag_out
          write(fates_log(),*) 'aresp_acc: ',site_mass%aresp_acc
          write(fates_log(),*) 'herbivory_flux_out: ',site_mass%herbivory_flux_out
          write(fates_log(),*) 'error=net_flux-dstock:', error
          write(fates_log(),*) 'biomass', biomass_stock
          write(fates_log(),*) 'litter',litter_stock
          write(fates_log(),*) 'seeds',seed_stock
          write(fates_log(),*) 'total stock', total_stock
          write(fates_log(),*) 'previous total',site_mass%old_stock
          write(fates_log(),*) 'lat lon',currentSite%lat,currentSite%lon

          ! If this is the first day of simulation, carbon balance reports but does not end the run
!          if(( hlm_current_year*10000 + hlm_current_month*100 + hlm_current_day).ne.hlm_reference_date) then

             currentPatch => currentSite%oldest_patch
             do while(associated(currentPatch))
                litt => currentPatch%litter(el)
                write(fates_log(),*) '---------------------------------------'
                write(fates_log(),*) 'patch area: ',currentPatch%area
                write(fates_log(),*) 'AG CWD: ', sum(litt%ag_cwd)
                write(fates_log(),*) 'BG CWD (by layer): ', sum(litt%bg_cwd,dim=1)
                write(fates_log(),*) 'leaf litter:',sum(litt%leaf_fines)
                write(fates_log(),*) 'root litter (by layer): ',sum(litt%root_fines,dim=1)
                write(fates_log(),*) 'land_use_label: ',currentPatch%land_use_label
                write(fates_log(),*) 'use_this_pft: ', currentSite%use_this_pft(:)
                if(print_cohorts)then
                    write(fates_log(),*) '---- Biomass by cohort and organ -----'
                    currentCohort => currentPatch%tallest
                    do while(associated(currentCohort))
                        write(fates_log(),*) 'pft: ',currentCohort%pft
                        write(fates_log(),*) 'dbh: ',currentCohort%dbh
                        leaf_m   = currentCohort%prt%GetState(leaf_organ,element_list(el))
                        struct_m = currentCohort%prt%GetState(struct_organ,element_list(el))
                        store_m  = currentCohort%prt%GetState(store_organ,element_list(el))
                        fnrt_m   = currentCohort%prt%GetState(fnrt_organ,element_list(el))
                        repro_m  = currentCohort%prt%GetState(repro_organ,element_list(el))
                        sapw_m   = currentCohort%prt%GetState(sapw_organ,element_list(el))
                        write(fates_log(),*) 'leaf: ',leaf_m,' structure: ',struct_m,' store: ',store_m
                        write(fates_log(),*) 'fineroot: ',fnrt_m,' repro: ',repro_m,' sapwood: ',sapw_m
                        write(fates_log(),*) 'num plant: ',currentCohort%n
                        write(fates_log(),*) 'resp excess: ',currentCohort%resp_excess_hold*currentCohort%n

                        if(element_list(el).eq.nitrogen_element) then
                           write(fates_log(),*) 'NH4 uptake: ',currentCohort%daily_nh4_uptake*currentCohort%n
                           write(fates_log(),*) 'NO3 uptake: ',currentCohort%daily_no3_uptake*currentCohort%n
                           write(fates_log(),*) 'N efflux: ',currentCohort%daily_n_efflux*currentCohort%n
                           write(fates_log(),*) 'N fixation: ',currentCohort%sym_nfix_daily*currentCohort%n
                        elseif(element_list(el).eq.phosphorus_element) then
                           write(fates_log(),*) 'P uptake: ',currentCohort%daily_p_gain*currentCohort%n
                           write(fates_log(),*) 'P efflux: ',currentCohort%daily_p_efflux*currentCohort%n
                        elseif(element_list(el).eq.carbon12_element) then
                           write(fates_log(),*) 'C efflux: ',currentCohort%daily_c_efflux*currentCohort%n
                        end if


                        currentCohort => currentCohort%shorter
                    enddo !end cohort loop
                end if
                currentPatch => currentPatch%younger
             enddo !end patch loop
             write(fates_log(),*) 'aborting on date:',hlm_current_year,hlm_current_month,hlm_current_day
             call endrun(msg=errMsg(sourcefile, __LINE__))
         !end if

      endif

      ! This is the last check of the sequence, where we update our total
      ! error check and the final fates stock
      if(call_index == final_check_id) then
          site_mass%old_stock = total_stock
          site_mass%err_fates = net_flux - change_in_stock
      end if

   end do
  end if ! not SP mode
  end subroutine TotalBalanceCheck

  ! =====================================================================================

  subroutine bypass_dynamics(currentSite, bc_out)

    ! ----------------------------------------------------------------------------------
    ! If dynamics are bypassed, various fluxes, rates and flags need to be set
    ! to trivial values.
    ! WARNING: Turning off things like dynamics is experimental. The setting of
    ! variables to trivial values may not be complete, use at your own risk.
    ! ----------------------------------------------------------------------------------

    ! Arguments
    type(ed_site_type)      , intent(inout) :: currentSite
    type(bc_out_type)       , intent(inout) :: bc_out
    
    ! Locals
    type(fates_patch_type), pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    bc_out%gpp_site = 0._r8
    bc_out%ar_site = 0._r8
    
    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))

          currentCohort%isnew=.false.

          currentCohort%resp_g_acc_hold = prt_params%grperc(currentCohort%pft) * &
               max(0._r8,(currentCohort%gpp_acc - currentCohort%resp_m_acc))*real(hlm_days_per_year,r8)
          
          currentCohort%npp_acc_hold  = currentCohort%npp_acc  * real(hlm_days_per_year,r8)
          currentCohort%gpp_acc_hold  = currentCohort%gpp_acc  * real(hlm_days_per_year,r8)
          currentCohort%resp_m_acc_hold = currentCohort%resp_m_acc * real(hlm_days_per_year,r8)
          
          currentCohort%resp_excess_hold = 0._r8
          currentCohort%npp_acc  = 0.0_r8
          currentCohort%gpp_acc  = 0.0_r8
          currentCohort%resp_m_acc = 0.0_r8

          ! No need to set the "net_art" terms to zero
          ! they are zeroed at the beginning of the daily step
          ! If DailyPRT, maintenance, and phenology are not called
          ! then these should stay zero.

          currentCohort%bmort = 0.0_r8
          currentCohort%hmort = 0.0_r8
          currentCohort%cmort = 0.0_r8
          currentCohort%frmort = 0.0_r8
          currentCohort%smort = 0.0_r8
          currentCohort%asmort = 0.0_r8
          currentCohort%dgmort = 0.0_r8

          currentCohort%dndt      = 0.0_r8
          currentCohort%dhdt      = 0.0_r8
          currentCohort%ddbhdt    = 0.0_r8

          ! Shouldn't need to zero any nutrient fluxes
          ! as they should just be zero, no uptake
          ! in ST3 mode.

          ! Passing 
          bc_out%gpp_site = bc_out%gpp_site + currentCohort%gpp_acc_hold * &
               AREA_INV * currentCohort%n / hlm_days_per_year / sec_per_day
          bc_out%ar_site = bc_out%ar_site + (currentCohort%resp_m_acc_hold + &
               currentCohort%resp_g_acc_hold + currentCohort%resp_excess_hold) * & 
               AREA_INV * currentCohort%n / hlm_days_per_year / sec_per_day

          currentCohort => currentCohort%taller
       enddo
       currentPatch => currentPatch%older
    enddo

 end subroutine bypass_dynamics

end module EDMainMod
