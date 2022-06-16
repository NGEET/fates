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
  use EDCohortDynamicsMod      , only : sort_cohorts
  use EDCohortDynamicsMod      , only : count_cohorts
  use EDCohortDynamicsMod      , only : EvaluateAndCorrectDBH
  use EDPatchDynamicsMod       , only : disturbance_rates
  use EDPatchDynamicsMod       , only : fuse_patches
  use EDPatchDynamicsMod       , only : spawn_patches
  use EDPatchDynamicsMod       , only : terminate_patches
  use EDPhysiologyMod          , only : phenology
  use EDPhysiologyMod          , only : satellite_phenology
  use EDPhysiologyMod          , only : recruitment
  use EDPhysiologyMod          , only : trim_canopy
  use EDPhysiologyMod          , only : SeedIn
  use EDPhysiologyMod          , only : ZeroAllocationRates
  use EDPhysiologyMod          , only : ZeroLitterFluxes
  use EDPhysiologyMod          , only : PreDisturbanceLitterFluxes
  use EDPhysiologyMod          , only : PreDisturbanceIntegrateLitter
  use FatesSoilBGCFluxMod      , only : FluxIntoLitterPools
  use EDCohortDynamicsMod      , only : UpdateCohortBioPhysRates
  use FatesSoilBGCFluxMod      , only : PrepNutrientAquisitionBCs
  use FatesSoilBGCFluxMod      , only : PrepCH4BCs
  use SFMainMod                , only : fire_model
  use FatesSizeAgeTypeIndicesMod, only : get_age_class_index
  use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index
  use FatesLitterMod           , only : litter_type
  use FatesLitterMod           , only : ncwd
  use EDtypesMod               , only : ed_site_type
  use EDtypesMod               , only : ed_patch_type
  use EDtypesMod               , only : ed_cohort_type
  use EDTypesMod               , only : AREA
  use EDTypesMod               , only : site_massbal_type
  use PRTGenericMod            , only : num_elements
  use PRTGenericMod            , only : element_list
  use PRTGenericMod            , only : element_pos
  use EDTypesMod               , only : phen_dstat_moiston
  use EDTypesMod               , only : phen_dstat_timeon
  use FatesConstantsMod        , only : itrue,ifalse
  use FatesConstantsMod        , only : primaryforest, secondaryforest
  use FatesConstantsMod        , only : nearzero
  use FatesConstantsMod        , only : m2_per_ha
  use FatesConstantsMod        , only : sec_per_day
  use FatesPlantHydraulicsMod  , only : do_growthrecruiteffects
  use FatesPlantHydraulicsMod  , only : UpdateSizeDepPlantHydProps
  use FatesPlantHydraulicsMod  , only : UpdateSizeDepPlantHydStates
  use FatesPlantHydraulicsMod  , only : InitPlantHydStates
  use FatesPlantHydraulicsMod  , only : UpdateSizeDepRhizHydProps
  use FatesPlantHydraulicsMod  , only : AccumulateMortalityWaterStorage
  use FatesAllometryMod        , only : h_allom,tree_sai,tree_lai
  use EDLoggingMortalityMod    , only : IsItLoggingTime
  use EDLoggingMortalityMod    , only : get_harvestable_carbon
  use EDPatchDynamicsMod       , only : get_frac_site_primary
  use FatesGlobals             , only : endrun => fates_endrun
  use ChecksBalancesMod        , only : SiteMassStock
  use EDMortalityFunctionsMod  , only : Mortality_Derivative
  use EDTypesMod               , only : AREA_INV
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
  use FatesHistoryInterfaceMod, only : ih_nh4uptake_si, ih_no3uptake_si, ih_puptake_si
  use FatesHistoryInterfaceMod, only : ih_nh4uptake_scpf, ih_no3uptake_scpf, ih_puptake_scpf
  use FatesHistoryInterfaceMod, only : fates_hist

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
    type(ed_patch_type), pointer :: currentPatch
    integer :: el                ! Loop counter for elements
    integer :: do_patch_dynamics ! for some modes, we turn off patch dynamics

    !-----------------------------------------------------------------------

    if ( hlm_masterproc==itrue ) write(fates_log(),'(A,I4,A,I2.2,A,I2.2)') 'FATES Dynamics: ',&
          hlm_current_year,'-',hlm_current_month,'-',hlm_current_day

    ! Consider moving this towards the end, because some of these
    ! are being integrated over the short time-step

    do el = 1,num_elements
       call currentSite%mass_balance(el)%ZeroMassBalFlux()
       call currentSite%flux_diags(el)%ZeroFluxDiags()
    end do

    ! zero dynamics (upfreq_in = 1) output history variables
    call fates_hist%zero_site_hvars(currentSite,upfreq_in=1)


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
       call fire_model(currentSite, bc_in)

       ! Calculate disturbance and mortality based on previous timestep vegetation.
       ! disturbance_rates calls logging mortality and other mortalities, Yi Xu
       call disturbance_rates(currentSite, bc_in)

       ! Integrate state variables from annual rates to daily timestep
       call ed_integrate_state_variables(currentSite, bc_in, bc_out )

    else
       ! ed_intergrate_state_variables is where the new cohort flag
       ! is set. This flag designates wether a cohort has
       ! experienced a day, and therefore has been populated with non-nonsense
       ! values.  If we aren't entering that sequence, we need to set the flag
       ! Make sure cohorts are marked as non-recruits

       call bypass_dynamics(currentSite)

    end if

    !******************************************************************************
    ! Reproduction, Recruitment and Cohort Dynamics : controls cohort organization
    !******************************************************************************

    if(hlm_use_ed_st3.eq.ifalse.and.hlm_use_sp.eq.ifalse) then
       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))

          ! adds small cohort of each PFT
          call recruitment(currentSite, currentPatch, bc_in)

          currentPatch => currentPatch%younger
       enddo

       call TotalBalanceCheck(currentSite,1)

       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))

          ! puts cohorts in right order
          call sort_cohorts(currentPatch)

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
       call terminate_patches(currentSite)
    end if

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
    use FatesInterfaceTypesMod, only : hlm_use_cohort_age_tracking
    use FatesConstantsMod, only : itrue
    ! !ARGUMENTS:

    type(ed_site_type)     , intent(inout) :: currentSite
    type(bc_in_type)        , intent(in)   :: bc_in
    type(bc_out_type)       , intent(inout)  :: bc_out

    !
    ! !LOCAL VARIABLES:
    type(site_massbal_type), pointer :: site_cmass
    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort

    integer  :: c                     ! Counter for litter size class
    integer  :: ft                    ! Counter for PFT
    integer  :: io_si                 ! global site index for history writing
    integer  :: iscpf                 ! index for the size-class x pft multiplexed bins
    integer  :: el                    ! Counter for element type (c,n,p,etc)
    real(r8) :: cohort_biomass_store  ! remembers the biomass in the cohort for balance checking
    real(r8) :: dbh_old               ! dbh of plant before daily PRT [cm]
    real(r8) :: hite_old              ! height of plant before daily PRT [m]
    logical  :: is_drought            ! logical for if the plant (site) is in a drought state
    real(r8) :: delta_dbh             ! correction for dbh
    real(r8) :: delta_hite            ! correction for hite

    real(r8) :: current_npp           ! place holder for calculating npp each year in prescribed physiology mode
    !-----------------------------------------------------------------------
    real(r8) :: frac_site_primary

    real(r8) :: harvestable_forest_c(hlm_num_lu_harvest_cats)
    real(r8) :: available_forest_c(hlm_num_lu_harvest_cats)
    integer  :: harvest_tag(hlm_num_lu_harvest_cats)


    call get_frac_site_primary(currentSite, frac_site_primary)

    ! Patch level biomass are required for C-based harvest
    call get_harvestable_carbon(currentSite, bc_in%site_area, bc_in%hlm_harvest_catnames, harvestable_forest_c, available_forest_c)

    ! Set a pointer to this sites carbon12 mass balance
    site_cmass => currentSite%mass_balance(element_pos(carbon12_element))

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
          call Mortality_Derivative( currentSite, currentCohort, bc_in, frac_site_primary, &
              harvestable_forest_c, available_forest_c, harvest_tag)

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
             currentCohort%resp_acc = 0._r8

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

          currentCohort%npp_acc_hold  = currentCohort%npp_acc  * real(hlm_days_per_year,r8)
          currentCohort%gpp_acc_hold  = currentCohort%gpp_acc  * real(hlm_days_per_year,r8)
          currentCohort%resp_acc_hold = currentCohort%resp_acc * real(hlm_days_per_year,r8)


          ! Conduct Maintenance Turnover (parteh)
          if(debug) call currentCohort%prt%CheckMassConservation(ft,3)
          if(any(currentSite%dstatus == [phen_dstat_moiston,phen_dstat_timeon])) then
             is_drought = .false.
          else
             is_drought = .true.
          end if
          call PRTMaintTurnover(currentCohort%prt,ft,is_drought)

          ! If the current diameter of a plant is somehow less than what is consistent
          ! with what is allometrically consistent with the stuctural biomass, then
          ! correct the dbh to match.

          call EvaluateAndCorrectDBH(currentCohort,delta_dbh,delta_hite)

          hite_old = currentCohort%hite
          dbh_old  = currentCohort%dbh

          ! -----------------------------------------------------------------------------
          ! Growth and Allocation (PARTEH)
          ! -----------------------------------------------------------------------------

          call currentCohort%prt%DailyPRT()


          ! Update the mass balance tracking for the daily nutrient uptake flux
          ! Then zero out the daily uptakes, they have been used
          ! -----------------------------------------------------------------------------

          if(hlm_parteh_mode .eq. prt_cnp_flex_allom_hyp ) then

             ! Mass balance for N uptake
             currentSite%mass_balance(element_pos(nitrogen_element))%net_root_uptake = &
                  currentSite%mass_balance(element_pos(nitrogen_element))%net_root_uptake + &
                  (currentCohort%daily_nh4_uptake+currentCohort%daily_no3_uptake- &
                  currentCohort%daily_n_efflux)*currentCohort%n

             ! Mass balance for P uptake
             currentSite%mass_balance(element_pos(phosphorus_element))%net_root_uptake = &
                  currentSite%mass_balance(element_pos(phosphorus_element))%net_root_uptake + &
                  (currentCohort%daily_p_uptake-currentCohort%daily_p_efflux)*currentCohort%n

             ! mass balance for C efflux (if any)
             currentSite%mass_balance(element_pos(carbon12_element))%net_root_uptake = &
                  currentSite%mass_balance(element_pos(carbon12_element))%net_root_uptake - &
                  currentCohort%daily_c_efflux*currentCohort%n

             ! size class index
             iscpf = currentCohort%size_by_pft_class

             ! Diagnostics for uptake, by size and pft, [kgX/ha/day]

             io_si  = currentSite%h_gid

             fates_hist%hvars(ih_nh4uptake_scpf)%r82d(io_si,iscpf) =           &
                  fates_hist%hvars(ih_nh4uptake_scpf)%r82d(io_si,iscpf) +      &
                  currentCohort%daily_nh4_uptake*currentCohort%n /             &
                  m2_per_ha / sec_per_day

             fates_hist%hvars(ih_no3uptake_scpf)%r82d(io_si,iscpf) =           &
                  fates_hist%hvars(ih_no3uptake_scpf)%r82d(io_si,iscpf) +      &
                  currentCohort%daily_no3_uptake*currentCohort%n /             &
                  m2_per_ha / sec_per_day

             fates_hist%hvars(ih_puptake_scpf)%r82d(io_si,iscpf) =             &
                  fates_hist%hvars(ih_puptake_scpf)%r82d(io_si,iscpf) +        &
                  currentCohort%daily_p_uptake*currentCohort%n /               &
                  m2_per_ha / sec_per_day

             fates_hist%hvars(ih_nh4uptake_si)%r81d(io_si) =                   &
                  fates_hist%hvars(ih_nh4uptake_si)%r81d(io_si)  +             &
                  currentCohort%daily_nh4_uptake*currentCohort%n /             &
                  m2_per_ha / sec_per_day

             fates_hist%hvars(ih_no3uptake_si)%r81d(io_si) =                   &
                  fates_hist%hvars(ih_no3uptake_si)%r81d(io_si)  +             &
                  currentCohort%daily_no3_uptake*currentCohort%n /             &
                  m2_per_ha / sec_per_day

             fates_hist%hvars(ih_puptake_si)%r81d(io_si) =                     &
                  fates_hist%hvars(ih_puptake_si)%r81d(io_si)  +               &
                  currentCohort%daily_p_uptake*currentCohort%n /               &
                  m2_per_ha / sec_per_day


             ! Diagnostics on efflux, size and pft [kgX/ha/day]
             currentSite%flux_diags(element_pos(nitrogen_element))%nutrient_efflux_scpf(iscpf) = &
                  currentSite%flux_diags(element_pos(nitrogen_element))%nutrient_efflux_scpf(iscpf) + &
                  currentCohort%daily_n_efflux*currentCohort%n

             currentSite%flux_diags(element_pos(phosphorus_element))%nutrient_efflux_scpf(iscpf) = &
                  currentSite%flux_diags(element_pos(phosphorus_element))%nutrient_efflux_scpf(iscpf) + &
                  currentCohort%daily_p_efflux*currentCohort%n

             currentSite%flux_diags(element_pos(carbon12_element))%nutrient_efflux_scpf(iscpf) = &
                  currentSite%flux_diags(element_pos(carbon12_element))%nutrient_efflux_scpf(iscpf) + &
                  currentCohort%daily_c_efflux*currentCohort%n

             ! Diagnostics on plant nutrient need
             currentSite%flux_diags(element_pos(nitrogen_element))%nutrient_need_scpf(iscpf) = &
                  currentSite%flux_diags(element_pos(nitrogen_element))%nutrient_need_scpf(iscpf) + &
                  currentCohort%daily_n_need*currentCohort%n

             currentSite%flux_diags(element_pos(phosphorus_element))%nutrient_need_scpf(iscpf) = &
                  currentSite%flux_diags(element_pos(phosphorus_element))%nutrient_need_scpf(iscpf) + &
                  currentCohort%daily_p_need*currentCohort%n

          end if

          ! And simultaneously add the input fluxes to mass balance accounting
          site_cmass%gpp_acc   = site_cmass%gpp_acc + &
                currentCohort%gpp_acc * currentCohort%n
          site_cmass%aresp_acc = site_cmass%aresp_acc + &
                currentCohort%resp_acc * currentCohort%n

          call currentCohort%prt%CheckMassConservation(ft,5)

          ! Update the leaf biophysical rates based on proportion of leaf
          ! mass in the different leaf age classes. Following growth
          ! and turnover, these proportions won't change again. This
          ! routine is also called following fusion
          call UpdateCohortBioPhysRates(currentCohort)

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
             call UpdateSizeDepPlantHydProps(currentSite,currentCohort, bc_in)
             call UpdateSizeDepPlantHydStates(currentSite,currentCohort)
          end if

          ! if we are in age-dependent mortality mode
          if (hlm_use_cohort_age_tracking .eq. itrue) then
             ! update cohort age
             currentCohort%coage = currentCohort%coage + hlm_freq_day
             if(currentCohort%coage < 0.0_r8)then
                write(fates_log(),*) 'negative cohort age?',currentCohort%coage
             end if

             ! update cohort age class and age x pft class
             call coagetype_class_index(currentCohort%coage, currentCohort%pft, &
                  currentCohort%coage_class,currentCohort%coage_by_pft_class)
          end if


          currentCohort => currentCohort%taller
      end do

       currentPatch => currentPatch%older
   end do


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

    call SeedIn(currentSite,bc_in)

    ! Calculate all other litter fluxes
    ! -----------------------------------------------------------------------------------

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))

       call PreDisturbanceLitterFluxes( currentSite, currentPatch, bc_in)


       call PreDisturbanceIntegrateLitter(currentPatch )

       currentPatch => currentPatch%older
    enddo


    ! Before we start messing with the patch areas, and before we start removing
    ! trees, this is a good time to pass fragmentation litter fluxes and
    ! plant-to-soil fluxes (such as efflux and fixation fluxes)

    call FluxIntoLitterPools(currentsite, bc_in, bc_out)


    ! Update cohort number.
    ! This needs to happen after the CWD_input and seed_input calculations as they
    ! assume the pre-mortality currentCohort%n.

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))
          currentCohort%n = max(0._r8,currentCohort%n + currentCohort%dndt * hlm_freq_day )
          currentCohort => currentCohort%taller
       enddo
       currentPatch => currentPatch%older
   enddo


   return
  end subroutine ed_integrate_state_variables

  !-------------------------------------------------------------------------------!
  subroutine ed_update_site( currentSite, bc_in, bc_out )
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
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type) , pointer :: currentPatch
    !-----------------------------------------------------------------------
    if(hlm_use_sp.eq.ifalse)then
      call canopy_spread(currentSite)
    end if

    call TotalBalanceCheck(currentSite,6)

    if(hlm_use_sp.eq.ifalse)then
       call canopy_structure(currentSite, bc_in)
    endif

    call TotalBalanceCheck(currentSite,final_check_id)

    currentSite%area_by_age(:) = 0._r8
    
    currentPatch => currentSite%oldest_patch
    do while(associated(currentPatch))

        ! Is termination really needed here?
        ! Canopy_structure just called it several times! (rgk)
        call terminate_cohorts(currentSite, currentPatch, 1, 11, bc_in)
        call terminate_cohorts(currentSite, currentPatch, 2, 11, bc_in)

        ! This cohort count is used in the photosynthesis loop
        call count_cohorts(currentPatch)

        ! Update the total area of by patch age class array 
        currentSite%area_by_age(currentPatch%age_class) = &
             currentSite%area_by_age(currentPatch%age_class) + currentPatch%area
        
        currentPatch => currentPatch%younger

    enddo

    ! The HLMs need to know about nutrient demand, and/or
    ! root mass and affinities
    call PrepNutrientAquisitionBCs(currentSite,bc_in,bc_out)

    ! The HLM methane module needs information about
    ! rooting mass, distributions, respiration rates and NPP
    call PrepCH4BCs(currentSite,bc_in,bc_out)


    ! FIX(RF,032414). This needs to be monthly, not annual
    ! If this is the second to last day of the year, then perform trimming
    if( hlm_day_of_year == hlm_days_per_year-1) then

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

    type(ed_patch_type)  , pointer :: currentPatch
    type(ed_cohort_type) , pointer :: currentCohort
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

       flux_out = site_mass%wood_product + &
                  site_mass%burn_flux_to_atm + &
                  site_mass%seed_out + &
                  site_mass%flux_generic_out + &
                  site_mass%frag_out + &
                  site_mass%aresp_acc

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
          write(fates_log(),*) 'wood_product: ',site_mass%wood_product
          write(fates_log(),*) 'error from patch resizing: ',site_mass%patch_resize_err
          write(fates_log(),*) 'burn_flux_to_atm: ',site_mass%burn_flux_to_atm
          write(fates_log(),*) 'seed_out: ',site_mass%seed_out
          write(fates_log(),*) 'flux_generic_out: ',site_mass%flux_generic_out
          write(fates_log(),*) 'frag_out: ',site_mass%frag_out
          write(fates_log(),*) 'aresp_acc: ',site_mass%aresp_acc
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
                write(fates_log(),*) 'anthro_disturbance_label: ',currentPatch%anthro_disturbance_label
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
                        write(fates_log(),*) 'resp m def: ',currentCohort%resp_m_def*currentCohort%n

                        if(element_list(el).eq.nitrogen_element) then
                           write(fates_log(),*) 'NH4 uptake: ',currentCohort%daily_nh4_uptake*currentCohort%n
                           write(fates_log(),*) 'NO3 uptake: ',currentCohort%daily_no3_uptake*currentCohort%n
                           write(fates_log(),*) 'N efflux: ',currentCohort%daily_n_efflux*currentCohort%n
                        elseif(element_list(el).eq.phosphorus_element) then
                           write(fates_log(),*) 'P uptake: ',currentCohort%daily_p_uptake*currentCohort%n
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
          currentCohort%smort = 0.0_r8
          currentCohort%asmort = 0.0_r8

          currentCohort%dndt      = 0.0_r8
          currentCohort%dhdt      = 0.0_r8
          currentCohort%ddbhdt    = 0.0_r8

          ! Shouldn't need to zero any nutrient fluxes
          ! as they should just be zero, no uptake
          ! in ST3 mode.

          currentCohort => currentCohort%taller
       enddo
       currentPatch => currentPatch%older
    enddo

 end subroutine bypass_dynamics

end module EDMainMod
