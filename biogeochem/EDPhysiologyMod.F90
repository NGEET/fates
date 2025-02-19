module EDPhysiologyMod

#include "shr_assert.h"

  ! ============================================================================
  ! Miscellaneous physiology routines from ED.
  ! ============================================================================

  use FatesGlobals, only         : fates_log
  use FatesInterfaceTypesMod, only    : hlm_days_per_year
  use FatesInterfaceTypesMod, only    : hlm_model_day
  use FatesInterfaceTypesMod, only    : hlm_freq_day
  use FatesInterfaceTypesMod, only    : hlm_day_of_year
  use FatesInterfaceTypesMod, only    : numpft
  use FatesInterfaceTypesMod, only    : nleafage
  use FatesInterfaceTypesMod, only    : nlevdamage
  use FatesInterfaceTypesMod, only    : hlm_use_planthydro
  use FatesInterfaceTypesMod, only    : hlm_parteh_mode
  use FatesInterfaceTypesMod, only    : hlm_use_fixed_biogeog
  use FatesInterfaceTypesMod, only    : hlm_use_nocomp
  use EDParamsMod           , only    : crop_lu_pft_vector     
  use FatesInterfaceTypesMod, only    : hlm_nitrogen_spec
  use FatesInterfaceTypesMod, only    : hlm_phosphorus_spec
  use FatesInterfaceTypesMod, only    : hlm_use_tree_damage
  use FatesInterfaceTypesMod, only : hlm_use_ed_prescribed_phys
  use FatesConstantsMod, only    : r8 => fates_r8
  use FatesConstantsMod, only    : nearzero
  use FatesConstantsMod, only    : sec_per_day
  use FatesConstantsMod, only    : default_regeneration
  use FatesConstantsMod, only    : TRS_regeneration
  use FatesConstantsMod, only    : TRS_no_seedling_dyn
  use FatesConstantsMod, only    : min_max_dbh_for_trees
  use FatesConstantsMod, only    : megajoules_per_joule
  use FatesConstantsMod, only    : mpa_per_mm_suction
  use FatesConstantsMod, only    : g_per_kg
  use FatesConstantsMod, only    : ndays_per_year
  use FatesConstantsMod, only    : nocomp_bareground
  use FatesConstantsMod, only    : nocomp_bareground_land
  use FatesConstantsMod, only    : is_crop
  use FatesConstantsMod, only    : area_error_2
  use EDPftvarcon      , only    : EDPftvarcon_inst
  use PRTParametersMod , only    : prt_params
  use EDPftvarcon      , only    : GetDecompyFrac
  use FatesInterfaceTypesMod, only    : bc_in_type
  use FatesInterfaceTypesMod, only    : bc_out_type
  use EDCohortDynamicsMod , only : create_cohort
  use EDCohortDynamicsMod , only : InitPRTObject
  use FatesAllometryMod   , only : tree_lai_sai
  use FatesAllometryMod   , only : leafc_from_treelai
  use FatesAllometryMod   , only : decay_coeff_vcmax
  use FatesLitterMod      , only : litter_type
  use EDTypesMod          , only : site_massbal_type
  use EDTypesMod          , only : numlevsoil_max
  use EDTypesMod          , only : numWaterMem
  use FatesFuelClassesMod , only : fuel_classes
  use EDTypesMod          , only : elem_diag_type
  use EDParamsMod         , only : dinc_vai, dlower_vai
  use EDTypesMod          , only : area_inv
  use EDTypesMod          , only : AREA
  use FatesLitterMod      , only : ncwd
  use FatesLitterMod      , only : ndcmpy
  use FatesLitterMod      , only : ilabile
  use FatesLitterMod      , only : ilignin
  use FatesLitterMod      , only : icellulose
  use FatesLitterMod      , only : adjust_SF_CWD_frac
  use EDParamsMod         , only : nclmax
  use EDTypesMod          , only : AREA,AREA_INV
  use FatesConstantsMod   , only : leaves_shedding
  use FatesConstantsMod   , only : ihard_stress_decid
  use FatesConstantsMod   , only : isemi_stress_decid
  use EDParamsMod         , only : nlevleaf
  use EDTypesMod          , only : num_vegtemp_mem
  use EDParamsMod         , only : maxpft
  use EDTypesMod          , only : ed_site_type
  use FatesPatchMod,        only : fates_patch_type
  use FatesCohortMod,       only : fates_cohort_type
  use FatesConstantsMod   , only : leaves_on
  use FatesConstantsMod   , only : leaves_off
  use EDTypesMod          , only : min_n_safemath
  use PRTGenericMod       , only : num_elements
  use PRTGenericMod       , only : element_list
  use PRTGenericMod       , only : element_pos
  use EDTypesMod          , only : site_fluxdiags_type
  use EDTypesMod          , only : phen_cstat_nevercold
  use EDTypesMod          , only : phen_cstat_iscold
  use EDTypesMod          , only : phen_cstat_notcold
  use EDTypesMod          , only : phen_dstat_timeoff
  use EDTypesMod          , only : phen_dstat_moistoff
  use EDTypesMod          , only : phen_dstat_moiston
  use EDTypesMod          , only : phen_dstat_timeon
  use EDTypesMod          , only : phen_dstat_pshed
  use EDTypesMod          , only : phen_dstat_pshed
  use EDTypesMod          , only : init_recruit_trim
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use EDParamsMod           , only : fates_mortality_disturbance_fraction
  use EDParamsMod           , only : q10_mr
  use EDParamsMod           , only : q10_froz
  use EDParamsMod           , only : logging_export_frac
  use EDParamsMod           , only : regeneration_model
  use EDParamsMod           , only : sdlng_mort_par_timescale
  use FatesPlantHydraulicsMod  , only : AccumulateMortalityWaterStorage
  use FatesConstantsMod     , only : itrue,ifalse
  use FatesConstantsMod     , only : area_error_3
  use FatesConstantsMod     , only : years_per_day
  use FatesAllometryMod  , only : h_allom
  use FatesAllometryMod  , only : h2d_allom
  use FatesAllometryMod  , only : bagw_allom
  use FatesAllometryMod  , only : bsap_allom
  use FatesAllometryMod  , only : bleaf
  use FatesAllometryMod  , only : bfineroot
  use FatesAllometryMod  , only : bdead_allom
  use FatesAllometryMod  , only : bstore_allom
  use FatesAllometryMod  , only : bbgw_allom
  use FatesAllometryMod  , only : carea_allom
  use FatesAllometryMod  , only : CheckIntegratedAllometries
  use FatesAllometryMod, only : set_root_fraction
  use PRTGenericMod, only : prt_carbon_allom_hyp
  use PRTGenericMod, only : prt_cnp_flex_allom_hyp
  use PRTGenericMod, only : prt_vartypes
  use PRTGenericMod, only : leaf_organ
  use PRTGenericMod, only : sapw_organ, struct_organ
  use PRTGenericMod, only : carbon12_element
  use PRTGenericMod, only : nitrogen_element
  use PRTGenericMod, only : phosphorus_element
  use PRTGenericMod, only : leaf_organ
  use PRTGenericMod, only : fnrt_organ
  use PRTGenericMod, only : sapw_organ
  use PRTGenericMod, only : store_organ
  use PRTGenericMod, only : repro_organ
  use PRTGenericMod, only : struct_organ
  use PRTGenericMod, only : SetState
  use PRTLossFluxesMod, only  : PRTPhenologyFlush
  use PRTLossFluxesMod, only  : PRTDeciduousTurnover
  use PRTLossFluxesMod, only  : PRTReproRelease
  use PRTLossFluxesMod, only  : PRTDamageLosses
  use PRTGenericMod, only     : StorageNutrientTarget
  use DamageMainMod, only     : damage_time
  use DamageMainMod, only     : GetCrownReduction
  use DamageMainMod, only     : GetDamageFrac
  use SFParamsMod, only       : SF_val_CWD_frac
  use FatesParameterDerivedMod, only : param_derived
  use FatesPlantHydraulicsMod, only : InitHydrCohort
  use PRTInitParamsFatesMod, only : NewRecruitTotalStoichiometry
  use FatesInterfaceTypesMod    , only : hlm_use_luh

  implicit none
  private

  public :: trim_canopy
  public :: phenology
  public :: satellite_phenology
  public :: assign_cohort_SP_properties
  public :: calculate_SP_properties
  public :: recruitment
  public :: ZeroLitterFluxes

  public :: ZeroAllocationRates
  public :: PreDisturbanceLitterFluxes
  public :: PreDisturbanceIntegrateLitter
  public :: GenerateDamageAndLitterFluxes
  public :: SeedUpdate
  public :: UpdateRecruitL2FR
  public :: UpdateRecruitStoicH
  public :: SetRecruitL2FR
  
  logical, parameter :: debug  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  integer :: istat           ! return status code
  character(len=255) :: smsg ! Message string for deallocation errors
  
  integer, parameter :: dleafon_drycheck = 100 ! Drought deciduous leaves max days on check parameter

  real(r8), parameter :: decid_leaf_long_max = 1.0_r8 ! Maximum leaf lifespan for
                                                      !    deciduous PFTs [years]

  integer, parameter :: min_daysoff_dforcedflush = 30 ! This is the number of days that must had elapsed
                                                      ! since leaves had dropped, in order to forcably
                                                      ! flush leaves again.  This does not impact flushing
                                                      ! due to real moisture constraints, and will prevent
                                                      ! drought deciduous in perennially wet environments
                                                      ! that have been forced to drop their leaves, from
                                                      ! flushing them back immediately.

  integer, parameter  :: dd_offon_toler = 30          ! When flushing or shedding leaves, we check that
                                                      ! the dates are near last year's dates. This controls
                                                      ! the tolerance for deviating from last year.

  real(r8), parameter :: elongf_min = 0.05_r8         ! Minimum elongation factor. If elongation factor
                                                      !    reaches or falls below elongf_min, we assume
                                                      !    complete abscission.  This avoids carrying out
                                                      !    a residual amount of leaves, which may create
                                                      !    computational problems. The current threshold
                                                      !    is the same used in ED-2.2.

  real(r8), parameter :: smp_lwr_bound = -1000000._r8 ! Imposed soil matric potential lower bound for 
                                                      !    frozen or excessively dry soils, used when
                                                      !    computing water stress.
  ! ============================================================================

contains

  subroutine ZeroLitterFluxes( currentSite )

    ! This routine loops through all patches in a site
    ! and zero's the flux terms for the litter pools.
    ! This is typically called at the beginning of the dynamics
    ! call sequence.


    ! !ARGUMENTS
    type(ed_site_type), intent(inout), target  :: currentSite
    type(fates_patch_type), pointer               :: currentPatch

    integer :: el

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       do el=1,num_elements
          call currentPatch%litter(el)%ZeroFlux()
       end do
       currentPatch => currentPatch%older
    end do


    return
  end subroutine ZeroLitterFluxes

  ! =====================================================================================

  subroutine ZeroAllocationRates( currentSite )

    ! !ARGUMENTS
    type(ed_site_type), intent(inout), target  :: currentSite
    type(fates_patch_type), pointer               :: currentPatch
    type(fates_cohort_type), pointer              :: currentCohort

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))

       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))

          ! This sets turnover and growth rates to zero
          call currentCohort%prt%ZeroRates()

          currentCohort => currentCohort%shorter
       enddo
       currentPatch => currentPatch%older
    end do

    return
  end subroutine ZeroAllocationRates

  ! ============================================================================
  
  subroutine GenerateDamageAndLitterFluxes( csite, cpatch, bc_in )

    ! Arguments
    type(ed_site_type)  :: csite
    type(fates_patch_type) :: cpatch
    type(bc_in_type), intent(in) :: bc_in
    

    ! Locals
    type(fates_cohort_type), pointer :: ccohort    ! Current cohort
    type(fates_cohort_type), pointer :: ndcohort   ! New damage-class cohort
    type(litter_type), pointer :: litt     ! Points to the litter object
    type(elem_diag_type), pointer :: elflux_diags ! pointer to site level flux diagnostics object
    integer  :: cd               ! Damage class index
    integer  :: el               ! Element index
    integer  :: dcmpy            ! Decomposition pool index
    integer  :: c                ! CWD pool index
    real(r8) :: cd_frac          ! Fraction of trees damaged in this class transition
    real(r8) :: num_trees_cd     ! Number of trees to spawn into the new damage class cohort
    real(r8) :: crown_loss_frac  ! Fraction of crown lost from one damage class to next
    real(r8) :: branch_loss_frac ! Fraction of sap, structure and storage lost in branch
                                 ! fall during damage
    real(r8) :: leaf_loss        ! Mass lost to each organ during damage [kg]
    real(r8) :: repro_loss       ! "" [kg]
    real(r8) :: sapw_loss        ! "" [kg]
    real(r8) :: store_loss       ! "" [kg]
    real(r8) :: struct_loss      ! "" [kg]       
    real(r8) :: dcmpy_frac       ! fraction of mass going to each decomposition pool
    real(r8) :: SF_val_CWD_frac_adj(4) !SF_val_CWD_frac adjusted based on cohort dbh 
    
    if(hlm_use_tree_damage .ne. itrue) return

    if(.not.damage_time) return

    ccohort => cpatch%tallest
    do while (associated(ccohort))

       ! Ignore damage to new plants and non-woody plants
       if(prt_params%woody(ccohort%pft)==ifalse  ) cycle
       if(ccohort%isnew ) cycle

       associate( ipft     => ccohort%pft, & 
                  agb_frac => prt_params%allom_agb_frac(ccohort%pft), &
                  branch_frac => param_derived%branch_frac(ccohort%pft))
         
       do_dclass: do cd = ccohort%crowndamage+1, nlevdamage
          
          call GetDamageFrac(ccohort%crowndamage, cd, ipft, cd_frac)

          ! now to get the number of damaged trees we multiply by damage frac
          num_trees_cd = ccohort%n * cd_frac

          ! if non negligable lets create a new cohort and generate some litter
          if_numtrees: if (num_trees_cd > nearzero ) then

             ! Create a new damaged cohort
             allocate(ndcohort)  ! new cohort surviving but damaged
             if(hlm_use_planthydro.eq.itrue) call InitHydrCohort(csite,ndcohort)
             
             ! Initialize the PARTEH object and point to the
             ! correct boundary condition fields
             ndcohort%prt => null()

             call InitPRTObject(ndcohort%prt)
             call ndcohort%InitPRTBoundaryConditions()
             call ndcohort%ZeroValues()
             
             ! nc_canopy_d is the new cohort that gets damaged 
             call ccohort%Copy(ndcohort)
             
             ! new number densities - we just do damaged cohort here -
             ! undamaged at the end of the cohort loop once we know how many damaged to
             ! subtract
             
             ndcohort%n = num_trees_cd
             ndcohort%crowndamage = cd

             ! Remove these trees from the donor cohort
             ccohort%n = ccohort%n - num_trees_cd
             
             ! update crown area here - for cohort fusion and canopy organisation below 
             call carea_allom(ndcohort%dbh, ndcohort%n, csite%spread, &
                  ipft, ndcohort%crowndamage, ndcohort%c_area)
             
             call GetCrownReduction(cd-ccohort%crowndamage, crown_loss_frac)

             do_element: do el = 1, num_elements
                
                litt => cpatch%litter(el)
                elflux_diags => csite%flux_diags%elem(el)
                
                ! Reduce the mass of the newly damaged cohort
                ! Fine-roots are not damaged as of yet
                ! only above-ground sapwood,structure and storage in
                ! branches is damaged/removed
                branch_loss_frac = crown_loss_frac * branch_frac * agb_frac
                
                leaf_loss = ndcohort%prt%GetState(leaf_organ,element_list(el))*crown_loss_frac
                repro_loss = ndcohort%prt%GetState(repro_organ,element_list(el))*crown_loss_frac
                sapw_loss = ndcohort%prt%GetState(sapw_organ,element_list(el))*branch_loss_frac
                store_loss = ndcohort%prt%GetState(store_organ,element_list(el))*branch_loss_frac
                struct_loss = ndcohort%prt%GetState(struct_organ,element_list(el))*branch_loss_frac

                ! ------------------------------------------------------
                ! Transfer the biomass from the cohort's
                ! damage to the litter input fluxes
                ! ------------------------------------------------------
    
                do dcmpy=1,ndcmpy
                   dcmpy_frac = GetDecompyFrac(ipft,leaf_organ,dcmpy)
                   litt%leaf_fines_in(dcmpy) = litt%leaf_fines_in(dcmpy) + &
                        (store_loss+leaf_loss+repro_loss) * &
                        ndcohort%n * dcmpy_frac / cpatch%area
                end do

                elflux_diags%surf_fine_litter_input(ipft) = &
                     elflux_diags%surf_fine_litter_input(ipft) +  &
                     (store_loss+leaf_loss+repro_loss) * ndcohort%n
                
                call adjust_SF_CWD_frac(ndcohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)

                do c = 1,ncwd
                   litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + &
                        (sapw_loss + struct_loss) * &
                        SF_val_CWD_frac_adj(c) * ndcohort%n / &
                        cpatch%area
                   
                   elflux_diags%cwd_ag_input(c)  = elflux_diags%cwd_ag_input(c) + &
                        (struct_loss + sapw_loss) * &
                        SF_val_CWD_frac_adj(c) * ndcohort%n
                end do
                
             end do do_element

             ! Applying the damage to the cohort, does not need to happen
             ! in the element loop, it will loop inside that call
             call PRTDamageLosses(ndcohort%prt, leaf_organ, crown_loss_frac)
             call PRTDamageLosses(ndcohort%prt, repro_organ, crown_loss_frac)
             call PRTDamageLosses(ndcohort%prt, sapw_organ, branch_loss_frac)
             call PRTDamageLosses(ndcohort%prt, store_organ, branch_loss_frac)
             call PRTDamageLosses(ndcohort%prt, struct_organ, branch_loss_frac)
                                
             
             !----------- Insert new cohort into the linked list
             ! This list is going tall to short, lets add this new
             ! cohort into a taller position so we don't hit it again
             ! as the loop traverses
             ! --------------------------------------------------------------!
             
             ndcohort%shorter => ccohort
             if(associated(ccohort%taller))then
                ndcohort%taller => ccohort%taller
                ccohort%taller%shorter => ndcohort
             else
                cpatch%tallest => ndcohort
                ndcohort%taller => null()
             endif
             ccohort%taller => ndcohort
             
          end if if_numtrees

       end do do_dclass

       end associate
       ccohort => ccohort%shorter
    enddo
    
    return
  end subroutine GenerateDamageAndLitterFluxes

  ! ============================================================================

  subroutine PreDisturbanceLitterFluxes( currentSite, currentPatch, bc_in )

    ! -----------------------------------------------------------------------------------
    !
    ! This subroutine calculates all of the different litter input and output fluxes
    ! associated with seed turnover, seed influx, litterfall from live and
    ! dead plants, germination, and fragmentation.
    !
    ! At this time we do not have explicit herbivory, and burning losses to litter
    ! are handled elsewhere.
    !
    ! Note: The processes conducted here DO NOT handle litter fluxes associated
    !       with disturbance.  Those fluxes are handled elsewhere (EDPatchDynamcisMod)
    !       because the fluxes are potentially cross patch, and also dealing
    !       patch areas that are changing.
    !
    ! -----------------------------------------------------------------------------------


    ! !ARGUMENTS
    type(ed_site_type), intent(inout)  :: currentSite
    type(fates_patch_type), intent(inout) :: currentPatch
    type(bc_in_type), intent(in)       :: bc_in

    !
    ! !LOCAL VARIABLES:

    integer :: el                          ! Litter element loop index
    integer :: nlev_eff_decomp             ! Number of active layers over which
    ! fragmentation fluxes are transfered
    !------------------------------------------------------------------------------------

    ! Calculate the fragmentation rates
    call fragmentation_scaler(currentPatch, bc_in)

    do el = 1, num_elements

       associate( litt => currentPatch%litter(el), &
                  site_mass => currentSite%mass_balance(el), &
                  diag => currentSite%flux_diags%elem(el))

         ! Calculate loss rate of viable seeds to litter
         call SeedDecay(litt, currentPatch, bc_in)
         
         ! Calculate seed germination rate, the status flags prevent
         ! germination from occuring when the site is in a drought
         ! (for drought deciduous) or too cold (for cold deciduous)
         call SeedGermination(litt, currentSite%cstatus, currentSite%dstatus(1:numpft), bc_in, currentPatch)
         
         ! Send fluxes from newly created litter into the litter pools
         ! This litter flux is from non-disturbance inducing mortality, as well
         ! as litter fluxes from live trees
         call CWDInput(currentSite, currentPatch, litt,bc_in)
         
         ! Only calculate fragmentation flux over layers that are active
         ! (RGK-Mar2019) SHOULD WE MAX THIS AT 1? DONT HAVE TO
         
         nlev_eff_decomp = max(bc_in%max_rooting_depth_index_col, 1)
         call CWDOut(litt,currentPatch%fragmentation_scaler,nlev_eff_decomp)
         
         ! Fragmentation flux to soil decomposition model [kg/site/day]
         site_mass%frag_out = site_mass%frag_out + currentPatch%area * &
              ( sum(litt%ag_cwd_frag) + sum(litt%bg_cwd_frag) + &
              sum(litt%leaf_fines_frag) + sum(litt%root_fines_frag) + &
              sum(litt%seed_decay) + sum(litt%seed_germ_decay))
         
         ! Track total seed decay diagnostic in [kg/m2/day]
         diag%tot_seed_turnover = diag%tot_seed_turnover + &
              (sum(litt%seed_decay) + sum(litt%seed_germ_decay))*currentPatch%area*area_inv

       end associate
    end do


    return
  end subroutine PreDisturbanceLitterFluxes

  ! =====================================================================================

  subroutine PreDisturbanceIntegrateLitter(currentPatch)

    ! -----------------------------------------------------------------------------------
    ! This step applies the litter fluxes to the prognostic state variables.
    ! This procedure is called in response to fluxes generated from:
    ! 1) seed rain,
    ! 2) non-disturbance generating turnover
    ! 3) litter fall from living plants
    ! 4) fragmentation
    !
    ! This routine does NOT accomodate the litter fluxes associated with
    ! disturbance generation.  That will happen after this call.
    ! Fluxes associated with FIRE also happen after this step.
    !
    ! All states are in units kg/m2
    ! All fluxes are in units kg/m2/day
    ! The integration step is 1 day, thus time is implied
    !
    ! -----------------------------------------------------------------------------------

    ! Arguments
    type(fates_patch_type),intent(inout),target :: currentPatch


    ! Locals
    type(litter_type), pointer :: litt
    integer :: el          ! Loop counter for litter element type
    integer :: pft         ! pft loop counter
    integer :: c           ! CWD loop counter
    integer :: nlevsoil    ! number of soil layers
    integer :: ilyr        ! soil layer loop counter
    integer :: dcmpy       ! decomposability index

    do el = 1, num_elements

       litt => currentPatch%litter(el)

       ! Update the bank of viable seeds
       ! -----------------------------------------------------------------------------------

       do pft = 1,numpft
          litt%seed(pft) = litt%seed(pft) + &
               litt%seed_in_local(pft) +   &
               litt%seed_in_extern(pft) -  &
               litt%seed_decay(pft) -      &
               litt%seed_germ_in(pft)

          ! Note that the recruitment scheme will use seed_germ
          ! for its construction costs.
          litt%seed_germ(pft) = litt%seed_germ(pft) + &
               litt%seed_germ_in(pft) - &
               litt%seed_germ_decay(pft)

       enddo

       ! Update the Coarse Woody Debris pools (above and below)
       ! -----------------------------------------------------------------------------------
       nlevsoil = size(litt%bg_cwd,dim=2)
       do c = 1,ncwd
          litt%ag_cwd(c) = litt%ag_cwd(c)  + litt%ag_cwd_in(c) - litt%ag_cwd_frag(c)
          do ilyr=1,nlevsoil
             litt%bg_cwd(c,ilyr) = litt%bg_cwd(c,ilyr) &
                  + litt%bg_cwd_in(c,ilyr) &
                  - litt%bg_cwd_frag(c,ilyr)
          enddo
       end do

       ! Update the fine litter pools from leaves and fine-roots
       ! -----------------------------------------------------------------------------------

       do dcmpy = 1,ndcmpy

          litt%leaf_fines(dcmpy) = litt%leaf_fines(dcmpy) &
               + litt%leaf_fines_in(dcmpy)              &
               - litt%leaf_fines_frag(dcmpy)
          do ilyr=1,nlevsoil
             litt%root_fines(dcmpy,ilyr) = litt%root_fines(dcmpy,ilyr) &
                  + litt%root_fines_in(dcmpy,ilyr)      &
                  - litt%root_fines_frag(dcmpy,ilyr)
          enddo

       end do

    end do     ! litter element loop

    return
  end subroutine PreDisturbanceIntegrateLitter



  ! ============================================================================

  subroutine trim_canopy( currentSite )
    !
    ! !DESCRIPTION:
    ! Canopy trimming / leaf optimisation. Removes leaves in negative annual carbon balance.
    !
    ! !USES:

    ! !ARGUMENTS
    type (ed_site_type),intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (fates_cohort_type) , pointer :: currentCohort
    type (fates_patch_type)  , pointer :: currentPatch

    integer  :: z                     ! leaf layer
    integer  :: ipft                  ! pft index
    logical  :: trimmed               ! was this layer trimmed in this year? If not expand the canopy.
    real(r8) :: tar_bl                ! target leaf biomass       (leaves flushed, trimmed)
    real(r8) :: tar_bfr               ! target fine-root biomass  (leaves flushed, trimmed)
    real(r8) :: bfr_per_bleaf         ! ratio of fine root per leaf biomass
    real(r8) :: sla_levleaf           ! sla at leaf level z
    real(r8) :: nscaler_levleaf       ! nscaler value at leaf level z
    integer  :: cl                    ! canopy layer index
    real(r8) :: kn                    ! nitrogen decay coefficient
    real(r8) :: sla_max               ! Observational constraint on how large sla (m2/gC) can become
    real(r8) :: leaf_c                ! leaf carbon [kg]
    real(r8) :: sapw_c                ! sapwood carbon [kg]
    real(r8) :: store_c               ! storage carbon [kg]
    real(r8) :: struct_c              ! structure carbon [kg]
    real(r8) :: leaf_inc              ! LAI-only portion of the vegetation increment of dinc_vai
    real(r8) :: lai_canopy_above      ! the LAI in the canopy layers above the layer of interest
    real(r8) :: lai_layers_above      ! the LAI in the leaf layers, within the current canopy,
    ! above the leaf layer of interest
    real(r8) :: lai_current           ! the LAI in the current leaf layer
    real(r8) :: cumulative_lai        ! whole canopy cumulative LAI, top down, to the leaf layer of interest
    real(r8) :: cumulative_lai_cohort ! cumulative LAI within the current cohort only

    ! Temporary diagnostic ouptut

    ! LAPACK linear least squares fit variables
    ! The standard equation for a linear fit, y = mx + b, is converted to a linear system, AX=B and has
    ! the form: [n  sum(x); sum(x)  sum(x^2)] * [b; m]  = [sum(y); sum(x*y)] where
    ! n is the number of leaf layers
    ! x is yearly_net_uptake minus the leaf cost aka the net-net uptake
    ! y is the cumulative lai for the current cohort
    ! b is the y-intercept i.e. the cumulative lai that has zero net-net uptake
    ! m is the slope of the linear fit
    integer :: nll = 3                    ! Number of leaf layers to fit a regression to for calculating the optimum lai
    character(1) :: trans = 'N'           ! Input matrix is not transposed

    integer, parameter :: m = 2, n = 2    ! Number of rows and columns, respectively, in matrix A
    integer, parameter :: nrhs = 1        ! Number of columns in matrix B and X
    integer, parameter :: workmax = 100   ! Maximum iterations to minimize work

    integer :: lda = m, ldb = n           ! Leading dimension of A and B, respectively
    integer :: lwork                      ! Dimension of work array
    integer :: info                       ! Procedure diagnostic ouput

    real(r8) :: nnu_clai_a(m,n)           ! LHS of linear least squares fit, A matrix
    real(r8) :: nnu_clai_b(m,nrhs)        ! RHS of linear least squares fit, B matrix
    real(r8) :: work(workmax)             ! work array

    real(r8) :: initial_trim              ! Initial trim
    real(r8) :: optimum_trim              ! Optimum trim value

    real(r8) :: target_c_area

    real(r8) :: pft_leaf_lifespan         ! Leaf lifespan of each PFT [years]
    real(r8) :: leaf_long                 ! temporary leaf lifespan before accounting for deciduousness 
    !----------------------------------------------------------------------

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))

       ! Add debug diagnstic output to determine which patch
       if (debug) then
          write(fates_log(),*) 'Current patch:', currentPatch%patchno
          write(fates_log(),*) 'Current patch cohorts:', currentPatch%num_cohorts
       endif

       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))

          ! Save off the incoming trim
          initial_trim = currentCohort%canopy_trim


          ! Add debug diagnostic output to determine which cohort
          if (debug) then
             write(fates_log(),*) 'Starting canopy trim:', initial_trim
          endif

          trimmed = .false.
          ipft = currentCohort%pft
          call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread,currentCohort%pft,&
               currentCohort%crowndamage, currentCohort%c_area)

          leaf_c   = currentCohort%prt%GetState(leaf_organ, carbon12_element)

          call  tree_lai_sai(leaf_c, currentCohort%pft, currentCohort%c_area, currentCohort%n,           &
               currentCohort%canopy_layer, currentPatch%canopy_layer_tlai, currentCohort%vcmax25top,   &
               currentCohort%dbh, currentCohort%crowndamage, currentCohort%canopy_trim, &
               currentCohort%efstem_coh, 0, currentCohort%treelai, currentCohort%treesai )

          currentCohort%nv      = count((currentCohort%treelai+currentCohort%treesai) .gt. dlower_vai(:)) + 1

          if (currentCohort%nv > nlevleaf)then
             write(fates_log(),*) 'nv > nlevleaf',currentCohort%nv, &
                  currentCohort%treelai,currentCohort%treesai, &
                  currentCohort%c_area,currentCohort%n,leaf_c
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

          ! Find target leaf biomass. Here we assume that leaves would be fully flushed 
          ! (elongation factor = 1)
          call bleaf(currentcohort%dbh,ipft,&
               currentCohort%crowndamage, currentcohort%canopy_trim,1.0_r8, tar_bl)

          if ( int(prt_params%allom_fmode(ipft)) .eq. 1 ) then
             ! only query fine root biomass if using a fine root allometric model that takes leaf trim into account
             call bfineroot(currentcohort%dbh,ipft,currentcohort%canopy_trim, &
                  currentcohort%l2fr,1.0_r8, tar_bfr)
             bfr_per_bleaf = tar_bfr/tar_bl
          endif

          ! Identify current canopy layer (cl)
          cl = currentCohort%canopy_layer

          ! Get leaf lifespan- depends on canopy layer
          if  (cl .eq. 1 ) then
             leaf_long = sum(prt_params%leaf_long(ipft,:))
          else
             leaf_long = sum(prt_params%leaf_long_ustory(ipft,:))
          end if
          

          ! PFT-level maximum SLA value, even if under a thick canopy (same units as slatop)
          sla_max = prt_params%slamax(ipft)

          ! Initialize nnu_clai_a
          nnu_clai_a(:,:) = 0._r8
          nnu_clai_b(:,:) = 0._r8

          !Leaf cost vs net uptake for each leaf layer.
          do z = 1, currentCohort%nv

             ! Calculate the cumulative total vegetation area index (no snow occlusion, stems and leaves)
             leaf_inc    = dinc_vai(z) * &
                  currentCohort%treelai/(currentCohort%treelai+currentCohort%treesai)
             
             ! Now calculate the cumulative top-down lai of the current layer's midpoint within the current cohort
             lai_layers_above      = (dlower_vai(z) - dinc_vai(z)) * &
                  currentCohort%treelai/(currentCohort%treelai+currentCohort%treesai)
             lai_current           = min(leaf_inc, currentCohort%treelai - lai_layers_above)
             cumulative_lai_cohort = lai_layers_above + 0.5*lai_current

             ! Now add in the lai above the current cohort for calculating the sla leaf level
             lai_canopy_above  = sum(currentPatch%canopy_layer_tlai(1:cl-1))
             cumulative_lai    = lai_canopy_above + cumulative_lai_cohort

             ! There was activity this year in this leaf layer.  This should only occur for bottom most leaf layer
             if (currentCohort%year_net_uptake(z) /= 999._r8)then

                ! Calculate sla_levleaf following the sla profile with overlying leaf area
                ! Scale for leaf nitrogen profile
                kn = decay_coeff_vcmax(currentCohort%vcmax25top, &
                     prt_params%leafn_vert_scaler_coeff1(ipft), &
                     prt_params%leafn_vert_scaler_coeff2(ipft))
                
                ! Nscaler value at leaf level z
                nscaler_levleaf = exp(-kn * cumulative_lai)
                ! Sla value at leaf level z after nitrogen profile scaling (m2/gC)
                sla_levleaf = min(sla_max,prt_params%slatop(ipft)/nscaler_levleaf)

                ! Find the realised leaf lifespan, depending on the leaf phenology.
                if (prt_params%season_decid(ipft) ==  itrue) then
                   ! Cold-deciduous costs. Assume time-span to be 1 year to be consistent
                   ! with FATES default
                   pft_leaf_lifespan = decid_leaf_long_max

                elseif (any(prt_params%stress_decid(ipft) == [ihard_stress_decid,isemi_stress_decid]) )then
                   ! Drought-decidous costs. Assume time-span to be the least between
                   !    1 year and the life span provided by the parameter file.
                   pft_leaf_lifespan = &
                      min(decid_leaf_long_max,leaf_long)

                else !evergreen costs
                   pft_leaf_lifespan = leaf_long
                end if

                ! Leaf cost at leaf level z (kgC m-2 year-1) accounting for sla profile
                ! (Convert from SLA in m2g-1 to m2kg-1)
                currentCohort%leaf_cost = &
                   1.0_r8/(sla_levleaf*pft_leaf_lifespan*g_per_kg)


                if ( int(prt_params%allom_fmode(ipft)) == 1 ) then
                   ! if using trimmed leaf for fine root biomass allometry, add the cost of the root increment
                   ! to the leaf increment; otherwise do not.
                   currentCohort%leaf_cost = currentCohort%leaf_cost + &
                        1.0_r8/(sla_levleaf*g_per_kg) * &
                        bfr_per_bleaf / prt_params%root_long(ipft)
                end if
                currentCohort%leaf_cost = currentCohort%leaf_cost * &
                     (prt_params%grperc(ipft) + 1._r8)

                ! Construct the arrays for a least square fit of the net_net_uptake versus the cumulative lai
                ! if at least nll leaf layers are present in the current cohort and only for the bottom nll
                ! leaf layers.
                if (currentCohort%nv > nll .and. currentCohort%nv - z < nll) then

                   ! Build the A matrix for the LHS of the linear system. A = [n  sum(x); sum(x)  sum(x^2)]
                   ! where n = nll and x = yearly_net_uptake-leafcost
                   nnu_clai_a(1,1) = nnu_clai_a(1,1) + 1 ! Increment for each layer used
                   nnu_clai_a(1,2) = nnu_clai_a(1,2) + currentCohort%year_net_uptake(z) - currentCohort%leaf_cost
                   nnu_clai_a(2,1) = nnu_clai_a(1,2)
                   nnu_clai_a(2,2) = nnu_clai_a(2,2) + (currentCohort%year_net_uptake(z) - currentCohort%leaf_cost)**2

                   ! Build the B matrix for the RHS of the linear system. B = [sum(y); sum(x*y)]
                   ! where x = yearly_net_uptake-leafcost and y = cumulative_lai_cohort
                   nnu_clai_b(1,1) = nnu_clai_b(1,1) + cumulative_lai_cohort
                   nnu_clai_b(2,1) = nnu_clai_b(2,1) + (cumulative_lai_cohort * &
                        (currentCohort%year_net_uptake(z) - currentCohort%leaf_cost))
                end if

                ! Check leaf cost against the yearly net uptake for that cohort leaf layer
                if (currentCohort%year_net_uptake(z) < currentCohort%leaf_cost) then
                   ! Make sure the cohort trim fraction is great than the pft trim limit
                   if (currentCohort%canopy_trim > EDPftvarcon_inst%trim_limit(ipft)) then

                      ! keep trimming until none of the canopy is in negative carbon balance.
                      if (currentCohort%height > EDPftvarcon_inst%hgt_min(ipft)) then
                         currentCohort%canopy_trim = currentCohort%canopy_trim - &
                              EDPftvarcon_inst%trim_inc(ipft)

                         trimmed = .true.

                      endif ! height check
                   endif ! trim limit check
                endif ! net uptake check
             endif ! leaf activity check
          enddo ! z, leaf layer loop

          ! Compute the optimal cumulative lai based on the cohort net-net uptake profile if at least 2 leaf layers
          if (nnu_clai_a(1,1) > 1) then

             ! Compute the optimum size of the work array
             lwork = -1 ! Ask sgels to compute optimal number of entries for work
             call dgels(trans, m, n, nrhs, nnu_clai_a, lda, nnu_clai_b, ldb, work, lwork, info)
             lwork = int(work(1)) ! Pick the optimum.  TBD, can work(1) come back with greater than work size?

             ! Compute the minimum of 2-norm of of the least squares fit to solve for X
             ! Note that dgels returns the solution by overwriting the nnu_clai_b array.
             ! The result has the form: X = [b; m]
             ! where b = y-intercept (i.e. the cohort lai that has zero yearly net-net uptake)
             ! and m is the slope of the linear fit
             call dgels(trans, m, n, nrhs, nnu_clai_a, lda, nnu_clai_b, ldb, work, lwork, info)

             if (info < 0) then
                write(fates_log(),*) 'LLSF optimium LAI calculation returned illegal value'
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif

             if (debug) then
                write(fates_log(),*) 'LLSF optimium LAI (intercept,slope):', nnu_clai_b
                write(fates_log(),*) 'LLSF optimium LAI:', nnu_clai_b(1,1)
                write(fates_log(),*) 'LLSF optimium LAI info:', info
                write(fates_log(),*) 'LAI fraction (optimum_lai/cumulative_lai):', nnu_clai_b(1,1) / cumulative_lai_cohort
             endif

             ! Calculate the optimum trim based on the initial canopy trim value
             if (cumulative_lai_cohort > 0._r8) then  ! Sometime cumulative_lai comes in at 0.0?

                !
                optimum_trim = (nnu_clai_b(1,1) / cumulative_lai_cohort) * initial_trim

                ! Determine if the optimum trim value makes sense.  The smallest cohorts tend to have unrealistic fits.
                if (optimum_trim > 0. .and. optimum_trim < 1.) then
                   currentCohort%canopy_trim = optimum_trim

                   trimmed = .true.

                endif
             endif
          endif

          ! Reset activity for the cohort for the start of the next year
          currentCohort%year_net_uptake(:) = 999.0_r8

          ! Add to trim fraction if cohort not trimmed at all
          if ( (.not.trimmed) .and.currentCohort%canopy_trim < 1.0_r8)then
             currentCohort%canopy_trim = currentCohort%canopy_trim + EDPftvarcon_inst%trim_inc(ipft)
          endif

          if ( debug ) then
             write(fates_log(),*) 'trimming:',currentCohort%canopy_trim
          endif

          ! currentCohort%canopy_trim = 1.0_r8 !FIX(RF,032414) this turns off ctrim for now.
          currentCohort => currentCohort%shorter
       enddo
       currentPatch => currentPatch%older

    enddo

  end subroutine trim_canopy

  ! ============================================================================
  subroutine phenology( currentSite, bc_in )
    !
    ! !DESCRIPTION:
    ! Phenology.
    !
    ! !USES:
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use EDParamsMod, only : ED_val_phen_a, ED_val_phen_b, ED_val_phen_c
    use EDParamsMod, only : ED_val_phen_chiltemp
    use EDParamsMod, only : ED_val_phen_mindayson
    use EDParamsMod, only : ED_val_phen_ncolddayslim
    use EDParamsMod, only : ED_val_phen_coldtemp
    use EDBtranMod, only  : check_layer_water
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    type(bc_in_type),   intent(in)            :: bc_in

    !
    ! !LOCAL VARIABLES:

    type(fates_patch_type),pointer :: cpatch
    integer  :: model_day_int     ! integer model day 1 - inf
    integer  :: ncolddays         ! no days underneath the threshold for leaf drop
    integer  :: i_wmem            ! Loop counter for water mem days
    integer  :: i_tmem            ! Loop counter for veg temp mem days
    integer  :: ipft              ! plant functional type index
    integer  :: j                 ! Soil layer index
    real(r8) :: mean_10day_liqvol ! mean soil liquid volume over last 10 days [m3/m3]
    real(r8) :: mean_10day_smp    ! mean soil matric potential over last 10 days [mm]
    real(r8) :: leaf_c            ! leaf carbon [kg]
    real(r8) :: fnrt_c            ! fineroot carbon [kg]
    real(r8) :: sapw_c            ! sapwood carbon [kg]
    real(r8) :: store_c           ! storage carbon [kg]
    real(r8) :: struct_c          ! structure carbon [kg]
    real(r8) :: gdd_threshold     ! GDD accumulation function,
    real(r8) :: rootfrac_notop    ! Total rooting fraction excluding the top soil layer
    integer  :: ncdstart          ! beginning of counting period for chilling degree days.
    integer  :: gddstart          ! beginning of counting period for growing degree days.
    integer  :: nlevroot          ! Number of rooting levels to consider
    real(r8) :: temp_in_C         ! daily averaged temperature in celsius
    real(r8) :: temp_wgt          ! canopy area weighting factor for daily average
                                  ! vegetation temperature calculation
    real(r8) :: elongf_prev       ! Elongation factor from previous time
    real(r8) :: elongf_1st        ! First guess for elongation factor
    integer  :: ndays_pft_leaf_lifespan ! PFT life span of drought deciduous [days].
                                        !    This is the shortest between the PFT leaf 
                                        !    lifespan and the maximum lifespan of drought 
                                        !    deciduous (see parameter decid_leaf_long_max
                                        !    at the beginning of this file).
     real(r8) :: phen_drought_threshold ! For drought hard-deciduous, this is the threshold
                                        !   below which plants will abscise leaves, and
                                        !   above which plants will flush leaves. For semi-
                                        !   deciduous plants, this is the threshold below
                                        !   which abscission will be complete. This depends
                                        !   on the sign. If positive, these are soil
                                        !   volumetric water content [m3/m3]. If negative,
                                        !   the values are soil matric potential [mm]. Not
                                        !   used for non-deciduous plants. Ignored for 
                                        !   non-deciduous plants.
     real(r8) :: phen_moist_threshold   ! For semi-deciduous, this is the threshold above 
                                        !    which flushing will be complete.  This depends
                                        !    on the sign. If positive, these are soil
                                        !    volumetric water content [m3/m3]. If negative,
                                        !    the values are soil matric potential [mm].
                                        !    Ignored for hard-deciduous and evergreen 
                                        !    plants.
     real(r8) :: phen_doff_time         ! Minimum number of days that plants must remain
                                        !   leafless before flushing leaves again.

    ! Logical tests to make code more readable
    logical  :: smoist_below_threshold   ! Is soil moisture below threshold?
    logical  :: recent_flush             ! Last full flushing event is still very recent.
    logical  :: recent_abscission        ! Last abscission event is still very recent.
    logical  :: exceed_min_on_period     ! Have leaves been flushed for a minimum period of time?
    logical  :: exceed_min_off_period    ! Have leaves been off for a minimum period of time?
    logical  :: prolonged_on_period      ! Has leaves been flushed for too long?
    logical  :: prolonged_off_period     ! Have leaves been abscissed for too long?
    logical  :: last_flush_long_ago      ! Has it been a very long time since last flushing?


    ! This is the integer model day. The first day of the simulation is 1, and it
    ! continues monotonically, indefinitely
    ! Advance it. (this should be a global, no reason
    ! for site level, but we don't have global scalars in the
    ! restart file)
    currentSite%phen_model_date = currentSite%phen_model_date + 1
    model_day_int = currentSite%phen_model_date


    ! Parameter of drought decid leaf loss in mm in top layer...FIX(RF,032414)
    ! - this is arbitrary and poorly understood. Needs work. ED_
    !Parameters: defaults from Botta et al. 2000 GCB,6 709-725
    !Parameters, default from from SDGVM model of senesence

    temp_in_C = 0._r8
    temp_wgt = 0._r8
    cpatch => CurrentSite%oldest_patch
    do while(associated(cpatch))
       temp_in_C = temp_in_C + cpatch%tveg24%GetMean()*cpatch%total_canopy_area
       temp_wgt = temp_wgt + cpatch%total_canopy_area
       cpatch => cpatch%younger
    end do
    if(temp_wgt>nearzero)then
       temp_in_C = temp_in_C/temp_wgt - tfrz
    else
       ! If there is no canopy area, we use the veg temperature
       ! of the first patch, which is the forcing air temperature
       ! as defined in CLM/ELM. The forcing air temperature
       ! should be the same among all patches. (Although
       ! it is unlikely there are more than 1 in this scenario)
       temp_in_C = CurrentSite%oldest_patch%tveg24%GetMean() - tfrz
    end if

    !-----------------Cold Phenology--------------------!

    !Zero growing degree and chilling day counters
    if (currentSite%lat > 0)then
       ncdstart = 270  !Northern Hemisphere begining November
       gddstart = 1    !Northern Hemisphere begining January
    else
       ncdstart = 120  !Southern Hemisphere beginning May
       gddstart = 181  !Northern Hemisphere begining July
    endif

    ! Count the number of chilling days over a seasonal window.
    ! For comparing against GDD, we start calculating chilling
    ! in the late autumn.
    ! This value is used to determine the GDD exceedance threshold
    if (hlm_day_of_year == ncdstart)then
       currentSite%nchilldays = 0
    endif

    !Accumulate growing/chilling days after start of counting period
    if (temp_in_C  <  ED_val_phen_chiltemp)then
       currentSite%nchilldays = currentSite%nchilldays + 1
    endif

    !GDD accumulation function, which also depends on chilling days.
    !  -68 + 638 * (-0.001 * ncd)
    gdd_threshold = ED_val_phen_a + ED_val_phen_b*exp(ED_val_phen_c*real(currentSite%nchilldays,r8))

    !Accumulate temperature of last 10 days.
    currentSite%vegtemp_memory(2:num_vegtemp_mem) = currentSite%vegtemp_memory(1:num_vegtemp_mem-1)
    currentSite%vegtemp_memory(1) = temp_in_C

    !count number of days for leaves off
    ncolddays = 0
    do i_tmem = 1,num_vegtemp_mem
       if (currentSite%vegtemp_memory(i_tmem) < ED_val_phen_coldtemp)then
          ncolddays = ncolddays + 1
       endif
    enddo

    ! Here is where we do the GDD accumulation calculation
    !
    ! reset GDD on set dates
    if (hlm_day_of_year == gddstart)then
       currentSite%grow_deg_days = 0._r8
    endif
    !
    ! accumulate the GDD using daily mean temperatures
    ! Don't accumulate GDD during the growing season (that wouldn't make sense)
    if (temp_in_C .gt. 0._r8 .and. currentSite%cstatus == phen_cstat_iscold) then
       currentSite%grow_deg_days = currentSite%grow_deg_days + temp_in_C
    endif

    !this logic is to prevent GDD accumulating after the leaves have fallen and before the
    ! beginnning of the accumulation period, to prevend erroneous autumn leaf flushing.
    if(model_day_int> ndays_per_year)then !only do this after the first year to prevent odd behaviour

       if(currentSite%lat .gt. 0.0_r8)then !Northern Hemisphere
          ! In the north, don't accumulate when we are past the leaf fall date.
          ! Accumulation starts on day 1 of year in NH.
          ! The 180 is to prevent going into an 'always off' state after initialization
          if( model_day_int .gt. currentSite%cleafoffdate.and.hlm_day_of_year.gt.180)then !
             currentSite%grow_deg_days = 0._r8
          endif
       else !Southern Hemisphere
          ! In the South, don't accumulate after the leaf off date, and before the start of
          ! the accumulation phase (day 181).
          if(model_day_int .gt. currentSite%cleafoffdate.and.hlm_day_of_year.lt.gddstart) then!
             currentSite%grow_deg_days = 0._r8
          endif
       endif
    endif !year1

    ! Calculate the number of days since the leaves last came on
    ! and off. If this is the beginning of the simulation, that day might
    ! not had occured yet, so set it to last year to get things rolling

    if (model_day_int < currentSite%cleafoffdate) then
       currentSite%cndaysleafoff = model_day_int - (currentSite%cleafoffdate - ndays_per_year)
    else
       currentSite%cndaysleafoff = model_day_int - currentSite%cleafoffdate
    end if

    if (model_day_int < currentSite%cleafondate) then
       currentSite%cndaysleafon = model_day_int - (currentSite%cleafondate - ndays_per_year)
    else
       currentSite%cndaysleafon = model_day_int - currentSite%cleafondate
    end if



    !LEAF ON: COLD DECIDUOUS. Needs to
    !1) have exceeded the growing degree day threshold
    !2) The leaves should not be on already
    !3) There should have been at least one chilling day in the counting period.
    !   this prevents tropical or warm climate plants that are "cold-deciduous"
    !   from ever re-flushing after they have reached their maximum age (thus
    !   preventing them from competing

    if ( any(currentSite%cstatus == [phen_cstat_iscold,phen_cstat_nevercold]) .and. &
         (currentSite%grow_deg_days > gdd_threshold) .and. &
         (currentSite%cndaysleafoff > ED_val_phen_mindayson) .and. &
         (currentSite%nchilldays >= 1)) then
       currentSite%cstatus = phen_cstat_notcold  ! Set to not-cold status (leaves can come on)
       currentSite%cleafondate = model_day_int
       currentSite%cndaysleafon = 0
       currentSite%grow_deg_days = 0._r8 ! zero GDD for the rest of the year until counting season begins.
       if ( debug ) write(fates_log(),*) 'leaves on'
    endif !GDD




    !LEAF OFF: COLD THRESHOLD
    !Needs to:
    !1) have exceeded the number of cold days threshold
    !2) have exceeded the minimum leafon time.
    !3) The leaves should not be off already
    !4) The day of simulation should be larger than the counting period.


    if ( (currentSite%cstatus == phen_cstat_notcold) .and. &
         (model_day_int > num_vegtemp_mem)      .and. &
         (ncolddays > ED_val_phen_ncolddayslim) .and. &
         (currentSite%cndaysleafon > ED_val_phen_mindayson) )then

       currentSite%grow_deg_days  = 0._r8          ! The equations for Botta et al
       ! are for calculations of
       ! first flush, but if we dont
       ! clear this value, it will cause
       ! leaves to flush later in the year
       currentSite%cstatus       = phen_cstat_iscold  ! alter status of site to 'leaves off'
       currentSite%cleafoffdate = model_day_int       ! record leaf off date
       currentSite%cndaysleafoff = 0

       if ( debug ) write(fates_log(),*) 'leaves off'
    endif

    ! LEAF OFF: COLD LIFESPAN THRESHOLD
    ! NOTE: Some areas of the planet will never generate a cold day
    ! and thus %nchilldays will never go from zero to 1.  The following logic
    ! when coupled with this fact will essentially prevent cold-deciduous
    ! plants from re-emerging in areas without at least some cold days
    
    if( (currentSite%cstatus == phen_cstat_notcold)  .and. &
        (currentSite%cndaysleafoff > 400)) then   ! remove leaves after a whole year,
                                                  ! when there is no 'off' period.
       currentSite%grow_deg_days  = 0._r8

       currentSite%cstatus = phen_cstat_nevercold  ! alter status of site to imply that this
       ! site is never really cold enough
       ! for cold deciduous
       currentSite%cleafoffdate = model_day_int    ! record leaf off date
       currentSite%cndaysleafoff = 0

       if ( debug ) write(fates_log(),*) 'leaves off'
    endif



    ! Loop through every PFT to assign the elongation factor. 
    ! Add PFT look to account for different PFT rooting depth profiles.
    pft_elong_loop: do ipft=1,numpft

       ! Copy values to a local variable to make code more legible.
       phen_drought_threshold = prt_params%phen_drought_threshold(ipft)
       phen_moist_threshold   = prt_params%phen_moist_threshold  (ipft)
       phen_doff_time         = prt_params%phen_doff_time        (ipft)


       ! Update soil moisture information memory (we always track the last 10 days)
       do i_wmem = numWaterMem,2,-1 !shift memory to previous day, to make room for current day
          currentSite%liqvol_memory(i_wmem,ipft) = currentSite%liqvol_memory(i_wmem-1,ipft)
          currentSite%smp_memory   (i_wmem,ipft) = currentSite%smp_memory   (i_wmem-1,ipft)
       end do

       ! Find the rooting depth distribution for PFT
       call set_root_fraction( currentSite%rootfrac_scr, ipft, currentSite%zi_soil, &
                               bc_in%max_rooting_depth_index_col )
       nlevroot = max(2,min(ubound(currentSite%zi_soil,1),bc_in%max_rooting_depth_index_col))

       ! The top most layer is typically very thin (~ 2cm) and dries rather quickly. Despite
       ! being thin, it can have a non-negligible rooting fraction (e.g., using 
       ! exponential_2p_root_profile with default parameters make the top layer to contain
       ! about 7% of the total fine root density).  To avoid overestimating dryness, we 
       ! ignore the top layer when calculating the memory.
       rootfrac_notop = sum(currentSite%rootfrac_scr(2:nlevroot))
       if ( rootfrac_notop <= nearzero ) then
          ! Unlikely, but just in case all roots are in the first layer, we use the second
          ! layer the second layer (to avoid FPE issues).
          currentSite%rootfrac_scr(2) = 1.0_r8
          rootfrac_notop              = 1.0_r8
       end if

       ! Set the memory to be the weighted average of the soil properties, using the
       ! root fraction of each layer (except the topmost one) as the weighting factor.

       currentSite%liqvol_memory(1,ipft) = sum( bc_in%h2o_liqvol_sl     (2:nlevroot) * &
                                                currentSite%rootfrac_scr(2:nlevroot) ) / &
                                                rootfrac_notop
       currentSite%smp_memory   (1,ipft)  = 0._r8
       do j = 2,nlevroot
          if(check_layer_water(bc_in%h2o_liqvol_sl(j),bc_in%tempk_sl(j)) ) then
             currentSite%smp_memory   (1,ipft) = currentSite%smp_memory   (1,ipft) + & 
                  bc_in%smp_sl            (j) * &
                  currentSite%rootfrac_scr(j)  / &
                  rootfrac_notop
          else
             ! Nominal extreme suction for frozen or unreasonably dry soil
             currentSite%smp_memory   (1,ipft) = currentSite%smp_memory   (1,ipft) + & 
                  smp_lwr_bound * &
                  currentSite%rootfrac_scr(j)  / &
                  rootfrac_notop
          end if
       end do

       ! Calculate the mean soil moisture ( liquid volume (m3/m3) and matric potential (mm))
       !    over the last 10 days
       mean_10day_liqvol = sum(currentSite%liqvol_memory(1:numWaterMem,ipft)) / &
                           real(numWaterMem,r8)
       mean_10day_smp    = sum(currentSite%smp_memory   (1:numWaterMem,ipft)) / &
                           real(numWaterMem,r8)

       ! Compare the moisture with the threshold.
       if ( phen_drought_threshold >= 0. ) then
          ! Liquid volume in reference layer (m3/m3)
          smoist_below_threshold = mean_10day_liqvol < phen_drought_threshold
       else
          ! Soil matric potential in reference layer (mm)
          smoist_below_threshold = mean_10day_smp    < phen_drought_threshold
       end if

       ! Calculate days since last flushing and shedding event, but make a provision
       ! for the first year of simulation, we have to assume leaf drop / leaf flush
       ! dates to start, so if that is in the future, set it to last year
       if (model_day_int < currentSite%dleafoffdate(ipft)) then
          currentSite%dndaysleafoff(ipft) = model_day_int - (currentSite%dleafoffdate(ipft)-ndays_per_year)
       else
          currentSite%dndaysleafoff(ipft) = model_day_int - currentSite%dleafoffdate(ipft)
       end if
       if (model_day_int < currentSite%dleafondate(ipft)) then
          currentSite%dndaysleafon(ipft) = model_day_int - (currentSite%dleafondate(ipft)-ndays_per_year)
       else
          currentSite%dndaysleafon(ipft) = model_day_int - currentSite%dleafondate(ipft)
       end if


       ! Elongation factor from the previous step.
       elongf_prev = currentSite%elong_factor(ipft)


       ! PFT leaf lifespan in days. This is the shortest between the leaf longevity
       ! (defined as a PFT parameter) and the maximum canopy leaf life span allowed
       ! for drought deciduous (local parameter). The sum term accounts for the
       ! total leaf life span of this cohort.
       ! Note we only use canopy leaf lifespan here and assume  that understory cohorts
       ! would  behave the same as canopy cohorts with regards to phenology. 
       ndays_pft_leaf_lifespan = &
          nint(ndays_per_year*min(decid_leaf_long_max,sum(prt_params%leaf_long(ipft,:))))


       !---~---
       !    Find elongation factors by comparing the moisture with the thresholds. For each
       ! tissue --- leaves, fine roots, and stems (sapwood+heartwood) --- elongation factor
       ! is the maximum fraction of biomass (relative to maximum biomass given allometry)
       ! that can be allocated to each tissue due to phenology. In this select case, we
       ! define the elongation factor based on the PFT-specific phenology strategy of this
       ! each PFT. Options are evergreen, "hard deciduous", or semi-deciduous:
       !  - Evergreen: elongation factors shall be 1 at all times (fully flushed tissues).
       !  - "Hard-deciduous": elongation factors are either 0 (fully abscised tissues) or
       !    1 (fully flushed tissues)
       !  - Semi-deciduous: elongation factors can be any value between 0 and 1 (including
       !    0 and 1). For example, if elongation factor for leaves of a cohort is 0.4, then
       !    the leaf biomass will be capped at 40% of the biomass the cohort would have if
       !    it were in well-watered conditions.
       !---~---
       case_drought_phen: select case (prt_params%stress_decid(ipft))
       case (ihard_stress_decid)
          !---~---
          !    Default ("hard") drought deciduous phenology. The decision on whether to 
          ! abscise (shed) or flush leaves is in principle defined by the soil moisture
          ! in the rooting zone.  However, we must also account the time since last 
          ! abscission or flushing event, to avoid excessive "flickering" of the leaf 
          ! elongation factor if soil moisture is right at the threshold.
          !
          ! (MLO thought: maybe we should define moisture equivalents of GDD and chilling
          ! days to simplify the cases a bit...)
          !---~---


          !---~---
          ! Save some conditions in logical variables to simplify code below
          !---~---
          ! Leaves have been "on" for longer than the minimum number of days.
          exceed_min_on_period     = &
             any( currentSite%dstatus(ipft) == [phen_dstat_timeon,phen_dstat_moiston] )   .and. &
             (currentSite%dndaysleafon(ipft) > dleafon_drycheck)
          ! Leaves have been "off" for longer than the minimum number of days.
          exceed_min_off_period    = &
             ( currentSite%dstatus(ipft)       == phen_dstat_timeoff       ) .and. &
             ( currentSite%dndaysleafoff(ipft) >  min_daysoff_dforcedflush )
          ! Leaves have been "on" for longer than the leaf lifetime.
          prolonged_on_period      = &
             any( currentSite%dstatus(ipft) == [phen_dstat_timeon,phen_dstat_moiston] )   .and. &
             ( currentSite%dndaysleafon(ipft) > ndays_pft_leaf_lifespan )
          ! Leaves have been "off" for a sufficiently long time and the last flushing
          ! was about one year ago (+/- tolerance).
          prolonged_off_period     = &
             any( currentSite%dstatus(ipft) == [phen_dstat_timeoff,phen_dstat_moistoff] ) .and. &
             ( currentSite%dndaysleafoff(ipft) > phen_doff_time     )                     .and. &
             ( currentSite%dndaysleafon(ipft) >= ndays_per_year-dd_offon_toler )                     .and. &
             ( currentSite%dndaysleafon(ipft) <= ndays_per_year+dd_offon_toler )
          ! Last flushing was a very long time ago.
          last_flush_long_ago      = &
             ( currentSite%dstatus(ipft)      == phen_dstat_moistoff            ) .and. &
             ( currentSite%dndaysleafon(ipft) >  ndays_per_year+dd_offon_toler  )
          !---~---


          !---~---
          ! Revision of the conditions, added an if/elseif/else structure to ensure only 
          ! up to one change occurs at any given time. Also, prevent changes until the
          ! soil moisture memory is populated (the outer if check).
          !---~---
          past_spinup_ifelse: if (model_day_int > numWaterMem) then
             drought_smoist_ifelse: if ( prolonged_off_period .and. &
                                         ( .not. smoist_below_threshold ) ) then
                ! LEAF ON: DROUGHT DECIDUOUS WETNESS
                ! Here, we used a window of oppurtunity to determine if we are
                ! close to the time when then leaves came on last year
                ! The following conditions must be met
                ! a) a year, plus or minus 1 month since we last had leaf-on?
                ! b) Has there also been at least a nominaly short amount of "leaf-off"?
                ! c) Is the soil moisture sufficiently high?
                currentSite%dstatus(ipft)      = phen_dstat_moiston  ! set status to leaf-on
                currentSite%dleafondate(ipft)  = model_day_int       ! save the model day we start flushing
                currentSite%dndaysleafon(ipft) = 0
                currentSite%elong_factor(ipft) = 1.

             elseif ( last_flush_long_ago ) then
                ! LEAF ON: DROUGHT DECIDUOUS TIME EXCEEDANCE
                ! If we still haven't done budburst by end of window, then force it

                ! If the status is "phen_dstat_moistoff", it means this site currently has
                ! leaves off due to actual moisture limitations.
                ! So we trigger bud-burst at the end of the month since
                ! last year's bud-burst.  If this is imposed, then we set the new
                ! status to indicate bud-burst was forced by timing
                currentSite%dstatus(ipft)      = phen_dstat_timeon ! force budburst!
                currentSite%dleafondate(ipft)  = model_day_int     ! record leaf on date
                currentSite%dndaysleafon(ipft) = 0
                currentSite%elong_factor(ipft) = 1.

             elseif ( exceed_min_off_period ) then
                ! LEAF ON: DROUGHT DECIDUOUS EXCEEDED MINIMUM OFF PERIOD
                ! Leaves were off due to time, not really moisture, so we allow them to
                ! flush again as soon as they exceed a minimum off time
                ! This typically occurs in a perennially wet system.
                currentSite%dstatus(ipft)      = phen_dstat_timeon    ! force budburst!
                currentSite%dleafondate(ipft)  = model_day_int        ! record leaf on date
                currentSite%dndaysleafon(ipft) = 0
                currentSite%elong_factor(ipft) = 1.

             elseif ( prolonged_on_period ) then
                ! LEAF OFF: DROUGHT DECIDUOUS LIFESPAN
                ! Are the leaves rouhgly at the end of their lives? If so, shed leaves 
                ! even if it is not dry.
                currentSite%dstatus(ipft)      = phen_dstat_timeoff    !alter status of site to 'leaves off'
                currentSite%dleafoffdate(ipft) = model_day_int         !record leaf on date
                currentSite%dndaysleafoff(ipft) = 0
                currentSite%elong_factor(ipft)  = 0.

             elseif ( exceed_min_on_period .and. smoist_below_threshold ) then
                ! LEAF OFF: DROUGHT DECIDUOUS DRYNESS - if the soil gets too dry,
                ! and the leaves have already been on a while...
                currentSite%dstatus(ipft) = phen_dstat_moistoff     ! alter status of site to 'leaves off'
                currentSite%dleafoffdate(ipft) = model_day_int      ! record leaf on date
                currentSite%dndaysleafoff(ipft) = 0
                currentSite%elong_factor(ipft)  = 0.
             end if drought_smoist_ifelse
          end if past_spinup_ifelse
          !---~---


       case (isemi_stress_decid)
          !---~---
          ! Semi-deciduous PFT, based on ED2.  We compare the moisture with the lower
          ! and upper thresholds. If the moisture is in between the thresholds, we must
          ! also check whether or not the drought is developing or regressing.
          !---~---


          !---~---
          !   First guess elongation factor, solely based on rooting-zone moisture.
          ! These values may be adjusted based on the time since last flushing and/or
          ! abscising event.
          !---~---
          if (phen_drought_threshold >= 0.) then
             elongf_1st = elongf_min + (1.0_r8 - elongf_min ) * &
                          ( mean_10day_liqvol    - phen_drought_threshold ) / &
                          ( phen_moist_threshold - phen_drought_threshold )
          else
             elongf_1st = elongf_min + (1.0_r8 - elongf_min ) * &
                          ( mean_10day_smp       - phen_drought_threshold ) / &
                          ( phen_moist_threshold - phen_drought_threshold )
          end if
          elongf_1st = max(0.0_r8,min(1.0_r8,elongf_1st))
          !---~---



          !---~---
          ! Save some conditions in logical variables to simplify code below
          !---~---
          !  Leaves have been flushing for a short period of time.
          recent_flush         = elongf_prev >= elongf_min .and. &
                                 ( currentSite%dndaysleafon(ipft) <= dleafon_drycheck )
          !  Leaves have been abscissing for a short period of time.
          recent_abscission    = elongf_prev <  elongf_min .and. &
                                 ( currentSite%dndaysleafoff(ipft) <=  min_daysoff_dforcedflush )
          !  Leaves have been flushing for longer than their time span.
          prolonged_on_period  = all( [elongf_prev,elongf_1st] >= elongf_min ) .and. &
                                 ( currentSite%dndaysleafon(ipft)  > ndays_pft_leaf_lifespan )
          !  It's been a long time since the plants had flushed their leaves.
          last_flush_long_ago  = all( [elongf_prev,elongf_1st] <  elongf_min ) .and. &
                                 ( currentSite%dndaysleafon(ipft) >  ndays_per_year+dd_offon_toler )
          !---~---


          ! Make sure elongation factor is bounded and check for special cases.
          drought_gradual_ifelse: if ( model_day_int <= numWaterMem ) then
             ! Too early in the simulation, keep the same elongation factor as the day before.
             currentSite%elong_factor(ipft) = elongf_prev

          elseif ( prolonged_on_period ) then
             ! Leaves have been on for too long and exceeded leaf lifespan. Force abscission
             currentSite%elong_factor(ipft)  = 0.0_r8             ! Force full budburst
             currentSite%dstatus(ipft)       = phen_dstat_timeoff ! Flag that this has been forced
             currentSite%dleafoffdate(ipft)  = model_day_int      ! Record leaf off date
             currentSite%dndaysleafoff(ipft) = 0                  ! Reset clock

          elseif ( last_flush_long_ago ) then
             ! Plant has not flushed at all for a very long time. Force flushing
             currentSite%elong_factor(ipft)  = elongf_min         ! Force minimum budburst
             currentSite%dstatus(ipft)       = phen_dstat_timeon  ! Flag that this has been forced
             currentSite%dleafondate(ipft)   = model_day_int      ! Record leaf on date
             currentSite%dndaysleafon(ipft)  = 0                  ! Reset clock

          elseif ( recent_flush .and. elongf_1st < elongf_prev ) then
             ! Leaves have only recently reached flushed status. Elongation factor cannot decrease
             currentSite%elong_factor(ipft) = elongf_prev       ! Elongation factor cannot decrease
             currentSite%dstatus(ipft)      = phen_dstat_timeon ! Flag that this has been forced

          elseif ( recent_abscission .and. elongf_1st > elongf_min ) then
             ! Leaves have only recently abscissed. Prevent plant to flush leaves.
             currentSite%elong_factor(ipft) = 0.0_r8             ! Elongation factor must remain 0.
             currentSite%dstatus(ipft)      = phen_dstat_timeoff ! Flag that this has been forced

          elseif ( elongf_1st < elongf_min ) then
             ! First guess of elongation factor below minimum. Impose full abscission.
             currentSite%elong_factor(ipft) = 0.0_r8

             if (elongf_prev >= elongf_min ) then
                ! This is the first day moisture fell below minimum. Flag change of status.
                currentSite%dstatus(ipft)       = phen_dstat_moistoff ! Flag that this has not been forced
                currentSite%dleafoffdate(ipft)  = model_day_int       ! Record leaf off date
                currentSite%dndaysleafoff(ipft) = 0                   ! Reset clock
             end if

          else
             ! First guess of elongation factor is valid, use it.
             currentSite%elong_factor(ipft) = elongf_1st


             if (elongf_prev < elongf_min ) then
                ! This is the first day moisture allows leaves to exist. Flag change of status.
                currentSite%dstatus(ipft)       = phen_dstat_moiston  ! Flag that this has not been forced
                currentSite%dleafondate(ipft)   = model_day_int       ! Record leaf on date
                currentSite%dndaysleafon(ipft)  = 0                   ! Reset clock
             elseif (elongf_1st < elongf_prev) then
                currentSite%dstatus(ipft)       = phen_dstat_pshed    ! Flag partial shedding,
                                                                    ! but do not reset the clock
             end if
          end if drought_gradual_ifelse


       case default
          !    Neither hard deciduous or semi-deciduous. For now we treat this as synonym
          ! of non-drought deciduous. In the future we may consider other drought deciduous
          ! strategies (e.g., abscission driven by moisture, flushing driven by photo-
          ! period).
          currentSite%dstatus(ipft)      = phen_dstat_moiston

          ! Assign elongation factors for non-drought deciduous PFTs, which will be used
          ! to define the cohort status.
          case_cold_phen: select case(prt_params%season_decid(ipft))
          case (ifalse)
             ! Evergreen, ensure that elongation factor is always one.
             currentSite%elong_factor(ipft) = 1.0_r8
          case (itrue)
             ! Cold-deciduous. Define elongation factor based on cold status
             select case (currentSite%cstatus)
             case (phen_cstat_nevercold,phen_cstat_iscold)
                currentSite%elong_factor(ipft) = 0.0_r8
             case (phen_cstat_notcold)
                currentSite%elong_factor(ipft) = 1.0_r8
             end select
          end select case_cold_phen

       end select case_drought_phen

    end do pft_elong_loop

    call phenology_leafonoff(currentSite)

    return
  end subroutine phenology


  ! ============================================================================
  subroutine phenology_leafonoff(currentSite)
    !
    ! !DESCRIPTION:
    ! Controls the leaf on and off economics
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type) , pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    real(r8) :: leaf_c                   ! leaf carbon [kg]
    real(r8) :: fnrt_c                   ! fine root carbon [kg]
    real(r8) :: sapw_c                   ! sapwood carbon [kg]
    real(r8) :: struct_c                 ! structural wood carbon [kg]
    real(r8) :: store_c                  ! storage carbon [kg]
    real(r8) :: store_c_transfer_frac    ! Fraction of storage carbon used to flush leaves

    real(r8) :: leaf_deficit_c           ! leaf carbon deficit (relative to target) [kg]
    real(r8) :: fnrt_deficit_c           ! fine root carbon deficit (relative to target) [kg]
    real(r8) :: sapw_deficit_c           ! sapwood carbon deficit (relative to target) [kg]
    real(r8) :: struct_deficit_c         ! structural wood carbon deficit (relative to target) [kg]
    real(r8) :: total_deficit_c          ! total carbon deficit (relative to target) [kg]

    real(r8) :: target_leaf_c            ! target leaf carbon (allometry scaled by elongation factor) [kg]
    real(r8) :: target_fnrt_c            ! target fine root carbon (allometry scaled by elongation factor) [kg]
    real(r8) :: target_sapw_c            ! target sapwood carbon (allometry scaled by elongation factor) [kg]
    real(r8) :: target_agw_c             ! target Above ground biomass [kgC]
    real(r8) :: target_bgw_c             ! target Below ground biomass [kgC]
    real(r8) :: target_struct_c          ! target structural wood carbon (allometry scaled by elongation factor) [kg]

    real(r8) :: sapw_area                ! Sapwood area

    real(r8) :: eff_leaf_drop_fraction   ! Effective leaf drop fraction
    real(r8) :: eff_fnrt_drop_fraction   ! Effective fine-root drop fraction
    real(r8) :: eff_sapw_drop_fraction   ! Effective sapwood drop fraction
    real(r8) :: eff_struct_drop_fraction ! Effective structural wood drop fraction

    logical  :: is_flushing_time         ! Time to flush leaves
    logical  :: is_shedding_time         ! Time to shed leaves

    real(r8) :: fnrt_drop_fraction       ! Fine root relative drop fraction (0 = no drop, 1 = as much as leaves)
    real(r8) :: stem_drop_fraction       ! Stem drop relative fraction (0 = no drop, 1 = as much as leaves)
    real(r8) :: l2fr                     ! Leaf to fineroot biomass multiplier 

    integer  :: ipft                     ! Plant functional type index
    real(r8), parameter :: leaf_drop_fraction  = 1.0_r8
    real(r8), parameter :: carbon_store_buffer = 0.10_r8
    !------------------------------------------------------------------------

    currentPatch => CurrentSite%oldest_patch

    patch_loop: do while(associated(currentPatch))
       currentCohort => currentPatch%tallest
       cohort_loop: do while(associated(currentCohort))

          ipft = currentCohort%pft

          ! Retrieve existing leaf and storage carbon

          if(debug) call currentCohort%prt%CheckMassConservation(ipft,0)

          store_c  = currentCohort%prt%GetState(store_organ , carbon12_element)
          leaf_c   = currentCohort%prt%GetState(leaf_organ  , carbon12_element)
          fnrt_c   = currentCohort%prt%GetState(fnrt_organ  , carbon12_element)
          sapw_c   = currentCohort%prt%GetState(sapw_organ  , carbon12_element)
          struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)

          fnrt_drop_fraction = prt_params%phen_fnrt_drop_fraction(ipft)
          stem_drop_fraction = prt_params%phen_stem_drop_fraction(ipft)
          l2fr               = prt_params%allom_l2fr(ipft)

          ! MLO. To avoid duplicating code for drought and cold deciduous PFTs, we first
          !      check whether or not it's time to flush or time to shed leaves, then
          !      use a common code for flushing or shedding leaves.
          is_time_block: if (prt_params%season_decid(ipft) == itrue) then ! Cold deciduous

             ! A. Is this the time for COLD LEAVES to switch to ON?
             is_flushing_time = ( currentSite%cstatus      == phen_cstat_notcold .and. & ! We just moved to leaves being on
                                  currentCohort%status_coh == leaves_off         )        ! Leaves are currently off
             ! B. Is this the time for COLD LEAVES to switch to OFF?
             is_shedding_time = any(currentSite%cstatus == [phen_cstat_nevercold,phen_cstat_iscold]) .and. & ! Past leaf drop day or too cold
                                currentCohort%status_coh == leaves_on                                .and. & ! Leaves have not dropped yet
                                ( currentCohort%dbh > EDPftvarcon_inst%phen_cold_size_threshold(ipft) .or. & ! Grasses are big enough or...
                                  prt_params%woody(ipft) == itrue                                     )      ! this is a woody PFT.

          elseif (any(prt_params%stress_decid(ipft) == [ihard_stress_decid,isemi_stress_decid]) ) then ! Drought deciduous

             ! A. Is this the time for DROUGHT LEAVES to switch to ON?
             is_flushing_time = any( currentSite%dstatus(ipft) == [phen_dstat_moiston,phen_dstat_timeon] ) .and.  & ! Leaf flushing time (moisture or time)
                                any( currentCohort%status_coh  == [leaves_off,leaves_shedding] )
             ! B. Is this the time for DROUGHT LEAVES to switch to OFF?
             !    This will be true when leaves are abscissing (partially or fully) due to moisture or time
             is_shedding_time = any( currentSite%dstatus(ipft) == [phen_dstat_moistoff,phen_dstat_timeoff,phen_dstat_pshed] ) .and. &
                                any( currentCohort%status_coh  == [leaves_on,leaves_shedding] )
          else
             ! This PFT is not deciduous.
             is_flushing_time         = .false.
             is_shedding_time         = .false.
          end if is_time_block



          ! Elongation factor for leaves is always the same as the site- and 
          ! PFT-dependent factor computed in subroutine phenology. For evergreen
          ! PFTs, this value should be always 1.0. 
          currentCohort%efleaf_coh = currentSite%elong_factor(ipft)

          ! Find the effective "elongation factor" for fine roots and stems. The effective elongation
          ! factor is a combination of the PFT leaf elongation factor (efleaf_coh) and the tissue drop 
          ! fraction relative to leaves (xxxx_drop_fraction). When xxxx_drop_fraction is 0, the biomass
          ! of tissue xxxx will not be impacted by phenology. If xxxx_drop_fraction is 1, the biomass 
          ! of tissue xxxx will be as impacted by phenology as leaf biomass. Intermediate values will
          ! allow a more moderate impact of phenology in tissue xxxx relative to leaves.
          currentCohort%effnrt_coh = 1.0_r8 - (1.0_r8 - currentCohort%efleaf_coh ) * fnrt_drop_fraction
          currentCohort%efstem_coh = 1.0_r8 - (1.0_r8 - currentCohort%efleaf_coh ) * stem_drop_fraction

          ! Find the target biomass for each tissue  when accounting for elongation
          ! factors. Note that the target works for both flushing and shedding leaves.
          call bleaf(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage, &
               currentCohort%canopy_trim,currentCohort%efleaf_coh,target_leaf_c)
          call bfineroot(currentCohort%dbh,currentCohort%pft, &
               currentCohort%canopy_trim,l2fr,currentCohort%effnrt_coh,target_fnrt_c)
          call bsap_allom(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage, &
               currentCohort%canopy_trim,currentCohort%efstem_coh,sapw_area,target_sapw_c)
          call bagw_allom(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage,&
               currentCohort%efstem_coh,target_agw_c)
          call bbgw_allom(currentCohort%dbh,currentCohort%pft,currentCohort%efstem_coh,target_bgw_c)
          call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, &
               currentCohort%pft, target_struct_c)


          ! A.  This is time to switch to (COLD or DROUGHT) LEAF ON
          flush_block: if (is_flushing_time) then
             currentCohort%status_coh = leaves_on ! Leaves are on, so change status to
                                                  ! stop flow of carbon out of bstore.

             ! Transfer carbon from storage to living tissues (only if there is any carbon in storage)
             transf_block: if ( store_c > nearzero ) then
                ! Find the total deficit.  We no longer distinguish between woody and non-woody
                ! PFTs here (as sapwmemory is the same as sapw_c if this is a woody tissue).
                leaf_deficit_c   = max(0.0_r8, target_leaf_c   - leaf_c  )
                fnrt_deficit_c   = max(0.0_r8, target_fnrt_c   - fnrt_c  )
                sapw_deficit_c   = max(0.0_r8, target_sapw_c   - sapw_c  )
                struct_deficit_c = max(0.0_r8, target_struct_c - struct_c)
                total_deficit_c  = leaf_deficit_c + fnrt_deficit_c + sapw_deficit_c + &
                                   struct_deficit_c

                ! Flush either the amount required from the memory, or -most- of the storage pool
                ! RF: added a criterion to stop the entire store pool emptying and triggering termination mortality
                ! n.b. this might not be necessary if we adopted a more gradual approach to leaf flushing...
                store_c_transfer_frac = min( EDPftvarcon_inst%phenflush_fraction(ipft) * &
                                             total_deficit_c / store_c, &
                                             1.0_r8 - carbon_store_buffer )

                ! This call will request that storage carbon will be transferred to
                ! each tissue. It is specified as a fraction of the available storage
                ! MLO - Just to be safe, skip steps in the unlikely case total_deficit is zero, to avoid FPE errors.
                if (total_deficit_c > nearzero) then
                   call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                                          store_c_transfer_frac*leaf_deficit_c/total_deficit_c)
                   call PRTPhenologyFlush(currentCohort%prt, ipft, fnrt_organ, &
                                          store_c_transfer_frac*fnrt_deficit_c/total_deficit_c)

                   ! MLO - stem_drop_fraction is a PFT parameter, do we really need this 
                   !       check for woody/non-woody PFT?
                   if ( prt_params%woody(ipft) == ifalse ) then
                      call PRTPhenologyFlush(currentCohort%prt, ipft, sapw_organ, &
                                             store_c_transfer_frac*sapw_deficit_c/total_deficit_c)
                      call PRTPhenologyFlush(currentCohort%prt, ipft, struct_organ, &
                                             store_c_transfer_frac*struct_deficit_c/total_deficit_c)
                   end if
                end if
             else
                ! Not enough carbon to flush any living tissue.
                store_c_transfer_frac = 0.0_r8
             end if transf_block
          end if flush_block



          ! B.  This is time to switch to (COLD or DROUGHT) LEAF OFF
          shed_block: if (is_shedding_time) then
             if ( currentCohort%efleaf_coh > 0.0_r8 ) then
                ! Partial shedding
                currentCohort%status_coh  = leaves_shedding
             else
                ! Complete abscission
                currentCohort%status_coh  = leaves_off
             end if


             ! Find the effective fraction to drop. This fraction must be calculated every time
             ! because we must account for partial abscission. The simplest approach is to simply
             ! use the ratio between the target and the original biomass of each pool. The 
             ! max(tissue_c,nearzero) is overly cautious, because leaf_c = 0 would imply that
             ! leaves are already off, and this wouldn't be considered shedding time.
             eff_leaf_drop_fraction   = max( 0.0_r8, min( 1.0_r8,1.0_r8 - target_leaf_c   / max( leaf_c  , nearzero ) ) )
             eff_fnrt_drop_fraction   = max( 0.0_r8, min( 1.0_r8,1.0_r8 - target_fnrt_c   / max( fnrt_c  , nearzero ) ) )
             eff_sapw_drop_fraction   = max( 0.0_r8, min( 1.0_r8,1.0_r8 - target_sapw_c   / max( sapw_c  , nearzero ) ) )
             eff_struct_drop_fraction = max( 0.0_r8, min( 1.0_r8,1.0_r8 - target_struct_c / max( struct_c, nearzero ) ) )

             ! Drop leaves
             call PRTDeciduousTurnover(currentCohort%prt,ipft, leaf_organ, eff_leaf_drop_fraction)

             ! Drop fine roots
             call PRTDeciduousTurnover(currentCohort%prt,ipft, fnrt_organ, eff_fnrt_drop_fraction)

             ! If plant is not woody, shed sapwood and heartwood (they may have a minimum amount of woody tissues for
             ! running plant hydraulics, and it makes sense to shed them along with leaves when they should be off).
             ! MLO - stem_drop_fraction is a PFT parameter, do we really need this check for woody/non-woody PFT?
             if ( prt_params%woody(ipft) == ifalse ) then
                ! Shed sapwood and heartwood.
                call PRTDeciduousTurnover(currentCohort%prt,ipft,sapw_organ  , eff_sapw_drop_fraction  )
                call PRTDeciduousTurnover(currentCohort%prt,ipft,struct_organ, eff_struct_drop_fraction)
             end if

          end if shed_block

          if(debug) call currentCohort%prt%CheckMassConservation(ipft,1)

          currentCohort => currentCohort%shorter
       end do cohort_loop

       currentPatch => currentPatch%younger

    end do patch_loop

  end subroutine phenology_leafonoff

  ! =====================================================================================

  subroutine satellite_phenology(currentSite, bc_in)

    ! -----------------------------------------------------------------------------------
    ! Takes the daily inputs of leaf area index, stem area index and canopy height and
    ! translates them into a FATES structure with one patch and one cohort per PFT
    ! The leaf area of the cohort is modified each day to match that asserted by the HLM
    ! -----------------------------------------------------------------------------------

    ! !USES:
    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    type(bc_in_type),   intent(in)            :: bc_in

    class(prt_vartypes), pointer :: prt

    ! !LOCAL VARIABLES:
    type(fates_patch_type) , pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort

    real(r8) ::  spread        ! dummy value of canopy spread to estimate c_area
    real(r8) ::  leaf_c        ! leaf carbon estimated to generate target tlai
    real(r8) :: check_treelai
    integer ::   fates_pft     ! fates pft numer for weighting loop
    integer  ::   hlm_pft      ! host land model pft number for weighting loop.
    integer ::   s             ! site index


    ! To Do in this routine.
    ! Get access to HLM input varialbes.
    ! Weight them by PFT
    ! Loop around patches, and for each single cohort in each patch
    ! call assign_cohort_SP_properties to determine cohort height, dbh, 'n',  area, leafc from drivers.

    currentSite%sp_tlai(:) = 0._r8
    currentSite%sp_tsai(:) = 0._r8
    currentSite%sp_htop(:) = 0._r8

    ! WEIGHTING OF FATES PFTs on to HLM_PFTs
    ! 1. Add up the area associated with each FATES PFT
    ! where pft_areafrac is the area of land in each HLM PFT and (from surface dataset)
    ! hlm_pft_map is the area of that land in each FATES PFT (from param file)

    ! 2. weight each fates PFT target for lai, sai and htop by the area of the
    ! contrbuting HLM PFTs.

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))

       fates_pft = currentPatch%nocomp_pft_label
       if(fates_pft.ne.0)then

          do hlm_pft = 1,size( EDPftvarcon_inst%hlm_pft_map,2)

             if(bc_in%pft_areafrac(hlm_pft) * EDPftvarcon_inst%hlm_pft_map(fates_pft,hlm_pft).gt.0.0_r8)then
                !leaf area index
                currentSite%sp_tlai(fates_pft) = currentSite%sp_tlai(fates_pft) + &
                     bc_in%hlm_sp_tlai(hlm_pft) * bc_in%pft_areafrac(hlm_pft) &
                     * EDPftvarcon_inst%hlm_pft_map(fates_pft,hlm_pft)
                !stem area index
                currentSite%sp_tsai(fates_pft) = currentSite%sp_tsai(fates_pft) + &
                     bc_in%hlm_sp_tsai(hlm_pft) *	bc_in%pft_areafrac(hlm_pft) &
                     * EDPftvarcon_inst%hlm_pft_map(fates_pft,hlm_pft)
                ! canopy height
                currentSite%sp_htop(fates_pft) = currentSite%sp_htop(fates_pft) + &
                     bc_in%hlm_sp_htop(hlm_pft) * bc_in%pft_areafrac(hlm_pft) &
                     * EDPftvarcon_inst%hlm_pft_map(fates_pft,hlm_pft)
             end if ! there is some area in this patch
          end do !hlm_pft

          ! weight for total area in each patch/fates_pft
          ! this is needed because the area of pft_areafrac does not need to sum to 1.0
          if(currentPatch%area.gt.0.0_r8)then
             currentSite%sp_tlai(fates_pft) = currentSite%sp_tlai(fates_pft) &
                  /(currentPatch%area/area)
             currentSite%sp_tsai(fates_pft) = currentSite%sp_tsai(fates_pft) &
                  /(currentPatch%area/area)
             currentSite%sp_htop(fates_pft) = currentSite%sp_htop(fates_pft) &
                  /(currentPatch%area/area)
          endif

       end if ! not bare patch
       currentPatch => currentPatch%younger
    end do ! patch loop

    ! ------------------------------------------------------------
    ! now we have the target lai, sai and htop for each PFT/patch
    ! find properties of the cohort that go along with that
    ! 1. Find canopy area from HTOP (height)
    ! 2. Find 'n' associated with canopy area, given a closed canopy
    ! 3. Find 'bleaf' associated with TLAI and canopy area.
    ! These things happen in  the catchily titled "assign_cohort_SP_properties" routine.
    ! ------------------------------------------------------------

    currentPatch => currentSite%oldest_patch
    do while (associated(currentPatch))

       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))

          ! FIRST SOME CHECKS.
          fates_pft =currentCohort%pft
          if(fates_pft.ne.currentPatch%nocomp_pft_label)then ! does this cohort belong in this PFT patch?
             write(fates_log(),*) 'wrong PFT label in cohort in SP mode',fates_pft,currentPatch%nocomp_pft_label
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          if(fates_pft.eq.0)then
             write(fates_log(),*) 'PFT0 in SP mode'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          ! Call routine to invert SP drivers into cohort properites.
          call assign_cohort_SP_properties(currentCohort, currentSite%sp_htop(fates_pft), currentSite%sp_tlai(fates_pft)     , currentSite%sp_tsai(fates_pft),currentPatch%area,ifalse,leaf_c)

          currentCohort => currentCohort%shorter
       end do !cohort loop
       currentPatch => currentPatch%younger
    end do ! patch loop

  end subroutine satellite_phenology

  ! ======================================================================================

  subroutine calculate_SP_properties(htop, tlai, tsai, parea, pft, crown_damage,         &
      canopy_layer, vcmax25top, leaf_c, dbh, cohort_n, c_area)
    !
    ! DESCRIPTION:
    !  Takes the daily inputs of leaf area index, stem area index and canopy height and
    !  translates them into a FATES structure with one patch and one cohort per PFT.
    !  The leaf area of the cohort is modified each day to match that asserted by the HLM
    !

    ! ARGUMENTS:
    real(r8), intent(in)  :: tlai         ! target leaf area index from SP inputs [m2 m-2]
    real(r8), intent(in)  :: tsai         ! target stem area index from SP inputs [m2 m-2]
    real(r8), intent(in)  :: htop         ! target tree height from SP inputs [m]
    real(r8), intent(in)  :: parea        ! patch area for this PFT [m2]
    real(r8), intent(in)  :: vcmax25top   ! maximum carboxylation at canopy top and 25degC [umol CO2/m2/s]
    integer,  intent(in)  :: pft          ! cohort PFT index
    integer,  intent(in)  :: crown_damage ! cohort crown damage status
    integer,  intent(in)  :: canopy_layer ! canopy status of cohort [1 = canopy, 2 = understorey, etc.]
    real(r8), intent(out) :: leaf_c       ! leaf carbon estimated to generate target tlai [kgC]
    real(r8), intent(out) :: dbh          ! cohort diameter at breast height [cm]
    real(r8), intent(out) :: cohort_n     ! cohort density [/m2]
    real(r8), intent(out) :: c_area

    ! LOCAL VARIABLES:
    real(r8) :: check_treelai       ! check tree LAI against input tlai [m2/m2]
    real(r8) :: dummy_treesai       ! dummy
    real(r8) :: canopylai(1:nclmax) ! canopy LAI [m2/m2]
    real(r8) :: oldcarea            ! save value of crown area [m2]

    ! calculate DBH from input height
    call h2d_allom(htop, pft, dbh)

    ! calculate canopy area, assuming n = 1.0 and spread = 1.0_r8
    call carea_allom(dbh, 1.0_r8, 1.0_r8, pft, crown_damage, c_area)

    ! calculate canopy N assuming patch area is full
    cohort_n = parea/c_area

    ! correct c_area for the new nplant, assuming spread = 1.0
    call carea_allom(dbh, cohort_n, 1.0_r8, pft, crown_damage, c_area)

    ! calculate leaf carbon from target treelai
    canopylai(:) = 0._r8
    leaf_c = leafc_from_treelai(tlai, tsai, pft, c_area, cohort_n, canopy_layer, vcmax25top)
    
    ! check that the inverse calculation of leafc from treelai is the same as the
    ! standard calculation of treelai from leafc. Maybe can delete eventually?

    call tree_lai_sai(leaf_c, pft, c_area, cohort_n, canopy_layer, canopylai, vcmax25top, &
                          dbh, crown_damage, 1.0_r8, 1.0_r8, 11, check_treelai, dummy_treesai)
    
    if (abs(tlai - check_treelai) > area_error_2) then !this is not as precise as nearzero
      write(fates_log(),*) 'error in validate treelai', tlai, check_treelai, tlai - check_treelai
      write(fates_log(),*) 'tree_lai inputs: ', pft, c_area, cohort_n, canopy_layer, vcmax25top
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! the carea_allom routine sometimes generates precision-tolerance level errors in the canopy area
    ! these mean that the canopy area does not exactly add up to the patch area, which causes chaos in
    ! the radiation routines.  Correct both the area and the 'n' to remove error, and don't use
    ! carea_allom in SP mode after this point.

    if (abs(c_area - parea) > nearzero) then ! there is an error
      if (abs(c_area - parea) < area_error_3) then ! correct this if it's a very small error
          oldcarea = c_area
          ! generate new cohort area
          c_area = c_area - (c_area - parea)
          cohort_n = cohort_n*(c_area/oldcarea)
          if (abs(c_area-parea) > nearzero) then
            write(fates_log(),*) 'SPassign, c_area still broken', c_area - parea, c_area - oldcarea
            call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       else
          write(fates_log(),*) 'SPassign, big error in c_area', c_area - parea, pft
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if ! still broken
    end if !small error

  end subroutine calculate_SP_properties

  ! ======================================================================================

  subroutine assign_cohort_SP_properties(currentCohort, htop, tlai, tsai, parea, init,   &
    leaf_c)
    !
    ! DESCRIPTION:
    !  Takes the daily inputs of leaf area index, stem area index and canopy height and
    !  translates them into a FATES structure with one patch and one cohort per PFT.
    !  The leaf area of the cohort is modified each day to match that asserted by the HLM

   
    ! ARGUMENTS
    type(fates_cohort_type), intent(inout), target :: currentCohort ! cohort object
    real(r8),                intent(in)            :: tlai          ! target leaf area index from SP inputs [m2/m2]
    real(r8),                intent(in)            :: tsai          ! target stem area index from SP inputs [m2/m2]
    real(r8),                intent(in)            :: htop          ! target tree height from SP inputs [m]
    real(r8),                intent(in)            :: parea         ! patch area for this PFT [m2]
    integer,                 intent(in)            :: init          ! are we in the initialization routine? if so do not set leaf_c
    real(r8),                intent(out)           :: leaf_c       ! leaf carbon estimated to generate target tlai [kgC]

    ! LOCAL VARIABLES
    real(r8) :: dbh      ! cohort dbh [cm]
    real(r8) :: cohort_n ! cohort density [/m2]
    real(r8) :: c_area   ! cohort canopy area [m2]

    if (associated(currentCohort%shorter)) then
      write(fates_log(),*) 'SP mode has >1 cohort'
      write(fates_log(),*) "SP mode >1 cohort: PFT", currentCohort%pft, currentCohort%shorter%pft
      write(fates_log(),*) "SP mode >1 cohort: CL", currentCohort%canopy_layer, currentCohort%shorter%canopy_layer
      call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    if (init .eq. itrue) then
      ! If we are initializing, the canopy layer has not been set yet, so just set to 1
      currentCohort%canopy_layer = 1
      ! We need to get the vcmax25top
      currentCohort%vcmax25top = EDPftvarcon_inst%vcmax25top(currentCohort%pft, 1)
    endif

    call calculate_SP_properties(htop, tlai, tsai, parea, currentCohort%pft,             &
      currentCohort%crowndamage, currentCohort%canopy_layer, currentCohort%vcmax25top,   &
      leaf_c, dbh, cohort_n, c_area)

    ! set allometric characteristics
    currentCohort%height = htop
    currentCohort%dbh = dbh
    currentCohort%n = cohort_n
    currentCohort%c_area = c_area
    currentCohort%treelai = tlai
    currentCohort%treesai = tsai

    if (init .eq. ifalse) then
      call SetState(currentCohort%prt, leaf_organ, carbon12_element, leaf_c, 1)
    endif

  end subroutine assign_cohort_SP_properties

  ! =====================================================================================
  
  subroutine SeedUpdate( currentSite )

    ! -----------------------------------------------------------------------------------
    ! Flux from plants into the seed pool.
    ! It is assumed that allocation to seed on living pools has already been calculated
    ! at the daily time step.
    ! Note: Some seed generation can occur during disturbance. It is assumed that
    !       some plants use their storage upon death to create seeds, but this in only
    !       triggered during non-fire and non-logging events.  See
    !       subroutine mortality_litter_fluxes() and DistributeSeeds(), look for
    !       parameter allom_frbstor_repro
    ! -----------------------------------------------------------------------------------


    ! !USES:
    use EDTypesMod, only : area
    use EDTypesMod, only : homogenize_seed_pfts
    use FatesInterfaceTypesMod, only : hlm_seeddisp_cadence
    use FatesInterfaceTypesMod, only : fates_dispersal_cadence_none
    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout), target  :: currentSite

    type(fates_patch_type), pointer     :: currentPatch
    type(litter_type), pointer       :: litt
    type(fates_cohort_type), pointer    :: currentCohort
    type(site_massbal_type), pointer :: site_mass

    integer  :: pft
    real(r8) :: store_m_to_repro       ! mass sent from storage to reproduction upon death [kg/plant]
    real(r8) :: site_seed_rain(numpft) ! This is the sum of seed-rain for the site [kg/site/day]
    real(r8) :: site_disp_frac(numpft) ! Fraction of seeds from prodeced in current grid cell to
                                       ! disperse out to other gridcells
    real(r8) :: seed_in_external       ! Mass of externally generated seeds [kg/m2/day]
    real(r8) :: seed_stoich            ! Mass ratio of nutrient per C12 in seeds [kg/kg]
    real(r8) :: seed_prod              ! Seed produced in this dynamics step [kg/day]
    integer  :: n_litt_types           ! number of litter element types (c,n,p, etc)
    integer  :: el                     ! loop counter for litter element types
    integer  :: element_id             ! element id consistent with parteh/PRTGenericMod.F90

    ! If the dispersal kernel is not turned on, keep the dispersal fraction at zero
    site_disp_frac(:) = 0._r8
    if (hlm_seeddisp_cadence .ne. fates_dispersal_cadence_none) then
      site_disp_frac(:) = EDPftvarcon_inst%seed_dispersal_fraction(:)
    end if

    el_loop: do el = 1, num_elements

       site_seed_rain(:) = 0._r8
       element_id = element_list(el)

       site_mass => currentSite%mass_balance(el)

       ! Loop over all patches and sum up the seed input for each PFT
       currentPatch => currentSite%oldest_patch
       seed_rain_loop: do while (associated(currentPatch))

          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))

             pft = currentCohort%pft

             ! a certain fraction of bstore might go to clonal reproduction when plants die
             ! (since this is only applied to the dying portion of the cohort
             !  we do not actually pair down the storage via PARTEH, instead
             !  we just make sure that we don't send a portion of the storage
             !  to the litter in CWDInput)
             ! units = [kg/ha/day] = [kg] * [fraction] * [plants/ha/year] * [year/day]
             store_m_to_repro = -currentCohort%prt%GetState(store_organ,element_id) * &
                  EDPftvarcon_inst%allom_frbstor_repro(pft)*currentCohort%dndt*years_per_day

             ! Transfer all reproductive tissues into seed production
             ! The following call to PRTReproRelease, will return the mass
             ! of seeds [kg] released by the plant, per the mass_fraction
             ! specified as input.  This routine will also remove the mass
             ! from the parteh state-variable.
             
             call PRTReproRelease(currentCohort%prt,repro_organ,element_id, &
                  1.0_r8, seed_prod)

             if(element_id==carbon12_element)then
                currentcohort%seed_prod = seed_prod
             end if


             site_seed_rain(pft) = site_seed_rain(pft) +  &
                  (seed_prod * currentCohort%n + store_m_to_repro) ![kg/site/day, kg/ha/day]

             currentCohort => currentCohort%shorter
          enddo !cohort loop

          currentPatch => currentPatch%younger
       enddo seed_rain_loop

       ! We can choose to homogenize seeds. This is simple, we just
       ! add up all the seed from each pft at the site level, and then
       ! equally distribute to the PFT pools
       if ( homogenize_seed_pfts ) then
          site_seed_rain(1:numpft) = sum(site_seed_rain(:))/real(numpft,r8)
       end if

       ! Loop over all patches again and disperse the mixed seeds into the input flux
       ! arrays
       ! Loop over all patches and sum up the seed input for each PFT
       currentPatch => currentSite%oldest_patch
       seed_in_loop: do while (associated(currentPatch))

          litt => currentPatch%litter(el)
          do pft = 1,numpft

             if(currentSite%use_this_pft(pft).eq.itrue)then

                ! Seed input from local sources (within site).  Note that a fraction of the
                ! internal seed rain is sent out to neighboring gridcells.
                litt%seed_in_local(pft) = litt%seed_in_local(pft) + site_seed_rain(pft)*(1.0_r8-site_disp_frac(pft))/area ![kg/m2/day]

                ! If we are using the Tree Recruitment Scheme (TRS) with or w/o seedling dynamics
                if ( any(regeneration_model == [TRS_regeneration, TRS_no_seedling_dyn]) .and. &
                     prt_params%allom_dbh_maxheight(pft) > min_max_dbh_for_trees) then
                   
                   ! Send a fraction of reproductive carbon to litter to account for 
                   ! non-seed reproductive carbon (e.g. flowers, fruit, etc.)
                   litt%seed_decay(pft) = litt%seed_in_local(pft) * (1.0_r8 - EDPftvarcon_inst%repro_frac_seed(pft)) 
                   
                   ! Note: The default regeneration scheme sends all reproductive carbon to seed
                end if !Use TRS
                
                ! If there is forced external seed rain, we calculate the input mass flux
                ! from the different elements, using the mean stoichiometry of new
                ! recruits for the current patch and lowest canopy position
                
                select case(element_id)
                case(carbon12_element)
                   seed_stoich = 1._r8
                case(nitrogen_element)
                   seed_stoich = currentPatch%nitr_repro_stoich(pft)
                case(phosphorus_element)
                   seed_stoich = currentPatch%phos_repro_stoich(pft)
                case default
                   write(fates_log(), *) 'undefined element specified'
                   write(fates_log(), *) 'while defining forced external seed mass flux'
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end select
                
                ! Seed input from external sources (user param seed rain, or dispersal model)
                ! Include both prescribed seed_suppl and seed_in dispersed from neighbouring gridcells
                seed_in_external = seed_stoich*(currentSite%seed_in(pft)/area + EDPftvarcon_inst%seed_suppl(pft)*years_per_day) ![kg/m2/day]
                litt%seed_in_extern(pft) = litt%seed_in_extern(pft) + seed_in_external
                
                ! Seeds entering externally [kg/site/day]
                site_mass%seed_in = site_mass%seed_in + seed_in_external*currentPatch%area
             end if !use this pft  
          enddo

          currentPatch => currentPatch%younger
       enddo seed_in_loop

       ! Determine the total site-level seed output for the current element and update the seed_out mass
       ! for each element loop since the site_seed_rain is resent and updated for each element loop iteration
       do pft = 1,numpft
          site_mass%seed_out = site_mass%seed_out + site_seed_rain(pft)*site_disp_frac(pft) ![kg/site/day]
          currentSite%seed_out(pft) = currentSite%seed_out(pft) + site_seed_rain(pft)*site_disp_frac(pft) ![kg/site/day]
       end do
 
    end do el_loop

    return
  end subroutine SeedUpdate

  ! ============================================================================

  subroutine SeedDecay( litt , currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! 1. Flux from seed pool into leaf litter pool
    ! 2. If the TRS with seedling dynamics is on (regeneration_model = 3)
    !    then we calculate seedling mortality here (i.e. flux from seedling pool
    !    (into leaf litter pool)   
    !
    ! !ARGUMENTS
    type(litter_type) :: litt
    type(fates_patch_type), intent(in) :: currentPatch ! ahb added this
    type(bc_in_type), intent(in) :: bc_in ! ahb added this    
    !
    ! !LOCAL VARIABLES:
    integer  ::  pft
    real(r8) ::  seedling_layer_par          ! cumulative sum of PAR at the seedling layer (MJ)
                                             ! over prior window of days defined by 
                                             ! fates_trs_seedling_mort_par_timescale
    real(r8) ::  seedling_light_mort_rate    ! daily seedling mortality rate from light stress
    real(r8) ::  seedling_h2o_mort_rate      ! daily seedling mortality rate from moisture stress
    real(r8) ::  seedling_mdds               ! moisture deficit days accumulated in the seedling layer
   
    !----------------------------------------------------------------------

    
    ! 1. Seed mortality (i.e. flux from seed bank to litter)
    
    ! default value from Liscke and Loffler 2006 ; making this a PFT-specific parameter
    ! decays the seed pool according to exponential model
    ! seed_decay_rate is in yr-1
    ! seed_decay is kg/day
    ! Assume that decay rates are same for all chemical species

    !=====================================================================================
    do pft = 1,numpft 
    
       ! If the TRS is switched off or the pft can't get big enough to be considered a tree 
       ! then use FATES default regeneration.
       if ( regeneration_model == default_regeneration .or. &
            prt_params%allom_dbh_maxheight(pft) < min_max_dbh_for_trees ) then

          ! Default seed decay (TRS is off)
          litt%seed_decay(pft) = litt%seed(pft) * &
               EDPftvarcon_inst%seed_decay_rate(pft)*years_per_day

       end if

       ! If the TRS is switched on and the pft is a tree then add non-seed reproductive biomass
       ! to the seed decay flux. This was added to litt%seed_decay in the previously called SeedIn 
       ! subroutine
       if ( any(regeneration_model == [TRS_regeneration, TRS_no_seedling_dyn]) .and. &
            prt_params%allom_dbh_maxheight(pft) > min_max_dbh_for_trees ) then
          
          litt%seed_decay(pft) = litt%seed_decay(pft) + &! From non-seed reproductive biomass (added in
               ! in the SeedIn subroutine.
               litt%seed(pft) * EDPftvarcon_inst%seed_decay_rate(pft)*years_per_day
          
       end if 


       ! If the TRS is switched on with seedling dynamics (regeneration_model = 2) 
       ! then calculate seedling mortality.
       if_trs_germ_decay: if ( regeneration_model == TRS_regeneration .and. &
            prt_params%allom_dbh_maxheight(pft) > min_max_dbh_for_trees ) then
          
          !----------------------------------------------------------------------
          ! Seedling mortality (flux from seedling pool to litter)
          ! Note: The TRS uses the litt%seed_germ data struture to track seedlings
          !
          ! Step 1. Calculate the daily seedling mortality rate from light stress
          !
          ! Calculate the cumulative light at the seedling layer over a prior number of 
          ! days determined by the "fates_tres_seedling_mort_par_timescale" parameter.

          seedling_layer_par = currentPatch%sdlng_mort_par%GetMean() * megajoules_per_joule * & 
               sec_per_day * sdlng_mort_par_timescale 
          
          ! Calculate daily seedling mortality rate from light
          seedling_light_mort_rate = exp( EDPftvarcon_inst%seedling_light_mort_a(pft) * &
               seedling_layer_par + EDPftvarcon_inst%seedling_light_mort_b(pft) ) 
        
          ! Step 2. Calculate the daily seedling mortality rate from moisture stress
          
          ! Get the current seedling moisture deficit days (tracked as a pft-specific exponential
          ! average)
          seedling_mdds = currentPatch%sdlng_mdd(pft)%p%GetMean()     
          
          ! Calculate seedling mortality as a function of moisture deficit days (mdd)
          ! If the seedling mmd value is below a critical threshold then moisture-based mortality is zero
          if (seedling_mdds < EDPftvarcon_inst%seedling_mdd_crit(pft)) then
             seedling_h2o_mort_rate = 0.0_r8
          else
             seedling_h2o_mort_rate = EDPftvarcon_inst%seedling_h2o_mort_a(pft) * seedling_mdds**2 + &
                  EDPftvarcon_inst%seedling_h2o_mort_b(pft) * seedling_mdds + &
                  EDPftvarcon_inst%seedling_h2o_mort_c(pft)
          end if ! mdd threshold check
          
          ! Step 3. Sum modes of mortality (including background mortality) and send dead seedlings
          ! to litter        
          litt%seed_germ_decay(pft) = (litt%seed_germ(pft) * seedling_light_mort_rate) + &
               (litt%seed_germ(pft) * seedling_h2o_mort_rate) + &
               (litt%seed_germ(pft) * EDPftvarcon_inst%background_seedling_mort(pft) &
               * years_per_day)
       
       else
          
          litt%seed_germ_decay(pft) = litt%seed_germ(pft) * &
               EDPftvarcon_inst%seed_decay_rate(pft)*years_per_day

       end if if_trs_germ_decay
       
    enddo
    
    return
  end subroutine SeedDecay

  ! ============================================================================
  subroutine SeedGermination( litt, cold_stat, drought_stat, bc_in, currentPatch )
    !
    ! !DESCRIPTION:
    !  Flux from seed bank into the seedling pool    
    !
    ! !USES:

    !
    ! !ARGUMENTS
    type(litter_type) :: litt
    integer                   , intent(in) :: cold_stat    ! Is the site in cold leaf-off status?
    integer, dimension(numpft), intent(in) :: drought_stat ! Is the site in drought leaf-off status?
    type(bc_in_type),           intent(in) :: bc_in
    type(fates_patch_type),        intent(in) :: currentPatch
    !
    ! !LOCAL VARIABLES:
    integer :: pft
    real(r8), parameter ::  max_germination = 1.0_r8 ! Cap on germination rates. 
                                                    ! KgC/m2/yr Lishcke et al. 2009

    !Light and moisture-sensitive seedling emergence variables (ahb)
    !------------------------------------------------------------------------------------------------------------
    integer  :: ilayer_seedling_root           ! the soil layer at seedling rooting depth
    real(r8) :: seedling_layer_smp             ! soil matric potential at seedling rooting depth
    real(r8) :: wetness_index                  ! a soil 'wetness index' (1 / - SoilMatricPotetial (MPa) )
    real(r8) :: seedling_layer_par             ! par at the seedling layer (MJ m-2 day-1)
    real(r8) :: slsmp_emerg                    ! temp
    real(r8) :: slparmort                      ! temp
    real(r8) :: slpartrans                     ! temp
    real(r8) :: photoblastic_germ_modifier     ! seedling emergence rate modifier for light-sensitive germination
    real(r8) :: seedling_emerg_rate            ! the fraction of the seed bank emerging in the current time step
    !-------------------------------------------------------------------------------------------------------------


    ! Turning of this cap? because the cap will impose changes on proportionality
    ! of nutrients. (RGK 02-2019)
    !real(r8), parameter :: max_germination = 1.e6_r8  ! Force to very high number

    !----------------------------------------------------------------------

    ! germination_rate is being pulled to PFT parameter; units are 1/yr
    ! thus the mortality rate of seed -> recruit (in units of carbon)
    ! is seed_decay_rate(p)/germination_rate(p)
    ! and thus the mortality rate (in units of individuals) is the product of
    ! that times the ratio of (hypothetical) seed mass to recruit biomass
    
    !==============================================================================================
    do pft = 1,numpft

       ! If the TRS's seedling dynamics is switched off, then we use FATES's default approach
       ! to germination 
       if_tfs_or_def: if ( regeneration_model == default_regeneration .or. &
            regeneration_model == TRS_no_seedling_dyn .or. & 
            prt_params%allom_dbh_maxheight(pft) < min_max_dbh_for_trees ) then

          litt%seed_germ_in(pft) =  min(litt%seed(pft) * EDPftvarcon_inst%germination_rate(pft), &  
               max_germination)*years_per_day

          ! If TRS seedling dynamics is switched on we calculate seedling emergence (i.e. germination)
          ! as a pft-specific function of understory light and soil moisture.
       else if ( regeneration_model == TRS_regeneration .and. &
            prt_params%allom_dbh_maxheight(pft) > min_max_dbh_for_trees ) then	    

          ! Step 1. Calculate how germination rate is modified by understory light
          ! This applies to photoblastic germinators (e.g. many tropical pioneers) 

          ! Calculate mean PAR at the seedling layer (MJ m-2 day-1) over the prior 24 hours
          seedling_layer_par = currentPatch%seedling_layer_par24%GetMean() * sec_per_day * megajoules_per_joule

          ! Calculate the photoblastic germination rate modifier (Eq. 3 Hanbury-Brown et al., 2022) 
          photoblastic_germ_modifier = seedling_layer_par / &
               (seedling_layer_par + EDPftvarcon_inst%par_crit_germ(pft))

          ! Step 2. Calculate how germination rate is modified by soil moisture in the rooting zone of
          ! the seedlings. This is a pft-specific running mean based on pft-specific seedling rooting
          ! depth.

          ! Get running mean of soil matric potential (mm of H2O suction) at the seedling rooting depth
          ! This running mean based on pft-specific seedling rooting depth.
          seedling_layer_smp = currentPatch%sdlng_emerg_smp(pft)%p%GetMean()    

          ! Calculate a soil wetness index (1 / -soil matric pontential (MPa) ) used by the TRS
          ! to calculate seedling mortality from moisture stress. 
          wetness_index = 1.0_r8 / (seedling_layer_smp * (-1.0_r8) * mpa_per_mm_suction)          

          ! Step 3. Calculate the seedling emergence rate based on soil moisture and germination
          ! rate modifier (Step 1). See Eq. 4 of Hanbury-Brown et al., 2022

          ! If SMP is below a pft-specific value, then no germination occurs
          if ( seedling_layer_smp .GE. EDPftvarcon_inst%seedling_psi_emerg(pft) ) then
             seedling_emerg_rate = photoblastic_germ_modifier * EDPftvarcon_inst%a_emerg(pft) * &
                  wetness_index**EDPftvarcon_inst%b_emerg(pft)
          else 

             seedling_emerg_rate = 0.0_r8

          end if ! End soil-moisture based seedling emergence rate

          ! Step 4. Calculate the amount of carbon germinating out of the seed bank
          litt%seed_germ_in(pft) = litt%seed(pft) * seedling_emerg_rate

       end if if_tfs_or_def
    
      !set the germination only under the growing season...c.xu

      if ((prt_params%season_decid(pft) == itrue ) .and. &
            (any(cold_stat == [phen_cstat_nevercold,phen_cstat_iscold]))) then
          ! no germination for all PFTs when cold
          litt%seed_germ_in(pft) = 0.0_r8
       endif

       ! Drought deciduous, halt germination when status is shedding, even leaves are not
       ! completely abscissed. MLO
       select case (prt_params%stress_decid(pft))
       case (ihard_stress_decid,isemi_stress_decid)
          if (any(drought_stat(pft) == [phen_dstat_timeoff,phen_dstat_moistoff,phen_dstat_pshed])) then
             litt%seed_germ_in(pft) = 0.0_r8
          end if
       end select

    end do

  end subroutine SeedGermination

  ! =====================================================================================

   subroutine recruitment(currentSite, currentPatch, bc_in)
      !
      ! DESCRIPTION:
      ! spawn new cohorts of juveniles of each PFT
      !

      ! ARGUMENTS:
      type(ed_site_type),     intent(inout)          :: currentSite
      type(fates_patch_type), intent(inout), pointer :: currentPatch
      type(bc_in_type),       intent(in)             :: bc_in

      ! LOCAL VARIABLES:
      class(prt_vartypes),      pointer :: prt                ! PARTEH object
      type(litter_type),        pointer :: litt               ! litter object (carbon right now)
      type(site_massbal_type),  pointer :: site_mass          ! for accounting total in-out mass fluxes
      integer                           :: ft                 ! loop counter for PFTs
      integer                           :: leaf_status        ! cohort phenology status [leaves on/off]
      integer                           :: el                 ! loop counter for element
      integer                           :: element_id         ! element index consistent with definitions in PRTGenericMod
      integer                           :: iage               ! age loop counter for leaf age bins
      integer                           :: crowndamage        ! crown damage class of the cohort [1 = undamaged, >1 = damaged]  
      real(r8)                          :: height             ! new cohort height [m]
      real(r8)                          :: dbh                ! new cohort DBH [cm]
      real(r8)                          :: cohort_n           ! new cohort density 
      real(r8)                          :: l2fr               ! leaf to fineroot biomass ratio [0-1]
      real(r8)                          :: c_leaf             ! target leaf biomass [kgC]
      real(r8)                          :: c_fnrt             ! target fine root biomass [kgC]
      real(r8)                          :: c_sapw             ! target sapwood biomass [kgC]
      real(r8)                          :: a_sapw             ! target sapwood cross section are [m2] (dummy)
      real(r8)                          :: c_agw              ! target Above ground biomass [kgC]
      real(r8)                          :: c_bgw              ! target Below ground biomass [kgC]
      real(r8)                          :: c_struct           ! target Structural biomass [kgc]
      real(r8)                          :: c_store            ! target Storage biomass [kgC]
      real(r8)                          :: m_leaf             ! leaf mass (element agnostic) [kg]
      real(r8)                          :: m_fnrt             ! fine-root mass (element agnostic) [kg]
      real(r8)                          :: m_sapw             ! sapwood mass (element agnostic) [kg]
      real(r8)                          :: m_agw              ! AG wood mass (element agnostic) [kg]
      real(r8)                          :: m_bgw              ! BG wood mass (element agnostic) [kg]
      real(r8)                          :: m_struct           ! structural mass (element agnostic) [kg]
      real(r8)                          :: m_store            ! storage mass (element agnostic) [kg]
      real(r8)                          :: m_repro            ! reproductive mass (element agnostic) [kg]
      real(r8)                          :: efleaf_coh         
      real(r8)                          :: effnrt_coh 
      real(r8)                          :: efstem_coh 
      real(r8)                          :: mass_avail         ! mass of each nutrient/carbon available in the seed_germination pool [kg]
      real(r8)                          :: mass_demand        ! total mass demanded by the plant to achieve the stoichiometric 
                                          !    targets of all the organs in the recruits. Used for both [kg per plant] and [kg per cohort] 
      real(r8)                          :: stem_drop_fraction ! 
      real(r8)                          :: fnrt_drop_fraction ! 
      real(r8)                          :: sdlng2sap_par      ! running mean of PAR at the seedling layer [MJ/m2/day]
      real(r8)                          :: seedling_layer_smp ! soil matric potential at seedling rooting depth [mm H2O suction]
      integer, parameter                :: recruitstatus = 1  ! whether the newly created cohorts are recruited or initialized
      integer                           :: ilayer_seedling_root ! the soil layer at seedling rooting depth
      logical                           :: use_this_pft         ! logical flag for whether or not to allow a given PFT to recruit
      !---------------------------------------------------------------------------

      do ft = 1, numpft

       ! The following if block is for the prescribed biogeography and/or nocomp modes and/or crop land use types
       ! Since currentSite%use_this_pft is a site-level quantity and thus only limits whether a given PFT
       ! is permitted on a given gridcell or not, it applies to the prescribed biogeography case only.
       ! If nocomp is enabled, then we must determine whether a given PFT is allowed on a given patch or not.
       ! Whether or not nocomp or prescribed biogeography is enabled, if land use change is enabled, then we only want to
       ! allow crop PFTs on patches with crop land use types

       use_this_pft = .false.
       if(currentSite%use_this_pft(ft).eq.itrue &
            .and. ((hlm_use_nocomp .eq. ifalse) .or. (ft .eq. currentPatch%nocomp_pft_label)))then
          use_this_pft = .true.
       end if

       if ( currentPatch%land_use_label .ne. nocomp_bareground_land ) then ! cdk
          if ((hlm_use_luh .eq. itrue) .and. (is_crop(currentPatch%land_use_label))) then
             if ( crop_lu_pft_vector(currentPatch%land_use_label) .eq. ft ) then
                use_this_pft = .true.
             else
                use_this_pft = .false.
             end if
          end if
       endif

       use_this_pft_if: if(use_this_pft) then
            height             = EDPftvarcon_inst%hgt_min(ft)
            stem_drop_fraction = prt_params%phen_stem_drop_fraction(ft)
            fnrt_drop_fraction = prt_params%phen_fnrt_drop_fraction(ft)
            l2fr               = currentSite%rec_l2fr(ft, currentPatch%NCL_p)
            crowndamage        = 1 ! new recruits are undamaged

            ! calculate DBH from initial height 
            call h2d_allom(height, ft, dbh)

            ! default assumption is that leaves are on
            efleaf_coh  = 1.0_r8
            effnrt_coh  = 1.0_r8
            efstem_coh  = 1.0_r8
            leaf_status = leaves_on

            ! but if the plant is seasonally (cold) deciduous, and the site status is flagged
            ! as "cold", then set the cohort's status to leaves_off, and remember the leaf biomass
            if ((prt_params%season_decid(ft) == itrue) .and.                   &
               (any(currentSite%cstatus == [phen_cstat_nevercold, phen_cstat_iscold]))) then
               efleaf_coh  = 0.0_r8
               effnrt_coh  = 1.0_r8 - fnrt_drop_fraction
               efstem_coh  = 1.0_r8 - stem_drop_fraction
               leaf_status = leaves_off
            end if 

            ! Or.. if the plant is drought deciduous, make sure leaf status is consistent with the
            ! leaf elongation factor.
            ! For tissues other than leaves, the actual drop fraction is a combination of the
            ! elongation factor (e) and the drop fraction (x), which will ensure that the remaining
            ! tissue biomass will be exactly e when x=1, and exactly the original biomass when x = 0.
            select case (prt_params%stress_decid(ft))
            case (ihard_stress_decid, isemi_stress_decid)
               efleaf_coh = currentSite%elong_factor(ft)
               effnrt_coh = 1.0_r8 - (1.0_r8 - efleaf_coh)*fnrt_drop_fraction
               efstem_coh = 1.0_r8 - (1.0_r8 - efleaf_coh)*stem_drop_fraction

               ! For the initial state, we always assume that leaves are flushing (instead of partially abscissing)
               ! whenever the elongation factor is non-zero.  If the elongation factor is zero, then leaves are in
               ! the "off" state.
               if (efleaf_coh > 0.0_r8) then
                  leaf_status = leaves_on 
               else 
                  leaf_status = leaves_off
               end if
            end select

            ! calculate live pools
            call bleaf(dbh, ft, crowndamage, init_recruit_trim, efleaf_coh,    &
               c_leaf)
            call bfineroot(dbh, ft, init_recruit_trim, l2fr, effnrt_coh, c_fnrt)
            call bsap_allom(dbh, ft, crowndamage, init_recruit_trim,           &
               efstem_coh, a_sapw, c_sapw)
            call bagw_allom(dbh, ft, crowndamage, efstem_coh, c_agw)
            call bbgw_allom(dbh, ft, efstem_coh, c_bgw)
            call bdead_allom(c_agw, c_bgw, c_sapw, ft, c_struct)
            call bstore_allom(dbh, ft, crowndamage, init_recruit_trim, c_store)

            ! cycle through available carbon and nutrients, find the limiting element
            ! to dictate the total number of plants that can be generated
            if_not_prescribed: if ((hlm_use_ed_prescribed_phys .eq. ifalse) .or.                  &
               (EDPftvarcon_inst%prescribed_recruitment(ft) .lt. 0._r8)) then

               cohort_n = 1.e20_r8

               do_elem: do el = 1, num_elements
                  element_id = element_list(el)
                  select case(element_id)
                  case(carbon12_element)
                     mass_demand = c_struct + c_leaf + c_fnrt + c_sapw + c_store
                  case(nitrogen_element)
                     mass_demand =                                                                     &
                     c_struct*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(struct_organ)) + &
                     c_leaf*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(leaf_organ))     + &
                     c_fnrt*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ))     + &
                     c_sapw*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(sapw_organ))     + &
                     StorageNutrientTarget(ft, element_id,                                             &
                     c_leaf*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(leaf_organ)),      &
                     c_fnrt*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ)),      &
                     c_sapw*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(sapw_organ)),      &
                     c_struct*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(struct_organ)))
                  case(phosphorus_element)
                     mass_demand =                                                                     &
                     c_struct*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(struct_organ)) + &
                     c_leaf*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(leaf_organ)) +     &
                     c_fnrt*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ)) +     &
                     c_sapw*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(sapw_organ)) +     &
                     StorageNutrientTarget(ft, element_id,                                             &
                     c_leaf*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(leaf_organ)),      &
                     c_fnrt*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ)),      &
                     c_sapw*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(sapw_organ)),      &
                     c_struct*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(struct_organ)))
                  case default
                     write(fates_log(),*) 'Undefined element type in recruitment'
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  end select

                  ! If TRS seedling dynamics is switched off then the available mass to make new recruits
                  ! is everything in the seed_germ pool.
                  if (regeneration_model == default_regeneration .or.          &
                     regeneration_model == TRS_no_seedling_dyn .or.            & 
                     prt_params%allom_dbh_maxheight(ft) < min_max_dbh_for_trees) then

                     mass_avail = currentPatch%area * currentPatch%litter(el)%seed_germ(ft)

                     ! If TRS seedling dynamics is on then calculate the available mass to make new recruits
                     ! as a pft-specific function of light and soil moisture in the seedling layer.
                  else if (regeneration_model == TRS_regeneration .and.        &
                     prt_params%allom_dbh_maxheight(ft) > min_max_dbh_for_trees) then

                     sdlng2sap_par = currentPatch%sdlng2sap_par%GetMean()*     &
                        sec_per_day*megajoules_per_joule

                     mass_avail = currentPatch%area*                           &
                        currentPatch%litter(el)%seed_germ(ft)*                 & 
                        EDPftvarcon_inst%seedling_light_rec_a(ft)*             &
                        sdlng2sap_par**EDPftvarcon_inst%seedling_light_rec_b(ft) 

                     ! If soil moisture is below pft-specific seedling  moisture stress threshold the 
                     ! recruitment does not occur.
                     ilayer_seedling_root = minloc(abs(bc_in%z_sisl(:) -       &
                        EDPftvarcon_inst%seedling_root_depth(ft)), dim=1)

                     seedling_layer_smp = bc_in%smp_sl(ilayer_seedling_root)

                     if (seedling_layer_smp < EDPftvarcon_inst%seedling_psi_crit(ft)) then
                        mass_avail = 0.0_r8
                     end if 

                  end if ! End use TRS with seedling dynamics

                  ! update number density if this is the limiting mass
                  cohort_n = min(cohort_n, mass_avail/mass_demand)

               end do do_elem

            else
               ! prescribed recruitment rates. number per sq. meter per year
               cohort_n = currentPatch%area*EDPftvarcon_inst%prescribed_recruitment(ft) *  &
                  hlm_freq_day
            endif if_not_prescribed

            ! Only bother allocating a new cohort if there is a reasonable amount of it
            any_recruits: if (cohort_n > min_n_safemath) then

               ! --------------------------------------------------------------------------------
               ! PART II.
               ! Initialize the PARTEH object, and determine the initial masses of all
               ! organs and elements.
               ! --------------------------------------------------------------------------------

               prt => null()
               call InitPRTObject(prt)

               do el = 1,num_elements

                  element_id = element_list(el)

                  ! If this is carbon12, then the initialization is straight forward
                  ! otherwise, we use stoichiometric ratios
                  select case(element_id)
                  case(carbon12_element)
                     m_struct = c_struct
                     m_leaf   = c_leaf
                     m_fnrt   = c_fnrt
                     m_sapw   = c_sapw
                     m_store  = c_store
                     m_repro  = 0._r8
                  case(nitrogen_element)
                     m_struct = c_struct*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(struct_organ))
                     m_leaf   = c_leaf*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(leaf_organ))
                     m_fnrt   = c_fnrt*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ))
                     m_sapw   = c_sapw*prt_params%nitr_stoich_p1(ft, prt_params%organ_param_id(sapw_organ))
                     m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
                     m_repro  = 0._r8
                  case(phosphorus_element)
                     m_struct = c_struct*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(struct_organ))
                     m_leaf   = c_leaf*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(leaf_organ))
                     m_fnrt   = c_fnrt*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(fnrt_organ))
                     m_sapw   = c_sapw*prt_params%phos_stoich_p1(ft, prt_params%organ_param_id(sapw_organ))
                     m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct)
                     m_repro  = 0._r8
                  end select

                  
                  select case(hlm_parteh_mode)
                  case (prt_carbon_allom_hyp, prt_cnp_flex_allom_hyp)

                     ! put all of the leaf mass into the first bin
                     call SetState(prt, leaf_organ, element_id, m_leaf, 1)
                     do iage = 2, nleafage
                        call SetState(prt,leaf_organ, element_id, 0._r8, iage)
                     end do

                     call SetState(prt, fnrt_organ, element_id, m_fnrt)
                     call SetState(prt, sapw_organ, element_id, m_sapw)
                     call SetState(prt, store_organ, element_id, m_store)
                     call SetState(prt, struct_organ, element_id, m_struct)
                     call SetState(prt, repro_organ, element_id, m_repro)

                  case default
                     write(fates_log(),*) 'Unspecified PARTEH module during create_cohort'
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  end select

                  site_mass => currentSite%mass_balance(el)

                  ! Remove mass from the germination pool. However, if we are use prescribed physiology,
                  ! AND the forced recruitment model, then we are not realling using the prognostic
                  ! seed_germination model, so we have to short circuit things.  We send all of the
                  ! seed germination mass to an outflux pool, and use an arbitrary generic input flux
                  ! to balance out the new recruits.
                  if ((hlm_use_ed_prescribed_phys .eq. itrue) .and.            &
                     (EDPftvarcon_inst%prescribed_recruitment(ft) .ge. 0._r8)) then

                     site_mass%flux_generic_in = site_mass%flux_generic_in +   &
                     cohort_n*(m_struct + m_leaf + m_fnrt + m_sapw + m_store + m_repro)

                     site_mass%flux_generic_out = site_mass%flux_generic_out + &
                     currentPatch%area * currentPatch%litter(el)%seed_germ(ft)

                     currentPatch%litter(el)%seed_germ(ft) = 0._r8
                  else
                     currentPatch%litter(el)%seed_germ(ft) =                   &
                     currentPatch%litter(el)%seed_germ(ft) - cohort_n / currentPatch%area *   &
                     (m_struct + m_leaf + m_fnrt + m_sapw + m_store + m_repro)
                  end if
                  
               end do

               ! cycle through the initial conditions, and makes sure that they are all initialized
               call prt%CheckInitialConditions()

               call create_cohort(currentSite, currentPatch, ft, cohort_n,     &
                  height, 0.0_r8, dbh, prt, efleaf_coh, effnrt_coh, efstem_coh,  &
                  leaf_status, recruitstatus, init_recruit_trim, 0.0_r8,       &
                  currentPatch%NCL_p, crowndamage, currentSite%spread, bc_in)

               ! Note that if hydraulics is on, the number of cohorts may have
               ! changed due to hydraulic constraints.
               ! This constaint is applied during "create_cohort" subroutine.

               ! keep track of how many individuals were recruited for passing to history
               currentSite%recruitment_rate(ft) = currentSite%recruitment_rate(ft) + cohort_n

            endif any_recruits
         endif use_this_pft_if
      enddo  !pft loop
      call currentPatch%ValidateCohorts()
   end subroutine recruitment

   ! ======================================================================================

  subroutine CWDInput( currentSite, currentPatch, litt, bc_in)

    !
    ! !DESCRIPTION:
    ! Generate litter fields from turnover.
    ! Note, that the when this is called, the number density of the plants
    ! has not been reduced from non-mortal turnover yet.
    ! Thus, we need to avoid double counting losses from dying trees
    ! and turnover in dying trees.
    !
    ! !USES:
    use EDParamsMod           , only : landuse_grazing_carbon_use_eff
    use EDParamsMod           , only : landuse_grazing_nitrogen_use_eff
    use EDParamsMod           , only : landuse_grazing_phosphorus_use_eff
    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout), target :: currentSite
    type(fates_patch_type),intent(inout), target :: currentPatch
    type(litter_type),intent(inout),target    :: litt
    type(bc_in_type),intent(in)               :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(fates_cohort_type), pointer      :: currentCohort
    type(elem_diag_type), pointer :: elflux_diags
    type(site_massbal_type), pointer   :: site_mass
    integer  :: c
    real(r8) :: dead_n          ! total understorey dead tree density
    real(r8) :: dead_n_dlogging ! direct logging understory dead-tree density
    real(r8) :: dead_n_ilogging ! indirect understory dead-tree density (logging)
    real(r8) :: dead_n_natural  ! understory dead density not associated
    ! with direct logging
    real(r8) :: leaf_m          ! mass of the element of interest in the
    ! leaf  [kg]
    real(r8) :: fnrt_m           ! fine-root [kg]
    real(r8) :: sapw_m    ! sapwood [kg]
    real(r8) :: struct_m    ! structural [kg]
    real(r8) :: store_m    ! storage [kg]
    real(r8) :: repro_m    ! reproductive [kg]
    real(r8) :: leaf_m_turnover ! leaf turnover [kg]
    real(r8) :: fnrt_m_turnover
    real(r8) :: sapw_m_turnover
    real(r8) :: struct_m_turnover
    real(r8) :: store_m_turnover
    real(r8) :: repro_m_turnover
    real(r8) :: dcmpy_frac        ! Fraction of mass sent to decomposability pool
    real(r8) :: plant_dens        ! Number of plants per m2
    real(r8) :: bg_cwd_tot        ! Total below-ground coarse woody debris
    ! input flux
    real(r8) :: root_fines_tot    ! Total below-ground fine root coarse
    ! woody debris
    integer  :: element_id        ! element id consistent with parteh/PRTGenericMod.F90

    real(r8) :: trunk_wood        ! carbon flux into trunk products kgC/day/site
    integer  :: ilyr
    integer  :: pft
    integer  :: dcmpy             ! decomposability pool index
    integer  :: numlevsoil        ! Actual number of soil layers

    real(r8) :: SF_val_CWD_frac_adj(4) !SF_val_CWD_frac adjusted based on cohort dbh
    real(r8) :: leaf_herbivory  ! leaf that is eaten by grazers [kg]
    real(r8) :: herbivory_element_use_efficiency   ! fraction of grazed biomass that is returned to litter pool versus atmosphere
    !----------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Other direct litter fluxes happen in phenology and in spawn_patches.
    ! -----------------------------------------------------------------------------------

    numlevsoil = currentSite%nlevsoil

    element_id = litt%element_id

    select case(element_id)
       case (carbon12_element)
          herbivory_element_use_efficiency = landuse_grazing_carbon_use_eff
       case (nitrogen_element)
          herbivory_element_use_efficiency = landuse_grazing_nitrogen_use_eff
       case (phosphorus_element)
          herbivory_element_use_efficiency = landuse_grazing_phosphorus_use_eff
    end select

    ! Object tracking flux diagnostics for each element
    elflux_diags => currentSite%flux_diags%elem(element_pos(element_id))

    ! Object tracking site level mass balance for each element
    site_mass => currentSite%mass_balance(element_pos(element_id))

    ! Transfer litter from turnover of living plants
    
    currentCohort => currentPatch%shortest
    do while(associated(currentCohort))

       pft = currentCohort%pft
       call set_root_fraction(currentSite%rootfrac_scr, pft, currentSite%zi_soil, &
           bc_in%max_rooting_depth_index_col)

       store_m_turnover  = currentCohort%prt%GetTurnover(store_organ,element_id)
       fnrt_m_turnover   = currentCohort%prt%GetTurnover(fnrt_organ,element_id)
       repro_m_turnover  = currentCohort%prt%GetTurnover(repro_organ,element_id)
       

       
       store_m         = currentCohort%prt%GetState(store_organ,element_id)
       fnrt_m          = currentCohort%prt%GetState(fnrt_organ,element_id)
       repro_m         = currentCohort%prt%GetState(repro_organ,element_id)

       leaf_herbivory  = currentCohort%prt%GetHerbivory(leaf_organ, element_id)

       if (prt_params%woody(currentCohort%pft) == itrue) then
          ! Assumption: for woody plants fluxes from deadwood and sapwood go together in CWD pool
          leaf_m_turnover   = currentCohort%prt%GetTurnover(leaf_organ,element_id)
          sapw_m_turnover   = currentCohort%prt%GetTurnover(sapw_organ,element_id)
          struct_m_turnover = currentCohort%prt%GetTurnover(struct_organ,element_id)

          leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id)
          sapw_m          = currentCohort%prt%GetState(sapw_organ,element_id)
          struct_m        = currentCohort%prt%GetState(struct_organ,element_id)
       else
          ! for non-woody plants all stem fluxes go into the same leaf litter pool
          leaf_m_turnover   = currentCohort%prt%GetTurnover(leaf_organ,element_id) + &
               currentCohort%prt%GetTurnover(sapw_organ,element_id) + &
               currentCohort%prt%GetTurnover(struct_organ,element_id)
          sapw_m_turnover   = 0._r8
          struct_m_turnover = 0._r8

          leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id) + &
               currentCohort%prt%GetState(sapw_organ,element_id) + &
               currentCohort%prt%GetState(struct_organ,element_id)
          sapw_m          = 0._r8
          struct_m        = 0._r8
       end if

       plant_dens =  currentCohort%n/currentPatch%area

       ! ---------------------------------------------------------------------------------
       ! PART 1 Litter fluxes from non-mortal tissue turnovers  Kg/m2/day
       !        Important note:  Turnover has already been removed from the cohorts.
       !        So, in the next part of this algorithm, when we send the biomass
       !        from dying trees to the litter pools, we don't have to worry
       !        about double counting.
       ! ---------------------------------------------------------------------------------

       elflux_diags%surf_fine_litter_input(pft) = &
            elflux_diags%surf_fine_litter_input(pft) +  &
            (leaf_m_turnover+repro_m_turnover) * currentCohort%n

       root_fines_tot = (fnrt_m_turnover + store_m_turnover ) * &
            plant_dens

       do dcmpy=1,ndcmpy
          dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)
          litt%leaf_fines_in(dcmpy) = litt%leaf_fines_in(dcmpy) + &
               (leaf_m_turnover+repro_m_turnover + &
               leaf_herbivory * herbivory_element_use_efficiency) * &
               plant_dens * dcmpy_frac

          dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
          do ilyr = 1, numlevsoil
             litt%root_fines_in(dcmpy,ilyr) = litt%root_fines_in(dcmpy,ilyr) + &
                  currentSite%rootfrac_scr(ilyr) * root_fines_tot * dcmpy_frac
          end do
       end do

       elflux_diags%root_litter_input(pft) = &
            elflux_diags%root_litter_input(pft) +  &
            (fnrt_m_turnover + store_m_turnover ) * currentCohort%n

       ! send the part of the herbivory flux that doesn't go to litter to the atmosphere

       site_mass%herbivory_flux_out = &
            site_mass%herbivory_flux_out + &
            leaf_herbivory * (1._r8 - herbivory_element_use_efficiency) * currentCohort%n

       ! Assumption: turnover from deadwood and sapwood are lumped together in CWD pool

       !update partitioning of stem wood (struct + sapw) to cwd based on cohort dbh
       call adjust_SF_CWD_frac(currentCohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)


       do c = 1,ncwd
          litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + &
               (sapw_m_turnover + struct_m_turnover) * &
               SF_val_CWD_frac_adj(c) * plant_dens * &
               prt_params%allom_agb_frac(pft)

          elflux_diags%cwd_ag_input(c)  = elflux_diags%cwd_ag_input(c) + &
               (struct_m_turnover + sapw_m_turnover) * SF_val_CWD_frac_adj(c) * &
               prt_params%allom_agb_frac(pft) * currentCohort%n

          bg_cwd_tot = (sapw_m_turnover + struct_m_turnover) * &
               SF_val_CWD_frac_adj(c) * plant_dens * &
               (1.0_r8-prt_params%allom_agb_frac(pft))

          do ilyr = 1, numlevsoil
             litt%bg_cwd_in(c,ilyr) = litt%bg_cwd_in(c,ilyr) + &
                  bg_cwd_tot * currentSite%rootfrac_scr(ilyr)
          end do

          elflux_diags%cwd_bg_input(c)  = elflux_diags%cwd_bg_input(c) + &
               bg_cwd_tot*currentPatch%area

       enddo


       ! ---------------------------------------------------------------------------------
       ! PART 2 Litter fluxes from non-disturbance inducing  mortality. Kg/m2/day
       ! ---------------------------------------------------------------------------------

       ! Total number of dead (n/m2/day)
       dead_n = -1.0_r8 * currentCohort%dndt/currentPatch%area*years_per_day

       if(currentCohort%canopy_layer > 1)then

          ! Total number of dead understory from direct logging
          ! (it is possible that large harvestable trees are in the understory)
          dead_n_dlogging = currentCohort%lmort_direct * &
               currentCohort%n/currentPatch%area

          ! Total number of dead understory from indirect logging
          dead_n_ilogging = (currentCohort%lmort_collateral + currentCohort%lmort_infra) * &
               currentCohort%n/currentPatch%area

       else

          ! All mortality from logging in the canopy is
          ! is disturbance generating

          dead_n_dlogging = 0._r8
          dead_n_ilogging = 0._r8

       end if

       dead_n_natural = dead_n - dead_n_dlogging - dead_n_ilogging

       ! %n has not been updated due to mortality yet, thus
       ! the litter flux has already been counted since it captured
       ! the losses of live trees and those flagged for death

       root_fines_tot =  dead_n * (fnrt_m + &
            store_m*(1._r8-EDPftvarcon_inst%allom_frbstor_repro(pft)) )

       do dcmpy=1,ndcmpy

          dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)
          litt%leaf_fines_in(dcmpy) = litt%leaf_fines_in(dcmpy) + &
               (leaf_m+repro_m) * dead_n * dcmpy_frac

          dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
          do ilyr = 1, numlevsoil
             litt%root_fines_in(dcmpy,ilyr) = litt%root_fines_in(dcmpy,ilyr) + &
                  root_fines_tot * currentSite%rootfrac_scr(ilyr) * dcmpy_frac
          end do
       end do

       elflux_diags%surf_fine_litter_input(pft) = &
            elflux_diags%surf_fine_litter_input(pft) +  &
            (leaf_m+repro_m) * dead_n*currentPatch%area

       elflux_diags%root_litter_input(pft) = &
            elflux_diags%root_litter_input(pft) +  &
            root_fines_tot*currentPatch%area

       ! Track CWD inputs from dead plants

       do c = 1,ncwd

          ! Below-ground

          bg_cwd_tot = (struct_m + sapw_m) * &
               SF_val_CWD_frac_adj(c) * dead_n * &
               (1.0_r8-prt_params%allom_agb_frac(pft))

          do ilyr = 1, numlevsoil
             litt%bg_cwd_in(c,ilyr) = litt%bg_cwd_in(c,ilyr) + &
                  currentSite%rootfrac_scr(ilyr) * bg_cwd_tot
          end do

          elflux_diags%cwd_bg_input(c)  = elflux_diags%cwd_bg_input(c) + &
               bg_cwd_tot * currentPatch%area

          ! Send AGB component of boles from logging activities into the litter.
          ! This includes fluxes from indirect modes of death, as well as the
          ! non-exported boles due to direct harvesting.

          if (c==ncwd) then


             trunk_wood =  (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * dead_n_dlogging * &
                  prt_params%allom_agb_frac(pft)

             site_mass%wood_product_harvest(pft) = site_mass%wood_product_harvest(pft) + &
                  trunk_wood * currentPatch%area * logging_export_frac

             ! Add AG wood to litter from the non-exported fraction of wood
             ! from direct anthro sources

             litt%ag_cwd_in(c) = litt%ag_cwd_in(c) +  &
                  trunk_wood * (1._r8-logging_export_frac)

             elflux_diags%cwd_ag_input(c)  = elflux_diags%cwd_ag_input(c) + &
                  trunk_wood * (1._r8-logging_export_frac) * currentPatch%area

             ! Add AG wood to litter from indirect anthro sources

             litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * (dead_n_natural+dead_n_ilogging)  * &
                  prt_params%allom_agb_frac(pft)

             elflux_diags%cwd_ag_input(c)  = elflux_diags%cwd_ag_input(c) + &
                  (struct_m + sapw_m) * SF_val_CWD_frac_adj(c) * (dead_n_natural+dead_n_ilogging) * &
                  currentPatch%area * prt_params%allom_agb_frac(pft)

          else

             litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * dead_n  * &
                  prt_params%allom_agb_frac(pft)

             elflux_diags%cwd_ag_input(c)  = elflux_diags%cwd_ag_input(c) + &
                  SF_val_CWD_frac_adj(c) * dead_n * (struct_m + sapw_m) * &
                  currentPatch%area * prt_params%allom_agb_frac(pft)

          end if

       end do


       ! Update diagnostics that track resource management

       if( element_id .eq. carbon12_element ) then

          currentSite%resources_management%delta_litter_stock  = &
               currentSite%resources_management%delta_litter_stock + &
               (leaf_m + fnrt_m + store_m ) * &
               (dead_n_ilogging+dead_n_dlogging) * currentPatch%area

          currentSite%resources_management%delta_biomass_stock = &
               currentSite%resources_management%delta_biomass_stock + &
               (leaf_m + fnrt_m + store_m ) * &
               (dead_n_ilogging+dead_n_dlogging) *currentPatch%area

          currentSite%resources_management%trunk_product_site = &
               currentSite%resources_management%trunk_product_site + &
               trunk_wood * logging_export_frac * currentPatch%area

          do c = 1,ncwd
             currentSite%resources_management%delta_litter_stock  = &
                  currentSite%resources_management%delta_litter_stock + &
                  (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * (dead_n_natural+dead_n_ilogging) * &
                  currentPatch%area

             currentSite%resources_management%delta_biomass_stock = &
                  currentSite%resources_management%delta_biomass_stock + &
                  (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * dead_n * currentPatch%area
          end do

          ! Update diagnostics that track resource management
          currentSite%resources_management%delta_individual    = &
               currentSite%resources_management%delta_individual + &
               (dead_n_dlogging+dead_n_ilogging) * hlm_freq_day * currentPatch%area
       end if


       currentCohort => currentCohort%taller
    enddo  ! end loop over cohorts


    return
  end subroutine CWDInput

  ! =====================================================================================


  subroutine fragmentation_scaler( currentPatch, bc_in)
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! FIX(SPM, 091914) this should be a function as it returns a value in
    ! currentPatch%fragmentation_scaler
    !
    ! !USES:
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : pi => pi_const
    !
    ! !ARGUMENTS
    type(fates_patch_type), intent(inout) :: currentPatch
    type(bc_in_type),    intent(in)    :: bc_in

    !
    ! !LOCAL VARIABLES:
    logical  :: use_century_tfunc = .false.
    logical  :: use_hlm_soil_scalar = .true. ! Use hlm input decomp fraction scalars
    integer  :: j
    real(r8) :: t_scalar                     ! temperature scalar
    real(r8) :: w_scalar                     ! moisture scalar
    real(r8) :: catanf                       ! hyperbolic temperature function from CENTURY
    real(r8) :: catanf_30                    ! hyperbolic temperature function from CENTURY
    real(r8) :: t1                           ! temperature argument
    !----------------------------------------------------------------------

    catanf(t1) = 11.75_r8 +(29.7_r8 / pi) * atan( pi * 0.031_r8  * ( t1 - 15.4_r8 ))
    catanf_30 = catanf(30._r8)

    if(currentPatch%nocomp_pft_label.ne.nocomp_bareground)then

       ! Use the hlm temp and moisture decomp fractions by default
       if ( use_hlm_soil_scalar ) then

         ! Calculate the fragmentation_scaler
         currentPatch%fragmentation_scaler =  min(1.0_r8,max(0.0_r8,bc_in%t_scalar_sisl * bc_in%w_scalar_sisl))

       else

         if ( .not. use_century_tfunc ) then
            !calculate rate constant scalar for soil temperature,assuming that the base rate constants
            !are assigned for non-moisture limiting conditions at 25C.
            if (currentPatch%tveg24%GetMean()  >=  tfrz) then
               t_scalar = q10_mr**((currentPatch%tveg24%GetMean()-(tfrz+25._r8))/10._r8)
               !  Q10**((t_soisno(c,j)-(tfrz+25._r8))/10._r8)
            else
               t_scalar = (q10_mr**(-25._r8/10._r8))*(q10_froz**((currentPatch%tveg24%GetMean()-tfrz)/10._r8))
               !  Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-tfrz)/10._r8)
            endif
         else
            ! original century uses an arctangent function to calculate the
            ! temperature dependence of decomposition
            t_scalar = max(catanf(currentPatch%tveg24%GetMean()-tfrz)/catanf_30,0.01_r8)
         endif

         !Moisture Limitations
         !BTRAN APPROACH - is quite simple, but max's out decomp at all unstressed
         !soil moisture values, which is not realistic.
         !litter decomp is proportional to water limitation on average...
         w_scalar = sum(currentPatch%btran_ft(1:numpft))/real(numpft,r8)

       ! Calculate the fragmentation_scaler
         currentPatch%fragmentation_scaler(:) =  min(1.0_r8,max(0.0_r8,t_scalar * w_scalar))

      endif ! scalar

    endif ! not bare ground

  end subroutine fragmentation_scaler

  ! ============================================================================

  subroutine CWDOut( litt, fragmentation_scaler, nlev_eff_decomp )
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! spawn new cohorts of juveniles of each PFT
    !
    ! !USES:
    use SFParamsMod, only : SF_val_max_decomp

    !
    ! !ARGUMENTS
    type(litter_type),intent(inout),target     :: litt
    real(r8),intent(in)                        :: fragmentation_scaler(:)

    ! This is not necessarily every soil layer, this is the number
    ! of effective layers that are active and can be sent
    ! to the soil decomposition model
    integer,intent(in)                         :: nlev_eff_decomp

    !
    ! !LOCAL VARIABLES:
    integer :: c                       ! Fuel size class index
    integer :: ilyr                    ! Soil layer index
    integer :: dcmpy                   ! Decomposibility pool indexer
    integer :: soil_layer_index = 1    ! Soil layer index associated with above ground litter
    !----------------------------------------------------------------------


    ! Above ground litters are associated with the top soil layer temperature and
    ! moisture scalars and fragmentation scalar associated with specified index value
    ! is used for ag_cwd_frag and root_fines_frag calculations.

    do c = 1,ncwd

       litt%ag_cwd_frag(c)   = litt%ag_cwd(c) * SF_val_max_decomp(c) * &
             years_per_day * fragmentation_scaler(soil_layer_index)

       do ilyr = 1,nlev_eff_decomp
           litt%bg_cwd_frag(c,ilyr) = litt%bg_cwd(c,ilyr) * SF_val_max_decomp(c) * &
                years_per_day * fragmentation_scaler(ilyr)
       enddo
    end do

    ! this is the rate at which dropped leaves stop being part of the burnable pool
    ! and begin to be part of the decomposing pool. This should probably be highly
    ! sensitive to moisture, but also to the type of leaf thick leaves can dry out
    ! before they are decomposed, for example. This section needs further scientific input.

    do dcmpy = 1,ndcmpy

       litt%leaf_fines_frag(dcmpy) = litt%leaf_fines(dcmpy) * &
             years_per_day * SF_val_max_decomp(fuel_classes%dead_leaves()) * fragmentation_scaler(soil_layer_index)

       do ilyr = 1,nlev_eff_decomp
           litt%root_fines_frag(dcmpy,ilyr) = litt%root_fines(dcmpy,ilyr) * &
                 years_per_day *  SF_val_max_decomp(fuel_classes%dead_leaves()) * fragmentation_scaler(ilyr)
       end do
    enddo

  end subroutine CWDOut
  
  subroutine UpdateRecruitL2FR(csite)
    

    ! When CNP is active, the l2fr (target leaf to fine-root biomass multiplier)
    ! is dynamic. We therefore update what the l2fr for recruits
    ! are, taking an exponential moving average of all plants that
    ! are within recruit size limitations (less than recruit size + delta)
    ! and less than the max_count cohort.
    
    type(ed_site_type) :: csite
    type(fates_patch_type), pointer :: cpatch
    type(fates_cohort_type), pointer :: ccohort

    real(r8) :: rec_n(maxpft,nclmax)     ! plant count
    real(r8) :: rec_l2fr0(maxpft,nclmax) ! mean l2fr for this day
    integer  :: rec_count(maxpft,nclmax) ! sample count
    integer  :: ft                       ! functional type index
    integer  :: cl                       ! canopy layer index
    real(r8) :: dbh_min                  ! the dbh of a recruit
    real(r8), parameter :: max_delta = 5.0_r8  ! dbh tolerance, cm, consituting a recruit
    real(r8), parameter :: smth_wgt = 1._r8/300.0_r8
    integer, parameter :: max_count = 3
    
    ! Difference in dbh (cm) to consider a plant was recruited fairly recently

    if(hlm_parteh_mode .ne. prt_cnp_flex_allom_hyp) return
    
    rec_n(1:numpft,1:nclmax) = 0._r8
    rec_l2fr0(1:numpft,1:nclmax) = 0._r8

    cpatch => csite%youngest_patch
    do while(associated(cpatch))

       rec_count(1:numpft,1:nclmax) = 0
       
       ccohort => cpatch%shortest
       cloop: do while(associated(ccohort))

          ft = ccohort%pft
          cl = ccohort%canopy_layer
          call h2d_allom(EDPftvarcon_inst%hgt_min(ft),ft,dbh_min)

          if( .not.ccohort%isnew ) then

             if(rec_count(ft,cl) <= max_count .and. &
                  ccohort%dbh-dbh_min < max_delta ) then
                rec_count(ft,cl) = rec_count(ft,cl) + 1
                rec_n(ft,cl) = rec_n(ft,cl) + ccohort%n
                rec_l2fr0(ft,cl) = rec_l2fr0(ft,cl) + ccohort%n*ccohort%l2fr
             end if

          end if

          ccohort => ccohort%taller
       end do cloop

       cpatch => cpatch%older
    end do

    ! Find the daily mean for each PFT weighted by number and add it to the running mean
    do cl = 1,nclmax
       do ft = 1,numpft
          if(rec_n(ft,cl)>nearzero)then
             rec_l2fr0(ft,cl) = rec_l2fr0(ft,cl) / rec_n(ft,cl)
             csite%rec_l2fr(ft,cl) = &
                  (1._r8-smth_wgt)*csite%rec_l2fr(ft,cl) + smth_wgt*rec_l2fr0(ft,cl)
          end if
       end do
    end do

    return
  end subroutine UpdateRecruitL2FR

  ! ======================================================================

  subroutine UpdateRecruitStoich(csite)

    type(ed_site_type) :: csite
    type(fates_patch_type), pointer :: cpatch
    type(fates_cohort_type), pointer :: ccohort
    integer  :: ft                       ! functional type index
    integer  :: cl                       ! canopy layer index
    real(r8) :: rec_l2fr_pft             ! Actual l2fr of a pft in it's patch
    
    ! Update the total plant stoichiometry of a new recruit, based on the updated
    ! L2FR values

    if(hlm_parteh_mode .ne. prt_cnp_flex_allom_hyp) return
    
    cpatch => csite%youngest_patch
    do while(associated(cpatch))
       cl = cpatch%ncl_p
       
       do ft = 1,numpft
          rec_l2fr_pft = csite%rec_l2fr(ft,cl)
          cpatch%nitr_repro_stoich(ft) = &
               NewRecruitTotalStoichiometry(ft,rec_l2fr_pft,nitrogen_element)
          cpatch%phos_repro_stoich(ft) = &
               NewRecruitTotalStoichiometry(ft,rec_l2fr_pft,phosphorus_element)
       end do

       ccohort => cpatch%shortest
       cloop: do while(associated(ccohort))
          rec_l2fr_pft = csite%rec_l2fr(ccohort%pft,cl)
          ccohort%nc_repro = NewRecruitTotalStoichiometry(ccohort%pft,rec_l2fr_pft,nitrogen_element)
          ccohort%pc_repro = NewRecruitTotalStoichiometry(ccohort%pft,rec_l2fr_pft,phosphorus_element)
          ccohort => ccohort%taller
       end do cloop
       
       cpatch => cpatch%older
    end do
       
    return
  end subroutine UpdateRecruitStoich

  ! ======================================================================
  
  subroutine SetRecruitL2FR(csite)


    type(ed_site_type) :: csite
    type(fates_patch_type), pointer :: cpatch
    type(fates_cohort_type), pointer :: ccohort
    integer :: ft,cl
    
    if(hlm_parteh_mode .ne. prt_cnp_flex_allom_hyp) return
    
    cpatch => csite%youngest_patch
    do while(associated(cpatch))
       ccohort => cpatch%shortest
       cloop: do while(associated(ccohort))

          if( ccohort%isnew ) then
             ft = ccohort%pft
             cl = ccohort%canopy_layer
             ccohort%l2fr = csite%rec_l2fr(ft,cl)
          end if

          ccohort => ccohort%taller
       end do cloop

       cpatch => cpatch%older
    end do
    
    return
  end subroutine SetRecruitL2FR

end module EDPhysiologyMod
