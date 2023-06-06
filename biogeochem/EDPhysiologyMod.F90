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
  use FatesInterfaceTypesMod, only    : hlm_nitrogen_spec
  use FatesInterfaceTypesMod, only    : hlm_phosphorus_spec
  use FatesInterfaceTypesMod, only    : hlm_use_tree_damage
  use FatesConstantsMod, only    : r8 => fates_r8
  use FatesConstantsMod, only    : nearzero
  use FatesConstantsMod, only    : nocomp_bareground
  use EDPftvarcon      , only    : EDPftvarcon_inst
  use PRTParametersMod , only    : prt_params
  use EDPftvarcon      , only    : GetDecompyFrac
  use FatesInterfaceTypesMod, only    : bc_in_type
  use FatesInterfaceTypesMod, only    : bc_out_type
  use EDCohortDynamicsMod , only : zero_cohort
  use EDCohortDynamicsMod , only : create_cohort, sort_cohorts
  use EDCohortDynamicsMod , only : InitPRTObject
  use EDCohortDynamicsMod , only : InitPRTBoundaryConditions
  use EDCohortDynamicsMod , only : copy_cohort
  use FatesAllometryMod   , only : tree_lai
  use FatesAllometryMod   , only : tree_sai
  use FatesAllometryMod   , only : leafc_from_treelai
  use FatesAllometryMod   , only : decay_coeff_kn
  use FatesLitterMod      , only : litter_type
  use EDTypesMod          , only : site_massbal_type
  use EDTypesMod          , only : numlevsoil_max
  use EDTypesMod          , only : numWaterMem
  use EDTypesMod          , only : dl_sf, dinc_vai, dlower_vai, area_inv
  use EDTypesMod          , only : AREA
  use FatesLitterMod      , only : ncwd
  use FatesLitterMod      , only : ndcmpy
  use FatesLitterMod      , only : ilabile
  use FatesLitterMod      , only : ilignin
  use FatesLitterMod      , only : icellulose
  use FatesLitterMod      , only : adjust_SF_CWD_frac
  use EDTypesMod          , only : nclmax
  use EDTypesMod          , only : AREA,AREA_INV
  use EDTypesMod          , only : nlevleaf
  use EDTypesMod          , only : num_vegtemp_mem
  use EDTypesMod          , only : maxpft
  use EDTypesMod          , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod          , only : leaves_on
  use EDTypesMod          , only : leaves_off
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
  use EDTypesMod          , only : init_recruit_trim
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use EDParamsMod           , only : fates_mortality_disturbance_fraction
  use EDParamsMod           , only : q10_mr
  use EDParamsMod           , only : q10_froz
  use EDParamsMod           , only : logging_export_frac
  use FatesPlantHydraulicsMod  , only : AccumulateMortalityWaterStorage
  use FatesConstantsMod     , only : itrue,ifalse
  use FatesConstantsMod     , only : calloc_abs_error
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
  
  implicit none
  private

  public :: trim_canopy
  public :: phenology
  public :: satellite_phenology
  public :: assign_cohort_SP_properties
  public :: recruitment
  public :: ZeroLitterFluxes

  public :: ZeroAllocationRates
  public :: PreDisturbanceLitterFluxes
  public :: PreDisturbanceIntegrateLitter
  public :: GenerateDamageAndLitterFluxes
  public :: SeedIn
  public :: UpdateRecruitL2FR
  public :: UpdateRecruitStoicH
  public :: SetRecruitL2FR
  
  logical, parameter :: debug  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  integer :: istat           ! return status code
  character(len=255) :: smsg ! Message string for deallocation errors
  integer, parameter :: dleafon_drycheck = 100 ! Drought deciduous leaves max days on check parameter

  ! ============================================================================

contains

  subroutine ZeroLitterFluxes( currentSite )

    ! This routine loops through all patches in a site
    ! and zero's the flux terms for the litter pools.
    ! This is typically called at the beginning of the dynamics
    ! call sequence.


    ! !ARGUMENTS
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), pointer               :: currentPatch

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
    type(ed_patch_type), pointer               :: currentPatch
    type(ed_cohort_type), pointer              :: currentCohort

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
    type(ed_patch_type) :: cpatch
    type(bc_in_type), intent(in) :: bc_in
    

    ! Locals
    type(ed_cohort_type), pointer :: ccohort    ! Current cohort
    type(ed_cohort_type), pointer :: ndcohort   ! New damage-class cohort
    type(litter_type), pointer :: litt     ! Points to the litter object
    type(site_fluxdiags_type), pointer :: flux_diags ! pointer to site level flux diagnostics object
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
             call InitPRTBoundaryConditions(ndcohort)
             call zero_cohort(ndcohort)
             
             ! nc_canopy_d is the new cohort that gets damaged 
             call copy_cohort(ccohort, ndcohort)
             
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
                flux_diags => csite%flux_diags(el)
                
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

                flux_diags%leaf_litter_input(ipft) = &
                     flux_diags%leaf_litter_input(ipft) +  &
                     (store_loss+leaf_loss+repro_loss) * ndcohort%n
                
                call adjust_SF_CWD_frac(ndcohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)

                do c = 1,ncwd
                   litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + &
                        (sapw_loss + struct_loss) * &
                        SF_val_CWD_frac_adj(c) * ndcohort%n / &
                        cpatch%area
                   
                   flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
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
    type(ed_patch_type), intent(inout) :: currentPatch
    type(bc_in_type), intent(in)       :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(site_massbal_type), pointer :: site_mass
    type(litter_type), pointer :: litt     ! Points to the litter object for
    ! the different element types
    integer :: el                          ! Litter element loop index
    integer :: nlev_eff_decomp             ! Number of active layers over which
    ! fragmentation fluxes are transfered
    !------------------------------------------------------------------------------------

    ! Calculate the fragmentation rates
    call fragmentation_scaler(currentPatch, bc_in)

    do el = 1, num_elements

       litt => currentPatch%litter(el)

       ! Calculate loss rate of viable seeds to litter
       call SeedDecay(litt)

       ! Calculate seed germination rate, the status flags prevent
       ! germination from occuring when the site is in a drought
       ! (for drought deciduous) or too cold (for cold deciduous)
       call SeedGermination(litt, currentSite%cstatus, currentSite%dstatus)

       ! Send fluxes from newly created litter into the litter pools
       ! This litter flux is from non-disturbance inducing mortality, as well
       ! as litter fluxes from live trees
       call CWDInput(currentSite, currentPatch, litt,bc_in)

       ! Only calculate fragmentation flux over layers that are active
       ! (RGK-Mar2019) SHOULD WE MAX THIS AT 1? DONT HAVE TO

       nlev_eff_decomp = max(bc_in%max_rooting_depth_index_col, 1)
       call CWDOut(litt,currentPatch%fragmentation_scaler,nlev_eff_decomp)

       site_mass => currentSite%mass_balance(el)

       ! Fragmentation flux to soil decomposition model [kg/site/day]
       site_mass%frag_out = site_mass%frag_out + currentPatch%area * &
            ( sum(litt%ag_cwd_frag) + sum(litt%bg_cwd_frag) + &
            sum(litt%leaf_fines_frag) + sum(litt%root_fines_frag) + &
            sum(litt%seed_decay) + sum(litt%seed_germ_decay))

    end do


    return
  end subroutine PreDisturbanceLitterFluxes

  ! =====================================================================================

  subroutine PreDisturbanceIntegrateLitter(currentPatch)

    ! -----------------------------------------------------------------------------------
    !
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
    type(ed_patch_type),intent(inout),target :: currentPatch


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
    type (ed_cohort_type) , pointer :: currentCohort
    type (ed_patch_type)  , pointer :: currentPatch

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
    integer :: ipatch
    integer :: icohort

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
    !----------------------------------------------------------------------

    ipatch = 1 ! Start counting patches

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))

       ! Add debug diagnstic output to determine which patch
       if (debug) then
          write(fates_log(),*) 'Current patch:', ipatch
          write(fates_log(),*) 'Current patch cohorts:', currentPatch%countcohorts
       endif

       icohort = 1

       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))

          ! Save off the incoming trim
          initial_trim = currentCohort%canopy_trim

          ! Add debug diagnstic output to determine which cohort
          if (debug) then
             write(fates_log(),*) 'Current cohort:', icohort
             write(fates_log(),*) 'Starting canopy trim:', initial_trim
          endif

          trimmed = .false.
          ipft = currentCohort%pft
          call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread,currentCohort%pft,&
               currentCohort%crowndamage, currentCohort%c_area)

          
          leaf_c   = currentCohort%prt%GetState(leaf_organ, carbon12_element)

          currentCohort%treelai = tree_lai(leaf_c, currentCohort%pft, currentCohort%c_area, &
               currentCohort%n, currentCohort%canopy_layer,               &
               currentPatch%canopy_layer_tlai,currentCohort%vcmax25top )

          ! We don't need check on sp mode here since we don't trim_canopy with sp mode
          currentCohort%treesai = tree_sai(currentCohort%pft, &
               currentCohort%dbh, currentCohort%crowndamage,  &
               currentCohort%canopy_trim, &
               currentCohort%c_area, currentCohort%n,currentCohort%canopy_layer,& 
               currentPatch%canopy_layer_tlai, currentCohort%treelai, &
               currentCohort%vcmax25top,0 )  

          currentCohort%nv      = count((currentCohort%treelai+currentCohort%treesai) .gt. dlower_vai(:)) + 1

          if (currentCohort%nv > nlevleaf)then
             write(fates_log(),*) 'nv > nlevleaf',currentCohort%nv, &
                  currentCohort%treelai,currentCohort%treesai, &
                  currentCohort%c_area,currentCohort%n,leaf_c
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

          call bleaf(currentcohort%dbh,ipft,&
               currentCohort%crowndamage, currentcohort%canopy_trim,tar_bl)

          if ( int(prt_params%allom_fmode(ipft)) .eq. 1 ) then
             ! only query fine root biomass if using a fine root allometric model that takes leaf trim into account
             call bfineroot(currentcohort%dbh,ipft,currentcohort%canopy_trim,currentcohort%l2fr,tar_bfr)
             bfr_per_bleaf = tar_bfr/tar_bl
          endif

          ! Identify current canopy layer (cl)
          cl = currentCohort%canopy_layer

          ! PFT-level maximum SLA value, even if under a thick canopy (same units as slatop)
          sla_max = prt_params%slamax(ipft)

          ! Initialize nnu_clai_a
          nnu_clai_a(:,:) = 0._r8
          nnu_clai_b(:,:) = 0._r8

          !Leaf cost vs netuptake for each leaf layer.
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
                kn = decay_coeff_kn(ipft,currentCohort%vcmax25top)
                ! Nscaler value at leaf level z
                nscaler_levleaf = exp(-kn * cumulative_lai)
                ! Sla value at leaf level z after nitrogen profile scaling (m2/gC)
                sla_levleaf = prt_params%slatop(ipft)/nscaler_levleaf

                if(sla_levleaf > sla_max)then
                   sla_levleaf = sla_max
                end if

                !Leaf Cost kgC/m2/year-1
                !decidous costs.
                if (prt_params%season_decid(ipft) ==  itrue .or. &
                     prt_params%stress_decid(ipft) == itrue )then

                   ! Leaf cost at leaf level z accounting for sla profile (kgC/m2)
                   currentCohort%leaf_cost =  1._r8/(sla_levleaf*1000.0_r8)

                   if ( int(prt_params%allom_fmode(ipft)) .eq. 1 ) then
                      ! if using trimmed leaf for fine root biomass allometry, add the cost of the root increment
                      ! to the leaf increment; otherwise do not.
                      currentCohort%leaf_cost = currentCohort%leaf_cost + &
                           1.0_r8/(sla_levleaf*1000.0_r8) * &
                           bfr_per_bleaf / prt_params%root_long(ipft)
                   endif

                   currentCohort%leaf_cost = currentCohort%leaf_cost * &
                        (prt_params%grperc(ipft) + 1._r8)
                else !evergreen costs

                   ! Leaf cost at leaf level z accounting for sla profile
                   currentCohort%leaf_cost = 1.0_r8/(sla_levleaf* &
                        sum(prt_params%leaf_long(ipft,:))*1000.0_r8) !convert from sla in m2g-1 to m2kg-1


                   if ( int(prt_params%allom_fmode(ipft)) .eq. 1 ) then
                      ! if using trimmed leaf for fine root biomass allometry, add the cost of the root increment
                      ! to the leaf increment; otherwise do not.
                      currentCohort%leaf_cost = currentCohort%leaf_cost + &
                           1.0_r8/(sla_levleaf*1000.0_r8) * &
                           bfr_per_bleaf / prt_params%root_long(ipft)
                   endif
                   currentCohort%leaf_cost = currentCohort%leaf_cost * &
                        (prt_params%grperc(ipft) + 1._r8)
                endif

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
                      if (currentCohort%hite > EDPftvarcon_inst%hgt_min(ipft)) then
                         currentCohort%canopy_trim = currentCohort%canopy_trim - &
                              EDPftvarcon_inst%trim_inc(ipft)

                         trimmed = .true.

                      endif ! hite check
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
          icohort = icohort + 1
       enddo
       currentPatch => currentPatch%older
       ipatch = ipatch + 1
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
    use EDParamsMod, only : ED_val_phen_drought_threshold, ED_val_phen_doff_time
    use EDParamsMod, only : ED_val_phen_a, ED_val_phen_b, ED_val_phen_c, ED_val_phen_chiltemp
    use EDParamsMod, only : ED_val_phen_mindayson, ED_val_phen_ncolddayslim, ED_val_phen_coldtemp

    !
    ! !ARGUMENTS:
    type(ed_site_type), intent(inout), target :: currentSite
    type(bc_in_type),   intent(in)            :: bc_in

    !
    ! !LOCAL VARIABLES:

    type(ed_patch_type),pointer :: cpatch
    integer  :: model_day_int     ! integer model day 1 - inf
    integer  :: ncolddays         ! no days underneath the threshold for leaf drop
    integer  :: i_wmem            ! Loop counter for water mem days
    integer  :: i_tmem            ! Loop counter for veg temp mem days
    integer  :: dayssincedleafon  ! Days since drought-decid leaf-on started
    integer  :: dayssincedleafoff ! Days since drought-decid leaf-off started
    integer  :: dayssincecleafon  ! Days since cold-decid leaf-on started
    integer  :: dayssincecleafoff ! Days since cold-decid leaf-off started
    real(r8) :: mean_10day_liqvol ! mean liquid volume (m3/m3) over last 10 days
    real(r8) :: leaf_c            ! leaf carbon [kg]
    real(r8) :: fnrt_c            ! fineroot carbon [kg]
    real(r8) :: sapw_c            ! sapwood carbon [kg]
    real(r8) :: store_c           ! storage carbon [kg]
    real(r8) :: struct_c          ! structure carbon [kg]
    real(r8) :: gdd_threshold     ! GDD accumulation function,
    integer  :: ilayer_swater     ! Layer index for soil water
    ! which also depends on chilling days.
    integer  :: ncdstart          ! beginning of counting period for chilling degree days.
    integer  :: gddstart          ! beginning of counting period for growing degree days.
    real(r8) :: temp_in_C         ! daily averaged temperature in celcius

    integer, parameter :: canopy_leaf_lifespan = 365    ! Maximum lifespan of drought decid leaves

    integer, parameter :: min_daysoff_dforcedflush = 30 ! THis is the number of days that must had elapsed
    ! since leaves had dropped, in order to forcably
    ! flush leaves again.  This does not impact flushing
    ! due to real moisture constraints, and will prevent
    ! drought deciduous in perennially wet environments
    ! that have been forced to drop their leaves, from
    ! flushing them back immediately.

    real(r8),parameter :: dphen_soil_depth = 0.1        ! Use liquid soil water that is
    ! closest to this depth [m]

    ! This is the integer model day. The first day of the simulation is 1, and it
    ! continues monotonically, indefinitely
    ! Advance it. (this should be a global, no reason
    ! for site level, but we don't have global scalars in the
    ! restart file)

    currentSite%phen_model_date = currentSite%phen_model_date + 1
    model_day_int = currentSite%phen_model_date

    ! Use the following layer index to calculate drought conditions
    ilayer_swater = minloc(abs(bc_in%z_sisl(:)-dphen_soil_depth),dim=1)


    ! Parameter of drought decid leaf loss in mm in top layer...FIX(RF,032414)
    ! - this is arbitrary and poorly understood. Needs work. ED_
    !Parameters: defaults from Botta et al. 2000 GCB,6 709-725
    !Parameters, default from from SDGVM model of senesence

    temp_in_C = 0._r8
    cpatch => CurrentSite%oldest_patch
    do while(associated(cpatch))
       temp_in_C = temp_in_C + cpatch%tveg24%GetMean()*cpatch%area
       cpatch => cpatch%younger
    end do
    temp_in_C = temp_in_C * area_inv - tfrz


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
    if(model_day_int>365)then !only do this after the first year to prevent odd behaviour

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
       dayssincecleafoff = model_day_int - (currentSite%cleafoffdate - 365)
    else
       dayssincecleafoff = model_day_int - currentSite%cleafoffdate
    end if

    if (model_day_int < currentSite%cleafondate) then
       dayssincecleafon = model_day_int - (currentSite%cleafondate - 365)
    else
       dayssincecleafon = model_day_int - currentSite%cleafondate
    end if



    !LEAF ON: COLD DECIDUOUS. Needs to
    !1) have exceeded the growing degree day threshold
    !2) The leaves should not be on already
    !3) There should have been at least one chilling day in the counting period.
    !   this prevents tropical or warm climate plants that are "cold-deciduous"
    !   from ever re-flushing after they have reached their maximum age (thus
    !   preventing them from competing

    if ( (currentSite%cstatus == phen_cstat_iscold .or. &
         currentSite%cstatus == phen_cstat_nevercold) .and. &
         (currentSite%grow_deg_days > gdd_threshold) .and. &
         (dayssincecleafoff > ED_val_phen_mindayson) .and. &
         (currentSite%nchilldays >= 1)) then
       currentSite%cstatus = phen_cstat_notcold  ! Set to not-cold status (leaves can come on)
       currentSite%cleafondate = model_day_int
       dayssincecleafon = 0
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
         (dayssincecleafon > ED_val_phen_mindayson) )then

       currentSite%grow_deg_days  = 0._r8          ! The equations for Botta et al
       ! are for calculations of
       ! first flush, but if we dont
       ! clear this value, it will cause
       ! leaves to flush later in the year
       currentSite%cstatus       = phen_cstat_iscold  ! alter status of site to 'leaves off'
       currentSite%cleafoffdate = model_day_int       ! record leaf off date

       if ( debug ) write(fates_log(),*) 'leaves off'
    endif

    ! LEAF OFF: COLD LIFESPAN THRESHOLD
    ! NOTE: Some areas of the planet will never generate a cold day
    ! and thus %nchilldays will never go from zero to 1.  The following logic
    ! when coupled with this fact will essentially prevent cold-deciduous
    ! plants from re-emerging in areas without at least some cold days
    
    if( (currentSite%cstatus == phen_cstat_notcold)  .and. &
         (dayssincecleafoff > 400)) then           ! remove leaves after a whole year,
       ! when there is no 'off' period.
       currentSite%grow_deg_days  = 0._r8

       currentSite%cstatus = phen_cstat_nevercold  ! alter status of site to imply that this
       ! site is never really cold enough
       ! for cold deciduous
       currentSite%cleafoffdate = model_day_int    ! record leaf off date

       if ( debug ) write(fates_log(),*) 'leaves off'
    endif

    !-----------------Drought Phenology--------------------!
    ! Principles of drought-deciduos phenology model...
    ! The 'is_drought' flag is false when leaves are on, and true when leaves area off.
    ! The following sets those site-level flags, which are acted on in phenology_deciduos.
    ! A* The leaves live for either the length of time the soil moisture is over the threshold
    ! or the lifetime of the leaves, whichever is shorter.
    ! B*: If the soil is only wet for a very short time, then the leaves stay on for 100 days
    ! C*: The leaves are only permitted to come ON for a 60 day window around when they last came on,
    ! to prevent 'flickering' on in response to wet season storms
    ! D*: We don't allow anything to happen in the first ten days to allow the water memory window
    ! to come into equlibirium.
    ! E*: If the soil is always wet, the leaves come on at the beginning of the window, and then
    ! last for their lifespan.
    ! ISSUES
    ! 1. It's not clear what water content we should track. Here we are tracking the top layer,
    ! but we probably should track something like BTRAN, but BTRAN is defined for each PFT,
    ! and there could potentially be more than one stress-dec PFT.... ?
    ! 2. In the beginning, the window is set at an arbitrary time of the year, so the leaves
    ! might come on in the dry season, using up stored reserves
    ! for the stress-dec plants, and potentially killing them. To get around this,
    ! we need to read in the 'leaf on' date from some kind of start-up file
    ! but we would need that to happen for every resolution, etc.
    ! 3. Will this methodology properly kill off the stress-dec trees where there is no
    ! water stress? What about where the wet period coincides with the warm period?
    ! We would just get them overlapping with the cold-dec trees, even though that isn't appropriate
    ! Why don't the drought deciduous trees grow in the North?
    ! Is cold decidousness maybe even the same as drought deciduosness there (and so does this
    ! distinction actually matter??)....

    ! Accumulate surface water memory of last 10 days.
    ! Liquid volume in ground layer (m3/m3)
    do i_wmem = 1,numWaterMem-1 !shift memory along one
       currentSite%water_memory(numWaterMem+1-i_wmem) = currentSite%water_memory(numWaterMem-i_wmem)
    enddo
    currentSite%water_memory(1) = bc_in%h2o_liqvol_sl(ilayer_swater)

    ! Calculate the mean water content over the last 10 days (m3/m3)
    mean_10day_liqvol = sum(currentSite%water_memory(1:numWaterMem))/real(numWaterMem,r8)

    ! In drought phenology, we often need to force the leaves to stay
    ! on or off as moisture fluctuates...

    ! Calculate days since leaves have come off, but make a provision
    ! for the first year of simulation, we have to assume a leaf drop
    ! date to start, so if that is in the future, set it to last year

    if (model_day_int < currentSite%dleafoffdate) then
       dayssincedleafoff = model_day_int - (currentSite%dleafoffdate-365)
    else
       dayssincedleafoff = model_day_int - currentSite%dleafoffdate
    endif

    ! the leaves are on. How long have they been on?
    if (model_day_int < currentSite%dleafondate) then
       dayssincedleafon = model_day_int - (currentSite%dleafondate-365)
    else
       dayssincedleafon = model_day_int - currentSite%dleafondate
    endif

    ! LEAF ON: DROUGHT DECIDUOUS WETNESS
    ! Here, we used a window of oppurtunity to determine if we are
    ! close to the time when then leaves came on last year

    ! Has it been ...
    ! a) a year, plus or minus 1 month since we last had leaf-on?
    ! b) Has there also been at least a nominaly short amount of "leaf-off"
    ! c) is the model day at least > 10 (let soil water spin-up)
    ! Note that cold-starts begin in the "leaf-on"
    ! status
    if ( (currentSite%dstatus == phen_dstat_timeoff .or. &
         currentSite%dstatus == phen_dstat_moistoff) .and. &
         (model_day_int > numWaterMem) .and. &
         (dayssincedleafon >= 365-30 .and. dayssincedleafon <= 365+30 ) .and. &
         (dayssincedleafoff > ED_val_phen_doff_time) ) then

       ! If leaves are off, and have been off for at least a few days
       ! and the time is consistent with the correct
       ! time window... test if the moisture conditions allow for leaf-on

       if ( mean_10day_liqvol >= ED_val_phen_drought_threshold ) then
          currentSite%dstatus     = phen_dstat_moiston  ! set status to leaf-on
          currentSite%dleafondate = model_day_int       ! save the model day we start flushing
          dayssincedleafon        = 0
       endif
    endif

    ! LEAF ON: DROUGHT DECIDUOUS TIME EXCEEDANCE
    ! If we still haven't done budburst by end of window, then force it

    ! If the status is "phen_dstat_moistoff", it means this site currently has
    ! leaves off due to actual moisture limitations.
    ! So we trigger bud-burst at the end of the month since
    ! last year's bud-burst.  If this is imposed, then we set the new
    ! status to indicate bud-burst was forced by timing

    if( currentSite%dstatus == phen_dstat_moistoff ) then
       if ( dayssincedleafon > 365+30 ) then
          currentSite%dstatus     = phen_dstat_timeon ! force budburst!
          currentSite%dleafondate = model_day_int     ! record leaf on date
          dayssincedleafon        = 0
       end if
    end if

    ! But if leaves are off due to time, then we enforce
    ! a longer cool-down (because this is a perrenially wet system)

    if(currentSite%dstatus == phen_dstat_timeoff ) then
       if (dayssincedleafoff > min_daysoff_dforcedflush) then
          currentSite%dstatus     = phen_dstat_timeon    ! force budburst!
          currentSite%dleafondate = model_day_int        ! record leaf on date
          dayssincedleafon        = 0
       end if
    end if

    ! LEAF OFF: DROUGHT DECIDUOUS LIFESPAN - if the leaf gets to
    ! the end of its useful life. A*, E*
    ! i.e. Are the leaves rouhgly at the end of their lives?

    if ( (currentSite%dstatus == phen_dstat_moiston .or. &
         currentSite%dstatus == phen_dstat_timeon ) .and. &
         (dayssincedleafon > canopy_leaf_lifespan) )then
       currentSite%dstatus      = phen_dstat_timeoff    !alter status of site to 'leaves off'
       currentSite%dleafoffdate = model_day_int         !record leaf on date
    endif

    ! LEAF OFF: DROUGHT DECIDUOUS DRYNESS - if the soil gets too dry,
    ! and the leaves have already been on a while...

    if ( (currentSite%dstatus == phen_dstat_moiston .or. &
         currentSite%dstatus == phen_dstat_timeon ) .and. &
         (model_day_int > numWaterMem) .and. &
         (mean_10day_liqvol <= ED_val_phen_drought_threshold) .and. &
         (dayssincedleafon > dleafon_drycheck ) ) then
       currentSite%dstatus = phen_dstat_moistoff     ! alter status of site to 'leaves off'
       currentSite%dleafoffdate = model_day_int      ! record leaf on date
    endif

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
    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

    real(r8) :: leaf_c                 ! leaf carbon [kg]
    real(r8) :: sapw_c                 ! sapwood carbon [kg]
    real(r8) :: struct_c               ! structural wood carbon [kg]
    real(r8) :: store_c                ! storage carbon [kg]
    real(r8) :: store_c_transfer_frac  ! Fraction of storage carbon used to flush leaves
    real(r8) :: deficit_c              ! Amount of C needed to get flushing pools "on-allometry"
    real(r8) :: target_leaf_c
    real(r8) :: target_sapw_c
    real(r8) :: target_agw_c, target_bgw_c, target_struct_c
    real(r8) :: sapw_area
    integer  :: ipft
    real(r8), parameter :: leaf_drop_fraction = 1.0_r8
    real(r8), parameter :: carbon_store_buffer = 0.10_r8
    real(r8) :: stem_drop_fraction
    !------------------------------------------------------------------------

    currentPatch => CurrentSite%oldest_patch

    do while(associated(currentPatch))
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))

          ipft = currentCohort%pft

          ! Retrieve existing leaf and storage carbon

          if(debug) call currentCohort%prt%CheckMassConservation(ipft,0)

          store_c = currentCohort%prt%GetState(store_organ, carbon12_element)
          leaf_c  = currentCohort%prt%GetState(leaf_organ, carbon12_element)
          sapw_c  = currentCohort%prt%GetState(sapw_organ, carbon12_element)
          struct_c  = currentCohort%prt%GetState(struct_organ, carbon12_element)

          stem_drop_fraction = EDPftvarcon_inst%phen_stem_drop_fraction(ipft)

          ! COLD LEAF ON
          ! The site level flags signify that it is no-longer too cold
          ! for leaves. Time to signal flushing

          if_colddec: if (prt_params%season_decid(ipft) == itrue)then
             if_notcold: if ( currentSite%cstatus == phen_cstat_notcold  )then                ! we have just moved to leaves being on .
                if_leaves_off: if (currentCohort%status_coh == leaves_off)then ! Are the leaves currently off?
                   currentCohort%status_coh = leaves_on         ! Leaves are on, so change status to
                   ! stop flow of carbon out of bstore.

                   call bleaf(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage, &
                        currentCohort%canopy_trim,target_leaf_c)
                   call bsap_allom(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage, &
                        currentCohort%canopy_trim,sapw_area,target_sapw_c)
                   call bagw_allom(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage,&
                        target_agw_c)
                   call bbgw_allom(currentCohort%dbh,currentCohort%pft,target_bgw_c)
                   call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, &
                        currentCohort%pft, target_struct_c)

                   if (stem_drop_fraction .gt. 0.0_r8) then
                      ! Note, this is only true for some grasses, woody plants don't
                      ! have a stem drop fraction
                      deficit_c = target_leaf_c + (target_sapw_c-sapw_c) + (target_struct_c-struct_c)
                   else
                      deficit_c = target_leaf_c
                   end if
                   
                   if(store_c>nearzero) then

                      ! flush either the amount to get to the target, or -most- of the storage pool
                      ! RF: added a criterion to stop the entire store pool emptying and triggering termination mortality
                      ! n.b. this might not be necessary if we adopted a more gradual approach to leaf flushing...
                      store_c_transfer_frac = min((EDPftvarcon_inst%phenflush_fraction(ipft)*deficit_c)/store_c, &
                                                  (1.0_r8-carbon_store_buffer))

                   else
                      store_c_transfer_frac = 0.0_r8
                   end if

                   ! This call will request that storage carbon will be transferred to
                   ! leaf tissues. It is specified as a fraction of the available storage
                   ! Check that the stem drop fraction is set to non-zero amount
                   ! otherwise flush all carbon store to leaves
                   if (stem_drop_fraction .gt. 0.0_r8) then
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                           store_c_transfer_frac*target_leaf_c/deficit_c)
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, sapw_organ, &
                           store_c_transfer_frac*(target_sapw_c-sapw_c)/deficit_c)
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, struct_organ, &
                           store_c_transfer_frac*(target_struct_c-struct_c)/deficit_c)
                      
                   else
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                           store_c_transfer_frac)
                      
                   end if

                endif if_leaves_off
             endif if_notcold

             !COLD LEAF OFF
             if_cold:  if (currentSite%cstatus == phen_cstat_nevercold .or. &
                  currentSite%cstatus == phen_cstat_iscold) then ! past leaf drop day? Leaves still on tree?
                
                if_leaves_on: if (currentCohort%status_coh == leaves_on) then ! leaves have not dropped
                   
                   ! leaf off occur on individuals bigger than specific size for grass
                   if (currentCohort%dbh > EDPftvarcon_inst%phen_cold_size_threshold(ipft) &
                        .or. prt_params%woody(ipft)==itrue) then

                      ! This sets the cohort to the "leaves off" flag
                      currentCohort%status_coh  = leaves_off

                      ! Drop Leaves (this routine will update the leaf state variables,
                      ! for carbon and any other element that are prognostic. It will
                      ! also track the turnover masses that will be sent to litter later on)

                      call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                           leaf_organ, leaf_drop_fraction)

                      if(prt_params%woody(ipft).ne.itrue)then

                         call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                              sapw_organ, stem_drop_fraction)

                         call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                              struct_organ, stem_drop_fraction)

                      endif	! woody plant check
                   endif ! individual dbh size check
                endif if_leaves_on !leaf status
             endif if_cold !currentSite status
          endif if_colddec  !season_decid

          ! DROUGHT LEAF ON
          ! Site level flag indicates it is no longer in drought condition
          ! deciduous plants can flush

          if (prt_params%stress_decid(ipft) == itrue )then

             if (currentSite%dstatus == phen_dstat_moiston .or. &
                  currentSite%dstatus == phen_dstat_timeon )then

                ! we have just moved to leaves being on .
                if (currentCohort%status_coh == leaves_off)then

                   !is it the leaf-on day? Are the leaves currently off?

                   currentCohort%status_coh = leaves_on    ! Leaves are on, so change status to
                   ! stop flow of carbon out of bstore.

                   call bleaf(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage,&
                        currentCohort%canopy_trim,target_leaf_c)
                   call bsap_allom(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage,&
                        currentCohort%canopy_trim,sapw_area,target_sapw_c)
                   call bagw_allom(currentCohort%dbh,currentCohort%pft,currentCohort%crowndamage,&
                        target_agw_c)
                   call bbgw_allom(currentCohort%dbh,currentCohort%pft,target_bgw_c)
                   call bdead_allom( target_agw_c, target_bgw_c, target_sapw_c, &
                        currentCohort%pft, target_struct_c)

                   if (stem_drop_fraction .gt. 0.0_r8) then
                      ! Note, this is only true for some grasses, woody plants don't
                      ! have a stem drop fraction
                      deficit_c = target_leaf_c + (target_sapw_c-sapw_c) + (target_struct_c-struct_c)
                   else
                      deficit_c = target_leaf_c
                   end if
                   
                   if(store_c>nearzero) then
                      
                      store_c_transfer_frac = &
                           min((EDPftvarcon_inst%phenflush_fraction(ipft)*deficit_c)/store_c, &
                           (1.0_r8-carbon_store_buffer))

                   else
                      store_c_transfer_frac = 0.0_r8
                   endif

                   ! This call will request that storage carbon will be transferred to
                   ! leaf tissues. It is specified as a fraction of the available storage
                   if (stem_drop_fraction .gt. 0.0_r8) then
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                           store_c_transfer_frac*target_leaf_c/deficit_c)
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, sapw_organ, &
                           store_c_transfer_frac*(target_sapw_c-sapw_c)/deficit_c)
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, struct_organ, &
                           store_c_transfer_frac*(target_struct_c-struct_c)/deficit_c)
                      
                   else
                      
                      call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                           store_c_transfer_frac)
                      
                   end if

                endif !currentCohort status again?
             endif   !currentSite status

             !DROUGHT LEAF OFF
             if (currentSite%dstatus == phen_dstat_moistoff .or. &
                  currentSite%dstatus == phen_dstat_timeoff) then

                if (currentCohort%status_coh == leaves_on) then ! leaves have not dropped

                   ! This sets the cohort to the "leaves off" flag
                   currentCohort%status_coh      = leaves_off

                   call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                        leaf_organ, leaf_drop_fraction)

                   if(prt_params%woody(ipft).ne.itrue)then

                      call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                           sapw_organ, stem_drop_fraction)

                      call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                           struct_organ, stem_drop_fraction)
                   endif

                endif
             endif !status
          endif !drought dec.

          if(debug) call currentCohort%prt%CheckMassConservation(ipft,1)

          currentCohort => currentCohort%shorter
       enddo !currentCohort

       currentPatch => currentPatch%younger

    enddo !currentPatch

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
    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort

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

  ! =====================================================================================

  subroutine assign_cohort_SP_properties(currentCohort,htop,tlai,tsai,parea,init,leaf_c)

    ! -----------------------------------------------------------------------------------!
    ! Takes the daily inputs of leaf area index, stem area index and canopy height and
    ! translates them into a FATES structure with one patch and one cohort per PFT
    ! The leaf area of the cohort is modified each day to match that asserted by the HLM
    ! -----------------------------------------------------------------------------------!
   

    type(ed_cohort_type), intent(inout), target :: currentCohort

    real(r8), intent(in) :: tlai ! target leaf area index from SP inputs
    real(r8), intent(in) :: tsai ! target stem area index from SP inputs
    real(r8), intent(in) :: htop ! target tree height from SP inputs
    real(r8), intent(in) :: parea ! patch area for this PFT
    integer, intent(in)  :: init ! are we in the initialization routine? if so do not set leaf_c
    real(r8), intent(out) ::  leaf_c        ! leaf carbon estimated to generate target tlai

    real(r8) :: dummy_n       ! set cohort n to a dummy value of 1.0
    integer  :: fates_pft     ! fates pft numer for weighting loop
    real(r8) :: spread        ! dummy value of canopy spread to estimate c_area
    real(r8) :: check_treelai
    real(r8) :: canopylai(1:nclmax)
    real(r8) :: fracerr
    real(r8) :: oldcarea

    ! Do some checks
    if(associated(currentCohort%shorter))then
       write(fates_log(),*) 'SP mode has >1 cohort'
       write(fates_log(),*) "SP mode >1 cohort: PFT",currentCohort%pft, currentCohort%shorter%pft
       write(fates_log(),*) "SP mode >1 cohort: CL",currentCohort%canopy_layer, currentCohort%shorter%canopy_layer
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    !------------------------------------------
    !  Calculate dbh from input height, and c_area from dbh
    !------------------------------------------
    currentCohort%hite = htop

    fates_pft = currentCohort%pft
    call h2d_allom(currentCohort%hite,fates_pft,currentCohort%dbh)

    dummy_n = 1.0_r8 ! make n=1 to get area of one tree.
    spread = 1.0_r8  ! fix this to 0 to remove dynamics of canopy closure, assuming a closed canopy.
    ! n.b. the value of this will only affect 'n', which isn't/shouldn't be a diagnostic in
    ! SP mode.
    call carea_allom(currentCohort%dbh,dummy_n,spread,currentCohort%pft,&
         currentCohort%crowndamage,currentCohort%c_area)

    !------------------------------------------
    !  Calculate canopy N assuming patch area is full
    !------------------------------------------
    currentCohort%n = parea / currentCohort%c_area

    ! correct c_area for the new nplant
    call carea_allom(currentCohort%dbh,currentCohort%n,spread,currentCohort%pft,&
         currentCohort%crowndamage,currentCohort%c_area)

    ! ------------------------------------------
    ! Calculate leaf carbon from target treelai
    ! ------------------------------------------
    currentCohort%treelai = tlai
    canopylai(:) = 0._r8
    if(init.eq.itrue)then
       ! If we are initializing, the canopy layer has not been set yet, so just set to 1
       currentCohort%canopy_layer = 1
       ! We need to get the vcmax25top
       currentCohort%vcmax25top = EDPftvarcon_inst%vcmax25top(currentCohort%pft,1)
    endif
    leaf_c = leafc_from_treelai( currentCohort%treelai, currentCohort%pft, currentCohort%c_area,&
         currentCohort%n, currentCohort%canopy_layer, currentCohort%vcmax25top)

    !check that the inverse calculation of leafc from treelai is the same as the
    ! standard calculation of treelai from leafc. Maybe can delete eventually?

    check_treelai = tree_lai(leaf_c, currentCohort%pft, currentCohort%c_area, &
         currentCohort%n, currentCohort%canopy_layer,               &
         canopylai,currentCohort%vcmax25top )

    if( abs(currentCohort%treelai-check_treelai).gt.1.0e-12)then !this is not as precise as nearzero
       write(fates_log(),*) 'error in validate treelai',currentCohort%treelai,check_treelai,currentCohort%treelai-check_treelai
       write(fates_log(),*) 'tree_lai inputs: ', currentCohort%pft, currentCohort%c_area, currentCohort%n, &
               currentCohort%canopy_layer, currentCohort%vcmax25top
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    ! the carea_allom routine sometimes generates precision-tolerance level errors in the canopy area
    ! these mean that the canopy area does not exactly add up to the patch area, which causes chaos in
    ! the radiation routines.  Correct both the area and the 'n' to remove error, and don't use
    !! carea_allom in SP mode after this point.

    if(abs(currentCohort%c_area-parea).gt.nearzero)then ! there is an error
       if(abs(currentCohort%c_area-parea).lt.10.e-9)then !correct this if it's a very small error
          oldcarea = currentCohort%c_area
          !generate new cohort area
          currentCohort%c_area = currentCohort%c_area - (currentCohort%c_area- parea)
          currentCohort%n = currentCohort%n * (currentCohort%c_area/oldcarea)
          if(abs(currentCohort%c_area-parea).gt.nearzero)then
             write(fates_log(),*) 'SPassign, c_area still broken',currentCohort%c_area-parea,currentCohort%c_area-oldcarea
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if
       else
          write(fates_log(),*) 'SPassign, big error in c_area',currentCohort%c_area-parea,currentCohort%pft
       end if ! still broken
    end if !small error

    if(init.eq.ifalse)then
       call SetState(currentCohort%prt, leaf_organ, carbon12_element, leaf_c, 1)
    endif

    ! assert sai
    currentCohort%treesai = tsai

  end subroutine assign_cohort_SP_properties

  ! =====================================================================================

  subroutine SeedIn( currentSite, bc_in )

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
    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout), target  :: currentSite
    type(bc_in_type), intent(in)               :: bc_in

    type(ed_patch_type), pointer     :: currentPatch
    type(litter_type), pointer       :: litt
    type(ed_cohort_type), pointer    :: currentCohort
    type(site_massbal_type), pointer :: site_mass

    integer  :: pft
    real(r8) :: store_m_to_repro       ! mass sent from storage to reproduction upon death [kg/plant]
    real(r8) :: site_seed_rain(maxpft) ! This is the sum of seed-rain for the site [kg/site/day]
    real(r8) :: seed_in_external       ! Mass of externally generated seeds [kg/m2/day]
    real(r8) :: seed_stoich            ! Mass ratio of nutrient per C12 in seeds [kg/kg]
    real(r8) :: seed_prod              ! Seed produced in this dynamics step [kg/day]
    integer  :: n_litt_types           ! number of litter element types (c,n,p, etc)
    integer  :: el                     ! loop counter for litter element types
    integer  :: element_id             ! element id consistent with parteh/PRTGenericMod.F90
    !------------------------------------------------------------------------------------

    do el = 1, num_elements

       site_seed_rain(:) = 0._r8

       element_id = element_list(el)

       site_mass => currentSite%mass_balance(el)

       ! Loop over all patches and sum up the seed input for each PFT
       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))

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
                  (seed_prod * currentCohort%n + store_m_to_repro)

             currentCohort => currentCohort%shorter
          enddo !cohort loop

          currentPatch => currentPatch%younger
       enddo

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
       do while (associated(currentPatch))

          litt => currentPatch%litter(el)
          do pft = 1,numpft

             if(currentSite%use_this_pft(pft).eq.itrue)then
                ! Seed input from local sources (within site)
                litt%seed_in_local(pft) = litt%seed_in_local(pft) + site_seed_rain(pft)/area

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
                seed_in_external =  seed_stoich*EDPftvarcon_inst%seed_suppl(pft)*years_per_day
                litt%seed_in_extern(pft) = litt%seed_in_extern(pft) + seed_in_external

                ! Seeds entering externally [kg/site/day]
                site_mass%seed_in = site_mass%seed_in + seed_in_external*currentPatch%area
             end if !use this pft
          enddo


          currentPatch => currentPatch%younger
       enddo

    end do

    return
  end subroutine SeedIn

  ! ============================================================================

  subroutine SeedDecay( litt )
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into leaf litter pool
    !
    ! !ARGUMENTS
    type(litter_type) :: litt
    !
    ! !LOCAL VARIABLES:
    integer  ::  pft
    !----------------------------------------------------------------------

    ! default value from Liscke and Loffler 2006 ; making this a PFT-specific parameter
    ! decays the seed pool according to exponential model
    ! seed_decay_rate is in yr-1
    ! seed_decay is kg/day
    ! Assume that decay rates are same for all chemical species

    do pft = 1,numpft
       litt%seed_decay(pft) = litt%seed(pft) * &
            EDPftvarcon_inst%seed_decay_rate(pft)*years_per_day

       litt%seed_germ_decay(pft) = litt%seed_germ(pft) * &
            EDPftvarcon_inst%seed_decay_rate(pft)*years_per_day

    enddo

    return
  end subroutine SeedDecay

  ! ============================================================================
  subroutine SeedGermination( litt, cold_stat, drought_stat )
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into sapling pool
    !
    ! !USES:

    !
    ! !ARGUMENTS
    type(litter_type) :: litt
    integer, intent(in) :: cold_stat    ! Is the site in cold leaf-off status?
    integer, intent(in) :: drought_stat ! Is the site in drought leaf-off status?
    !
    ! !LOCAL VARIABLES:
    integer :: pft


    real(r8), parameter ::  max_germination = 1.0_r8 ! Cap on germination rates.
    ! KgC/m2/yr Lishcke et al. 2009

    ! Turning of this cap? because the cap will impose changes on proportionality
    ! of nutrients. (RGK 02-2019)
    !real(r8), parameter :: max_germination = 1.e6_r8  ! Force to very high number

    !----------------------------------------------------------------------

    ! germination_rate is being pulled to PFT parameter; units are 1/yr
    ! thus the mortality rate of seed -> recruit (in units of carbon)
    ! is seed_decay_rate(p)/germination_rate(p)
    ! and thus the mortality rate (in units of individuals) is the product of
    ! that times the ratio of (hypothetical) seed mass to recruit biomass

    do pft = 1,numpft
       litt%seed_germ_in(pft) =  min(litt%seed(pft) * EDPftvarcon_inst%germination_rate(pft), &
            max_germination)*years_per_day

       !set the germination only under the growing season...c.xu

       if ((prt_params%season_decid(pft) == itrue ) .and. &
            (any(cold_stat == [phen_cstat_nevercold,phen_cstat_iscold]))) then
          litt%seed_germ_in(pft) = 0.0_r8
       endif
       if ((prt_params%stress_decid(pft) == itrue ) .and. &
            (any(drought_stat == [phen_dstat_timeoff,phen_dstat_moistoff]))) then
          litt%seed_germ_in(pft) = 0.0_r8
       end if


    enddo

  end subroutine SeedGermination

  ! =====================================================================================





  ! =====================================================================================

  subroutine recruitment( currentSite, currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! spawn new cohorts of juveniles of each PFT
    !
    ! !USES:
    use FatesInterfaceTypesMod, only : hlm_use_ed_prescribed_phys
    use FatesLitterMod   , only : ncwd
    
    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout)  :: currentSite
    type(ed_patch_type), intent(inout),pointer :: currentPatch
    type(bc_in_type), intent(in)       :: bc_in
    !
    ! !LOCAL VARIABLES:
    class(prt_vartypes), pointer :: prt
    integer :: ft
    integer :: c 
    type (ed_cohort_type) , pointer :: temp_cohort
    type (litter_type), pointer     :: litt          ! The litter object (carbon right now)
    type(site_massbal_type), pointer :: site_mass    ! For accounting total in-out mass fluxes
    integer :: cohortstatus
    integer :: el          ! loop counter for element
    integer :: element_id  ! element index consistent with definitions in PRTGenericMod
    integer :: iage        ! age loop counter for leaf age bins
    integer :: crowndamage
    integer,parameter :: recruitstatus = 1 !weather it the new created cohorts is recruited or initialized
    real(r8) :: c_leaf      ! target leaf biomass [kgC]
    real(r8) :: c_fnrt      ! target fine root biomass [kgC]
    real(r8) :: c_sapw      ! target sapwood biomass [kgC]
    real(r8) :: a_sapw      ! target sapwood cross section are [m2] (dummy)
    real(r8) :: c_agw       ! target Above ground biomass [kgC]
    real(r8) :: c_bgw       ! target Below ground biomass [kgC]
    real(r8) :: c_struct    ! target Structural biomass [kgc]
    real(r8) :: c_store     ! target Storage biomass [kgC]
    real(r8) :: m_leaf      ! leaf mass (element agnostic) [kg]
    real(r8) :: m_fnrt      ! fine-root mass (element agnostic) [kg]
    real(r8) :: m_sapw      ! sapwood mass (element agnostic) [kg]
    real(r8) :: m_agw       ! AG wood mass (element agnostic) [kg]
    real(r8) :: m_bgw       ! BG wood mass (element agnostic) [kg]
    real(r8) :: m_struct    ! structural mass (element agnostic) [kg]
    real(r8) :: m_store     ! storage mass (element agnostic) [kg]
    real(r8) :: m_repro     ! reproductive mass (element agnostic) [kg]
    real(r8) :: mass_avail  ! The mass of each nutrient/carbon available in the seed_germination pool [kg]
    real(r8) :: mass_demand ! Total mass demanded by the plant to achieve the stoichiometric targets
    ! of all the organs in the recruits. Used for both [kg per plant] and [kg per cohort]
    real(r8) :: stem_drop_fraction

    !----------------------------------------------------------------------

    allocate(temp_cohort) ! create temporary cohort
    call zero_cohort(temp_cohort)


    do ft = 1,numpft

       ! The following if block is for the prescribed biogeography and/or nocomp modes.
       ! Since currentSite%use_this_pft is a site-level quantity and thus only limits whether a given PFT
       ! is permitted on a given gridcell or not, it applies to the prescribed biogeography case only.
       ! If nocomp is enabled, then we must determine whether a given PFT is allowed on a given patch or not.

       if(currentSite%use_this_pft(ft).eq.itrue &
            .and. ((hlm_use_nocomp .eq. ifalse) .or. (ft .eq. currentPatch%nocomp_pft_label)))then

          temp_cohort%canopy_trim = init_recruit_trim
          temp_cohort%pft         = ft
          temp_cohort%hite        = EDPftvarcon_inst%hgt_min(ft)
          temp_cohort%coage       = 0.0_r8
          stem_drop_fraction      = EDPftvarcon_inst%phen_stem_drop_fraction(ft)
          temp_cohort%l2fr        = currentSite%rec_l2fr(ft,currentPatch%NCL_p)
          temp_cohort%crowndamage = 1       ! new recruits are undamaged
          
          call h2d_allom(temp_cohort%hite,ft,temp_cohort%dbh)

       
          ! Initialize live pools
          call bleaf(temp_cohort%dbh,ft,temp_cohort%crowndamage,&
               temp_cohort%canopy_trim,c_leaf)
          call bfineroot(temp_cohort%dbh,ft,temp_cohort%canopy_trim,temp_cohort%l2fr,c_fnrt)
          call bsap_allom(temp_cohort%dbh,ft,temp_cohort%crowndamage, &
               temp_cohort%canopy_trim,a_sapw, c_sapw)
          call bagw_allom(temp_cohort%dbh,ft,temp_cohort%crowndamage, c_agw)
          call bbgw_allom(temp_cohort%dbh,ft,c_bgw)
          call bdead_allom(c_agw,c_bgw,c_sapw,ft,c_struct)
          call bstore_allom(temp_cohort%dbh,ft, temp_cohort%crowndamage, &
               temp_cohort%canopy_trim,c_store)

          ! Default assumption is that leaves are on
          cohortstatus = leaves_on

          ! But if the plant is seasonally (cold) deciduous, and the site status is flagged
          ! as "cold", then set the cohort's status to leaves_off, and remember the leaf biomass
          if ((prt_params%season_decid(ft) == itrue) .and. &
               (any(currentSite%cstatus == [phen_cstat_nevercold,phen_cstat_iscold]))) then
             c_leaf = 0.0_r8

             ! If plant is not woody then set sapwood and structural biomass as well
             if (prt_params%woody(ft).ne.itrue) then
                c_sapw = (1.0_r8 - stem_drop_fraction) * c_sapw
                c_struct = (1.0_r8 - stem_drop_fraction) * c_struct
             endif
             cohortstatus = leaves_off
          endif

          ! Or.. if the plant is drought deciduous, and the site status is flagged as
          ! "in a drought", then likewise, set the cohort's status to leaves_off, and remember leaf
          ! biomass
          if ((prt_params%stress_decid(ft) == itrue) .and. &
               (any(currentSite%dstatus == [phen_dstat_timeoff,phen_dstat_moistoff]))) then
             c_leaf = 0.0_r8

             ! If plant is not woody then set sapwood and structural biomass as well
             if(prt_params%woody(ft).ne.itrue)then
                c_sapw = (1.0_r8 - stem_drop_fraction) * c_sapw
                c_struct = (1.0_r8 - stem_drop_fraction) * c_struct
             endif
             cohortstatus = leaves_off
          endif


          ! Cycle through available carbon and nutrients, find the limiting element
          ! to dictate the total number of plants that can be generated

          if ( (hlm_use_ed_prescribed_phys .eq. ifalse) .or. &
               (EDPftvarcon_inst%prescribed_recruitment(ft) .lt. 0._r8) ) then

           temp_cohort%n = 1.e20_r8

             do el = 1,num_elements

                element_id = element_list(el)
                select case(element_id)
                case(carbon12_element)

                  mass_demand = c_struct+c_leaf+c_fnrt+c_sapw+c_store

                case(nitrogen_element)

                     mass_demand = &
                          c_struct*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(struct_organ)) + &
                          c_leaf*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)) + &
                          c_fnrt*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)) + &
                          c_sapw*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)) + &
                          StorageNutrientTarget(ft, element_id, &
                          c_leaf*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)), &
                          c_fnrt*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)), &
                          c_sapw*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)), &
                          c_struct*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(struct_organ)))

                case(phosphorus_element)

                  mass_demand = &
                       c_struct*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(struct_organ)) + &
                       c_leaf*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)) + &
                       c_fnrt*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)) + &
                       c_sapw*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)) + &
                       StorageNutrientTarget(ft, element_id, &
                       c_leaf*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(leaf_organ)), &
                       c_fnrt*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ)), &
                       c_sapw*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(sapw_organ)), &
                       c_struct*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(struct_organ)))

                case default
                   write(fates_log(),*) 'Undefined element type in recruitment'
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end select

                mass_avail = currentPatch%area * currentPatch%litter(el)%seed_germ(ft)

                ! ------------------------------------------------------------------------
                ! Update number density if this is the limiting mass
                ! ------------------------------------------------------------------------

                temp_cohort%n = min(temp_cohort%n, mass_avail/mass_demand)

             end do


          else
             ! prescribed recruitment rates. number per sq. meter per year
             temp_cohort%n  = currentPatch%area * &
                  EDPftvarcon_inst%prescribed_recruitment(ft) * &
                  hlm_freq_day
          endif

          ! Only bother allocating a new cohort if there is a reasonable amount of it
          any_recruits: if (temp_cohort%n > min_n_safemath )then

             ! -----------------------------------------------------------------------------
             ! PART II.
             ! Initialize the PARTEH object, and determine the initial masses of all
             ! organs and elements.
             ! -----------------------------------------------------------------------------
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

                 m_struct = c_struct*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(struct_organ))
                 m_leaf   = c_leaf*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(leaf_organ))
                 m_fnrt   = c_fnrt*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ))
                 m_sapw   = c_sapw*prt_params%nitr_stoich_p1(ft,prt_params%organ_param_id(sapw_organ))
                 m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct )
                   m_repro  = 0._r8

                case(phosphorus_element)

                 m_struct = c_struct*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(struct_organ))
                 m_leaf   = c_leaf*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(leaf_organ))
                 m_fnrt   = c_fnrt*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(fnrt_organ))
                 m_sapw   = c_sapw*prt_params%phos_stoich_p1(ft,prt_params%organ_param_id(sapw_organ))
                 m_store  = StorageNutrientTarget(ft, element_id, m_leaf, m_fnrt, m_sapw, m_struct )
                   m_repro  = 0._r8

                end select

                select case(hlm_parteh_mode)
                case (prt_carbon_allom_hyp,prt_cnp_flex_allom_hyp )

                   ! Put all of the leaf mass into the first bin
                   call SetState(prt,leaf_organ, element_id,m_leaf,1)
                   do iage = 2,nleafage
                      call SetState(prt,leaf_organ, element_id,0._r8,iage)
                   end do

                   call SetState(prt,fnrt_organ, element_id, m_fnrt)
                   call SetState(prt,sapw_organ, element_id, m_sapw)
                   call SetState(prt,store_organ, element_id, m_store)
                   call SetState(prt,struct_organ, element_id, m_struct)
                   call SetState(prt,repro_organ, element_id, m_repro)

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

                if ( (hlm_use_ed_prescribed_phys .eq. itrue ) .and. &
                     (EDPftvarcon_inst%prescribed_recruitment(ft) .ge. 0._r8 )) then

                   site_mass%flux_generic_in = site_mass%flux_generic_in + &
                        temp_cohort%n*(m_struct + m_leaf + m_fnrt + m_sapw + m_store + m_repro)

                   site_mass%flux_generic_out = site_mass%flux_generic_out + &
                        currentPatch%area * currentPatch%litter(el)%seed_germ(ft)

                   currentPatch%litter(el)%seed_germ(ft) = 0._r8


                else

                   currentPatch%litter(el)%seed_germ(ft) = currentPatch%litter(el)%seed_germ(ft) - &
                        temp_cohort%n / currentPatch%area * &
                        (m_struct + m_leaf + m_fnrt + m_sapw + m_store + m_repro)

                end if



             end do

             ! This call cycles through the initial conditions, and makes sure that they
             ! are all initialized.
             ! -----------------------------------------------------------------------------------

             call prt%CheckInitialConditions()

             ! This initializes the cohort

             call create_cohort(currentSite,currentPatch, temp_cohort%pft, temp_cohort%n, &
                  temp_cohort%hite, temp_cohort%coage, temp_cohort%dbh, prt, &
                  cohortstatus, recruitstatus, &
                  temp_cohort%canopy_trim,temp_cohort%c_area, &
                  currentPatch%NCL_p, &
                  temp_cohort%crowndamage, &
                  currentSite%spread, bc_in)

             ! Note that if hydraulics is on, the number of cohorts may had
             ! changed due to hydraulic constraints.
             ! This constaint is applied during "create_cohort" subroutine.

             ! keep track of how many individuals were recruited for passing to history
             currentSite%recruitment_rate(ft) = currentSite%recruitment_rate(ft) + temp_cohort%n


        endif any_recruits
       endif !use_this_pft
    enddo  !pft loop

    deallocate(temp_cohort, stat=istat, errmsg=smsg)
    if (istat/=0) then
       write(fates_log(),*) 'dealloc013: fail on deallocate(temp_cohort):'//trim(smsg)
       call endrun(msg=errMsg(sourcefile, __LINE__))
    endif

    
  end subroutine recruitment

  ! ============================================================================

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

    !
    ! !ARGUMENTS
    type(ed_site_type), intent(inout), target :: currentSite
    type(ed_patch_type),intent(inout), target :: currentPatch
    type(litter_type),intent(inout),target    :: litt
    type(bc_in_type),intent(in)               :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer      :: currentCohort
    type(site_fluxdiags_type), pointer :: flux_diags
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
    !----------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Other direct litter fluxes happen in phenology and in spawn_patches.
    ! -----------------------------------------------------------------------------------

    numlevsoil = currentSite%nlevsoil

    element_id = litt%element_id

    ! Object tracking flux diagnostics for each element
    flux_diags => currentSite%flux_diags(element_pos(element_id))

    ! Object tracking site level mass balance for each element
    site_mass => currentSite%mass_balance(element_pos(element_id))

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

       flux_diags%leaf_litter_input(pft) = &
            flux_diags%leaf_litter_input(pft) +  &
            leaf_m_turnover * currentCohort%n

       root_fines_tot = (fnrt_m_turnover + store_m_turnover ) * &
            plant_dens

       do dcmpy=1,ndcmpy
          dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)
          litt%leaf_fines_in(dcmpy) = litt%leaf_fines_in(dcmpy) + &
               (leaf_m_turnover+repro_m_turnover) * plant_dens * dcmpy_frac

          dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
          do ilyr = 1, numlevsoil
             litt%root_fines_in(dcmpy,ilyr) = litt%root_fines_in(dcmpy,ilyr) + &
                  currentSite%rootfrac_scr(ilyr) * root_fines_tot * dcmpy_frac
          end do
       end do

       flux_diags%root_litter_input(pft) = &
            flux_diags%root_litter_input(pft) +  &
            (fnrt_m_turnover + store_m_turnover ) * currentCohort%n


       ! Assumption: turnover from deadwood and sapwood are lumped together in CWD pool

       !update partitioning of stem wood (struct + sapw) to cwd based on cohort dbh
       call adjust_SF_CWD_frac(currentCohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj)


       do c = 1,ncwd
          litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + &
               (sapw_m_turnover + struct_m_turnover) * &
               SF_val_CWD_frac_adj(c) * plant_dens * &
               prt_params%allom_agb_frac(pft)

          flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
               (struct_m_turnover + sapw_m_turnover) * SF_val_CWD_frac_adj(c) * &
               prt_params%allom_agb_frac(pft) * currentCohort%n

          bg_cwd_tot = (sapw_m_turnover + struct_m_turnover) * &
               SF_val_CWD_frac_adj(c) * plant_dens * &
               (1.0_r8-prt_params%allom_agb_frac(pft))

          do ilyr = 1, numlevsoil
             litt%bg_cwd_in(c,ilyr) = litt%bg_cwd_in(c,ilyr) + &
                  bg_cwd_tot * currentSite%rootfrac_scr(ilyr)
          end do

          flux_diags%cwd_bg_input(c)  = flux_diags%cwd_bg_input(c) + &
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


       flux_diags%leaf_litter_input(pft) = &
            flux_diags%leaf_litter_input(pft) +  &
            leaf_m * dead_n*currentPatch%area


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

       flux_diags%root_litter_input(pft) = &
            flux_diags%root_litter_input(pft) +  &
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

          flux_diags%cwd_bg_input(c)  = flux_diags%cwd_bg_input(c) + &
               bg_cwd_tot * currentPatch%area

          ! Send AGB component of boles from logging activities into the litter.
          ! This includes fluxes from indirect modes of death, as well as the
          ! non-exported boles due to direct harvesting.

          if (c==ncwd) then


             trunk_wood =  (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * dead_n_dlogging * &
                  prt_params%allom_agb_frac(pft)

             site_mass%wood_product = site_mass%wood_product + &
                  trunk_wood * currentPatch%area * logging_export_frac

             ! Add AG wood to litter from the non-exported fraction of wood
             ! from direct anthro sources

             litt%ag_cwd_in(c) = litt%ag_cwd_in(c) +  &
                  trunk_wood * (1._r8-logging_export_frac)

             flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
                  trunk_wood * (1._r8-logging_export_frac) * currentPatch%area

             ! Add AG wood to litter from indirect anthro sources

             litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * (dead_n_natural+dead_n_ilogging)  * &
                  prt_params%allom_agb_frac(pft)

             flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
                  SF_val_CWD_frac_adj(c) * (dead_n_natural+dead_n_ilogging) * &
                  currentPatch%area * prt_params%allom_agb_frac(pft)

          else

             litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + (struct_m + sapw_m) * &
                  SF_val_CWD_frac_adj(c) * dead_n  * &
                  prt_params%allom_agb_frac(pft)

             flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
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

    use FatesSynchronizedParamsMod  , only : FatesSynchronizedParamsInst
    use FatesConstantsMod, only : tfrz => t_water_freeze_k_1atm
    use FatesConstantsMod, only : pi => pi_const
    !
    ! !ARGUMENTS
    type(ed_patch_type), intent(inout) :: currentPatch
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
             years_per_day * SF_val_max_decomp(dl_sf) * fragmentation_scaler(soil_layer_index)

       do ilyr = 1,nlev_eff_decomp
           litt%root_fines_frag(dcmpy,ilyr) = litt%root_fines(dcmpy,ilyr) * &
                 years_per_day *  SF_val_max_decomp(dl_sf) * fragmentation_scaler(ilyr)
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
    type(ed_patch_type), pointer :: cpatch
    type(ed_cohort_type), pointer :: ccohort

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
    type(ed_patch_type), pointer :: cpatch
    type(ed_cohort_type), pointer :: ccohort
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
    type(ed_patch_type), pointer :: cpatch
    type(ed_cohort_type), pointer :: ccohort
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
