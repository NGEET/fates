module EDPhysiologyMod

#include "shr_assert.h"

  ! ============================================================================
  ! Miscellaneous physiology routines from ED. 
  ! ============================================================================

  use FatesGlobals, only         : fates_log
  use FatesInterfaceMod, only    : hlm_days_per_year
  use FatesInterfaceMod, only    : hlm_model_day
  use FatesInterfaceMod, only    : hlm_freq_day
  use FatesInterfaceMod, only    : hlm_day_of_year
  use FatesInterfaceMod, only    : numpft
  use FatesInterfaceMod, only    : nleafage
  use FatesInterfaceMod, only    : hlm_use_planthydro
  use FatesInterfaceMod, only    : hlm_parteh_mode
  use FatesConstantsMod, only    : r8 => fates_r8
  use FatesConstantsMod, only    : nearzero
  use FatesConstantsMod, only    : g_per_kg
  use FatesConstantsMod, only    : days_per_sec
  use EDPftvarcon      , only    : EDPftvarcon_inst
  use EDPftvarcon      , only    : GetDecompyFrac
  use FatesInterfaceMod, only    : bc_in_type
  use EDCohortDynamicsMod , only : zero_cohort
  use EDCohortDynamicsMod , only : create_cohort, sort_cohorts
  use EDCohortDynamicsMod , only : InitPRTObject
  use FatesAllometryMod   , only : tree_lai
  use FatesAllometryMod   , only : tree_sai
  use FatesAllometryMod   , only : decay_coeff_kn
  use FatesLitterMod      , only : litter_type
  use EDTypesMod          , only : site_massbal_type
  use EDTypesMod          , only : numlevsoil_max
  use EDTypesMod          , only : numWaterMem
  use EDTypesMod          , only : dl_sf, dinc_ed, area_inv
  use FatesLitterMod      , only : ncwd
  use FatesLitterMod      , only : ndcmpy
  use FatesLitterMod      , only : ilabile
  use FatesLitterMod      , only : ilignin
  use FatesLitterMod      , only : icellulose
  use EDTypesMod          , only : nlevleaf
  use EDTypesMod          , only : num_vegtemp_mem
  use EDTypesMod          , only : maxpft
  use EDTypesMod          , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod          , only : leaves_on
  use EDTypesMod          , only : leaves_off
  use EDTypesMod          , only : min_n_safemath
  use EDTypesMod          , only : num_elements
  use EDTypesMod          , only : element_list
  use EDTypesMod          , only : element_pos
  use EDTypesMod          , only : site_fluxdiags_type
  use EDTypesMod          , only : phen_cstat_nevercold
  use EDTypesMod          , only : phen_cstat_iscold
  use EDTypesMod          , only : phen_cstat_notcold
  use EDTypesMod          , only : phen_dstat_timeoff
  use EDTypesMod          , only : phen_dstat_moistoff
  use EDTypesMod          , only : phen_dstat_moiston
  use EDTypesMod          , only : phen_dstat_timeon

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
  use FatesAllometryMod, only : i_biomass_rootprof_context 
  
  use PRTGenericMod, only : prt_carbon_allom_hyp
  use PRTGenericMod, only : prt_cnp_flex_allom_hyp
  use PRTGenericMod, only : prt_vartypes
  use PRTGenericMod, only : leaf_organ
  use PRTGenericMod, only : sapw_organ, struct_organ
  use PRTGenericMod, only : all_carbon_elements
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
  use PRTLossFluxesMod, only : PRTPhenologyFlush
  use PRTLossFluxesMod, only : PRTDeciduousTurnover
  use PRTLossFluxesMod, only : PRTReproRelease

  implicit none
  private

  public :: trim_canopy
  public :: phenology
  public :: recruitment
  public :: ZeroLitterFluxes
  public :: FluxIntoLitterPools
  public :: ZeroAllocationRates
  public :: PreDisturbanceLitterFluxes
  public :: PreDisturbanceIntegrateLitter  
  public :: SeedIn
  
  logical, parameter :: debug  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
        __FILE__

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
       
       ! Send those decaying seeds in the previous call 
       ! to the litter input flux
       call SeedDecayToFines(litt)
       
       ! Calculate seed germination rate, the status flags prevent
       ! germination from occuring when the site is in a drought 
       ! (for drought deciduous) or too cold (for cold deciduous)
       call SeedGermination(litt, currentSite%cstatus, currentSite%dstatus)
       
       ! Send fluxes from newly created litter into the litter pools
       ! This litter flux is from non-disturbance inducing mortality, as well
       ! as litter fluxes from live trees
       call CWDInput(currentSite, currentPatch, litt)


       ! Only calculate fragmentation flux over layers that are active
       ! (RGK-Mar2019) SHOULD WE MAX THIS AT 1? DONT HAVE TO

       nlev_eff_decomp = max(bc_in%max_rooting_depth_index_col, 1)
       call CWDOut(litt,currentPatch%fragmentation_scaler,nlev_eff_decomp)


       site_mass => currentSite%mass_balance(el)
       
       ! Fragmentation flux to soil decomposition model [kg/site/day]
       site_mass%frag_out = site_mass%frag_out + currentPatch%area * &
            ( sum(litt%ag_cwd_frag) + sum(litt%bg_cwd_frag) + &
            sum(litt%leaf_fines_frag) + sum(litt%root_fines_frag))
       
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

    integer  :: z                ! leaf layer
    integer  :: ipft             ! pft index
    logical  :: trimmed          ! was this layer trimmed in this year? If not expand the canopy. 
    real(r8) :: tar_bl           ! target leaf biomass       (leaves flushed, trimmed)
    real(r8) :: tar_bfr          ! target fine-root biomass  (leaves flushed, trimmed)
    real(r8) :: bfr_per_bleaf    ! ratio of fine root per leaf biomass
    real(r8) :: sla_levleaf      ! sla at leaf level z
    real(r8) :: nscaler_levleaf  ! nscaler value at leaf level z
    integer  :: cl               ! canopy layer index
    real(r8) :: kn               ! nitrogen decay coefficient
    real(r8) :: sla_max          ! Observational constraint on how large sla (m2/gC) can become
    real(r8) :: leaf_c           ! leaf carbon [kg]
    real(r8) :: sapw_c           ! sapwood carbon [kg]
    real(r8) :: store_c          ! storage carbon [kg]
    real(r8) :: struct_c         ! structure carbon [kg]
    real(r8) :: leaf_inc         ! LAI-only portion of the vegetation increment of dinc_ed
    real(r8) :: lai_canopy_above ! the LAI in the canopy layers above the layer of interest
    real(r8) :: lai_layers_above ! the LAI in the leaf layers, within the current canopy, 
                                 ! above the leaf layer of interest
    real(r8) :: lai_current      ! the LAI in the current leaf layer
    real(r8) :: cumulative_lai   ! the cumulative LAI, top down, to the leaf layer of interest

    !----------------------------------------------------------------------

    currentPatch => currentSite%youngest_patch
    do while(associated(currentPatch))
       
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort)) 

          trimmed = .false.
          ipft = currentCohort%pft
          call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread,currentCohort%pft,currentCohort%c_area)

          leaf_c   = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)

          currentCohort%treelai = tree_lai(leaf_c, currentCohort%pft, currentCohort%c_area, &
                                           currentCohort%n, currentCohort%canopy_layer,               &
                                           currentPatch%canopy_layer_tlai,currentCohort%vcmax25top )    

          currentCohort%treesai = tree_sai(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_trim, &
                                           currentCohort%c_area, currentCohort%n, currentCohort%canopy_layer, &
                                           currentPatch%canopy_layer_tlai, currentCohort%treelai, &
                                           currentCohort%vcmax25top,0 )  

          currentCohort%nv      = ceiling((currentCohort%treelai+currentCohort%treesai)/dinc_ed)

          if (currentCohort%nv > nlevleaf)then
             write(fates_log(),*) 'nv > nlevleaf',currentCohort%nv, &
                   currentCohort%treelai,currentCohort%treesai, &
                   currentCohort%c_area,currentCohort%n,leaf_c
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

          call bleaf(currentcohort%dbh,ipft,currentcohort%canopy_trim,tar_bl)

          if ( int(EDPftvarcon_inst%allom_fmode(ipft)) .eq. 1 ) then
             ! only query fine root biomass if using a fine root allometric model that takes leaf trim into account
             call bfineroot(currentcohort%dbh,ipft,currentcohort%canopy_trim,tar_bfr)
             bfr_per_bleaf = tar_bfr/tar_bl
          endif

          ! Identify current canopy layer (cl)
          cl = currentCohort%canopy_layer
          
          ! PFT-level maximum SLA value, even if under a thick canopy (same units as slatop)
          sla_max = EDPftvarcon_inst%slamax(ipft)

          !Leaf cost vs netuptake for each leaf layer. 
          do z = 1, currentCohort%nv

             ! Calculate the cumulative total vegetation area index (no snow occlusion, stems and leaves)

             leaf_inc    = dinc_ed * &
                   currentCohort%treelai/(currentCohort%treelai+currentCohort%treesai)
             
             ! Now calculate the cumulative top-down lai of the current layer's midpoint
             lai_canopy_above  = sum(currentPatch%canopy_layer_tlai(1:cl-1)) 
             lai_layers_above  = leaf_inc * (z-1)
             lai_current       = min(leaf_inc, currentCohort%treelai - lai_layers_above)
             cumulative_lai    = lai_canopy_above + lai_layers_above + 0.5*lai_current
             
             if (currentCohort%year_net_uptake(z) /= 999._r8)then !there was activity this year in this leaf layer.
             
                   
                ! Calculate sla_levleaf following the sla profile with overlying leaf area
                ! Scale for leaf nitrogen profile
                kn = decay_coeff_kn(ipft,currentCohort%vcmax25top)
                ! Nscaler value at leaf level z
                nscaler_levleaf = exp(-kn * cumulative_lai)
                ! Sla value at leaf level z after nitrogen profile scaling (m2/gC)
                sla_levleaf = EDPftvarcon_inst%slatop(ipft)/nscaler_levleaf

                if(sla_levleaf > sla_max)then
                   sla_levleaf = sla_max
                end if
                   
                !Leaf Cost kgC/m2/year-1
                !decidous costs. 
                if (EDPftvarcon_inst%season_decid(ipft) ==  itrue .or. &
                     EDPftvarcon_inst%stress_decid(ipft) == itrue )then 

                   ! Leaf cost at leaf level z accounting for sla profile (kgC/m2)
                   currentCohort%leaf_cost =  1._r8/(sla_levleaf*1000.0_r8)

                   if ( int(EDPftvarcon_inst%allom_fmode(ipft)) .eq. 1 ) then
                      ! if using trimmed leaf for fine root biomass allometry, add the cost of the root increment
                      ! to the leaf increment; otherwise do not.
                      currentCohort%leaf_cost = currentCohort%leaf_cost + &
                           1.0_r8/(sla_levleaf*1000.0_r8) * &
                           bfr_per_bleaf / EDPftvarcon_inst%root_long(ipft)
                   endif

                   currentCohort%leaf_cost = currentCohort%leaf_cost * &
                         (EDPftvarcon_inst%grperc(ipft) + 1._r8)
                else !evergreen costs

                   ! Leaf cost at leaf level z accounting for sla profile
                   currentCohort%leaf_cost = 1.0_r8/(sla_levleaf* &
                        sum(EDPftvarcon_inst%leaf_long(ipft,:))*1000.0_r8) !convert from sla in m2g-1 to m2kg-1
                   
                   
                   if ( int(EDPftvarcon_inst%allom_fmode(ipft)) .eq. 1 ) then
                      ! if using trimmed leaf for fine root biomass allometry, add the cost of the root increment
                      ! to the leaf increment; otherwise do not.
                      currentCohort%leaf_cost = currentCohort%leaf_cost + &
                           1.0_r8/(sla_levleaf*1000.0_r8) * &
                           bfr_per_bleaf / EDPftvarcon_inst%root_long(ipft)
                   endif
                   currentCohort%leaf_cost = currentCohort%leaf_cost * &
                         (EDPftvarcon_inst%grperc(ipft) + 1._r8)
                endif
                if (currentCohort%year_net_uptake(z) < currentCohort%leaf_cost)then
                   if (currentCohort%canopy_trim > EDPftvarcon_inst%trim_limit(ipft))then

                      if ( debug ) then
                         write(fates_log(),*) 'trimming leaves', &
                               currentCohort%canopy_trim,currentCohort%leaf_cost
                      endif

                      ! keep trimming until none of the canopy is in negative carbon balance.              
                      if (currentCohort%hite > EDPftvarcon_inst%hgt_min(ipft))then
                         currentCohort%canopy_trim = currentCohort%canopy_trim - &
                               EDPftvarcon_inst%trim_inc(ipft)
                         if (EDPftvarcon_inst%evergreen(ipft) /= 1)then
                            currentCohort%laimemory = currentCohort%laimemory * &
                                  (1.0_r8 - EDPftvarcon_inst%trim_inc(ipft)) 
                         endif
                         trimmed = .true.
                      endif
                   endif
                endif
             endif !leaf activity? 
          enddo !z

          currentCohort%year_net_uptake(:) = 999.0_r8
          if ( (.not.trimmed) .and.currentCohort%canopy_trim < 1.0_r8)then
             currentCohort%canopy_trim = currentCohort%canopy_trim + EDPftvarcon_inst%trim_inc(ipft)
          endif 

          if ( debug ) then
             write(fates_log(),*) 'trimming',currentCohort%canopy_trim
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
    model_day_int = nint(hlm_model_day)


    ! Use the following layer index to calculate drought conditions
    ilayer_swater = minloc(abs(bc_in%z_sisl(:)-dphen_soil_depth),dim=1)


    ! Parameter of drought decid leaf loss in mm in top layer...FIX(RF,032414) 
    ! - this is arbitrary and poorly understood. Needs work. ED_
    !Parameters: defaults from Botta et al. 2000 GCB,6 709-725 
    !Parameters, default from from SDGVM model of senesence

    temp_in_C = 0._r8
    cpatch => CurrentSite%oldest_patch   
    do while(associated(cpatch))    
       temp_in_C = temp_in_C + bc_in%t_veg24_pa(cpatch%patchno)*cpatch%area
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
       dayssincecleafon = model_day_int - (currentSite%cleafondate-365)
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
        (dayssincecleafoff > 400)) then           ! remove leaves after a whole year 
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
    real(r8) :: totalmemory            ! total memory of carbon [kg]
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

          store_c = currentCohort%prt%GetState(store_organ, all_carbon_elements)
          leaf_c  = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)
          sapw_c  = currentCohort%prt%GetState(sapw_organ, all_carbon_elements)
          struct_c  = currentCohort%prt%GetState(struct_organ, all_carbon_elements)
	  
          stem_drop_fraction = EDPftvarcon_inst%phen_stem_drop_fraction(ipft)
	  
          ! COLD LEAF ON
          ! The site level flags signify that it is no-longer too cold
          ! for leaves. Time to signal flushing

          if (EDPftvarcon_inst%season_decid(ipft) == itrue)then
             if ( currentSite%cstatus == phen_cstat_notcold  )then                ! we have just moved to leaves being on . 
                if (currentCohort%status_coh == leaves_off)then ! Are the leaves currently off?        
                   currentCohort%status_coh = leaves_on         ! Leaves are on, so change status to 
                                                                ! stop flow of carbon out of bstore. 
                   
                   if(store_c>nearzero) then
                   ! flush either the amount required from the laimemory, or -most- of the storage pool
                   ! RF: added a criterion to stop the entire store pool emptying and triggering termination mortality
                   ! n.b. this might not be necessary if we adopted a more gradual approach to leaf flushing... 
                     store_c_transfer_frac =  min((EDPftvarcon_inst%phenflush_fraction(ipft)* &
                     currentCohort%laimemory)/store_c,(1.0_r8-carbon_store_buffer))

                     if(EDPftvarcon_inst%woody(ipft).ne.itrue)then
                        totalmemory=currentCohort%laimemory+currentCohort%sapwmemory+currentCohort%structmemory
                        store_c_transfer_frac = min((EDPftvarcon_inst%phenflush_fraction(ipft)* &
                                                totalmemory)/store_c, (1.0_r8-carbon_store_buffer))
                     endif
		     
                   else
                      store_c_transfer_frac = 0.0_r8
                   end if

                   ! This call will request that storage carbon will be transferred to 
                   ! leaf tissues. It is specified as a fraction of the available storage
                  if(EDPftvarcon_inst%woody(ipft) == itrue) then

                     call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, store_c_transfer_frac)
                     currentCohort%laimemory = 0.0_r8		   

                  else
                  
                     ! Check that the stem drop fraction is set to non-zero amount otherwise flush all carbon store to leaves
                     if (stem_drop_fraction .gt. 0.0_r8) then

                        call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                        store_c_transfer_frac*currentCohort%laimemory/totalmemory)		   

                        call PRTPhenologyFlush(currentCohort%prt, ipft, sapw_organ, &
                        store_c_transfer_frac*currentCohort%sapwmemory/totalmemory)

                        call PRTPhenologyFlush(currentCohort%prt, ipft, struct_organ, & 
                        store_c_transfer_frac*currentCohort%structmemory/totalmemory)

                     else 

                        call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                        store_c_transfer_frac)                        

                     end if
               
                     currentCohort%laimemory = 0.0_r8
                     currentCohort%structmemory = 0.0_r8
                     currentCohort%sapwmemory = 0.0_r8
                  
                  endif		   
                endif !pft phenology
             endif ! growing season 

             !COLD LEAF OFF
             if (currentSite%cstatus == phen_cstat_nevercold .or. &
                 currentSite%cstatus == phen_cstat_iscold) then ! past leaf drop day? Leaves still on tree?  

                if (currentCohort%status_coh == leaves_on) then ! leaves have not dropped

		            ! leaf off occur on individuals bigger than specific size for grass
                  if (currentCohort%dbh > EDPftvarcon_inst%phen_cold_size_threshold(ipft) &
                  .or. EDPftvarcon_inst%woody(ipft)==itrue) then 
                     
                     ! This sets the cohort to the "leaves off" flag
                     currentCohort%status_coh  = leaves_off

                     ! Remember what the lai was (leaf mass actually) was for next year
                     ! the same amount back on in the spring...

                     currentCohort%laimemory   = leaf_c

                     ! Drop Leaves (this routine will update the leaf state variables,
                     ! for carbon and any other element that are prognostic. It will
                     ! also track the turnover masses that will be sent to litter later on)

                     call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                        leaf_organ, leaf_drop_fraction)
         
                     if(EDPftvarcon_inst%woody(ipft).ne.itrue)then
               
                           currentCohort%sapwmemory   = sapw_c * stem_drop_fraction
               
                           currentCohort%structmemory   = struct_c * stem_drop_fraction			 

                           call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                              sapw_organ, stem_drop_fraction)

                           call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                              struct_organ, stem_drop_fraction)

                     endif	! woody plant check
                  endif ! individual dbh size check
               endif !leaf status
            endif !currentSite status
         endif  !season_decid

          ! DROUGHT LEAF ON
          ! Site level flag indicates it is no longer in drought condition
          ! deciduous plants can flush

          if (EDPftvarcon_inst%stress_decid(ipft) == itrue )then
             
            if (currentSite%dstatus == phen_dstat_moiston .or. &
                 currentSite%dstatus == phen_dstat_timeon )then 

               ! we have just moved to leaves being on . 
               if (currentCohort%status_coh == leaves_off)then    

                   !is it the leaf-on day? Are the leaves currently off?    

                  currentCohort%status_coh = leaves_on    ! Leaves are on, so change status to 
                                                           ! stop flow of carbon out of bstore. 

                  if(store_c>nearzero) then

                     store_c_transfer_frac = &
                            min(EDPftvarcon_inst%phenflush_fraction(ipft)*currentCohort%laimemory, store_c)/store_c

                     if(EDPftvarcon_inst%woody(ipft).ne.itrue)then
                     
                        totalmemory=currentCohort%laimemory+currentCohort%sapwmemory+currentCohort%structmemory
                        store_c_transfer_frac = min(EDPftvarcon_inst%phenflush_fraction(ipft)* &
                                                totalmemory, store_c)/store_c

                     endif

                  else
                     store_c_transfer_frac = 0.0_r8
                  endif
                   
                   ! This call will request that storage carbon will be transferred to 
                   ! leaf tissues. It is specified as a fraction of the available storage
                  if(EDPftvarcon_inst%woody(ipft) == itrue) then
                  
                     call PRTPhenologyFlush(currentCohort%prt, ipft, &
                        leaf_organ, store_c_transfer_frac)

                     currentCohort%laimemory = 0.0_r8
               
                  else

                     ! Check that the stem drop fraction is set to non-zero amount otherwise flush all carbon store to leaves
                     if (stem_drop_fraction .gt. 0.0_r8) then

                        call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                        store_c_transfer_frac*currentCohort%laimemory/totalmemory)		   

                        call PRTPhenologyFlush(currentCohort%prt, ipft, sapw_organ, &
                        store_c_transfer_frac*currentCohort%sapwmemory/totalmemory)

                        call PRTPhenologyFlush(currentCohort%prt, ipft, struct_organ, & 
                        store_c_transfer_frac*currentCohort%structmemory/totalmemory)

                     else

                        call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, &
                        store_c_transfer_frac)	

                     end if
                        
                     currentCohort%laimemory = 0.0_r8
                     currentCohort%structmemory = 0.0_r8
                     currentCohort%sapwmemory = 0.0_r8
		     
		            endif ! woody plant check
               endif !currentCohort status again?
            endif   !currentSite status

            !DROUGHT LEAF OFF
            if (currentSite%dstatus == phen_dstat_moistoff .or. &
               currentSite%dstatus == phen_dstat_timeoff) then        

               if (currentCohort%status_coh == leaves_on) then ! leaves have not dropped

                  ! This sets the cohort to the "leaves off" flag
                  currentCohort%status_coh      = leaves_off
                  
                  ! Remember what the lai (leaf mass actually) was for next year
                  currentCohort%laimemory   = leaf_c
      
                  call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                        leaf_organ, leaf_drop_fraction)
			 
                  if(EDPftvarcon_inst%woody(ipft).ne.itrue)then
            
                     currentCohort%sapwmemory   = sapw_c * stem_drop_fraction
                     currentCohort%structmemory   = struct_c * stem_drop_fraction			 

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

       ! If there is forced external seed rain, we calculate the input mass flux
       ! from the different elements, usung the seed optimal stoichiometry
       ! for non-carbon
       select case(element_id)
       case(carbon12_element)
          seed_stoich = 1._r8
       case(nitrogen_element)
          seed_stoich = EDPftvarcon_inst%prt_nitr_stoich_p2(pft,repro_organ)
       case(phosphorus_element)
          seed_stoich = EDPftvarcon_inst%prt_phos_stoich_p2(pft,repro_organ)
       case default
          write(fates_log(), *) 'undefined element specified'
          write(fates_log(), *) 'while defining forced external seed mass flux'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end select
          

       ! Loop over all patches and sum up the seed input for each PFT
       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))

          litt => currentPatch%litter(el)
          do pft = 1,numpft
             
             ! Seed input from local sources (within site)
             litt%seed_in_local(pft) = litt%seed_in_local(pft) + site_seed_rain(pft)/area
             
             ! Seed input from external sources (user param seed rain, or dispersal model)
             seed_in_external =  seed_stoich*EDPftvarcon_inst%seed_suppl(pft)*years_per_day
             
             litt%seed_in_extern(pft) = litt%seed_in_extern(pft) + seed_in_external

             ! Seeds entering externally [kg/site/day]
             site_mass%seed_in = site_mass%seed_in + seed_in_external*currentPatch%area

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

       if ((EDPftvarcon_inst%season_decid(pft) == itrue ) .and. &
             (any(cold_stat == [phen_cstat_nevercold,phen_cstat_iscold]))) then
           litt%seed_germ_in(pft) = 0.0_r8
       endif
       if ((EDPftvarcon_inst%stress_decid(pft) == itrue ) .and. &
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
    use FatesInterfaceMod, only : hlm_use_ed_prescribed_phys
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target   :: currentSite
    type(ed_patch_type), intent(inout), pointer :: currentPatch
    type(bc_in_type), intent(in)                :: bc_in
    !
    ! !LOCAL VARIABLES:
    class(prt_vartypes), pointer :: prt
    integer :: ft
    type (ed_cohort_type) , pointer :: temp_cohort
    type (litter_type), pointer     :: litt          ! The litter object (carbon right now)
    type(site_massbal_type), pointer :: site_mass    ! For accounting total in-out mass fluxes
    integer :: cohortstatus
    integer :: el          ! loop counter for element
    integer :: element_id  ! element index consistent with definitions in PRTGenericMod
    integer :: iage        ! age loop counter for leaf age bins
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

       temp_cohort%canopy_trim = 0.8_r8  !starting with the canopy not fully expanded 
       temp_cohort%pft         = ft
       temp_cohort%hite        = EDPftvarcon_inst%hgt_min(ft)
       temp_cohort%coage       = 0.0_r8
       stem_drop_fraction = EDPftvarcon_inst%phen_stem_drop_fraction(ft)

       call h2d_allom(temp_cohort%hite,ft,temp_cohort%dbh)

       ! Initialize live pools
       call bleaf(temp_cohort%dbh,ft,temp_cohort%canopy_trim,c_leaf)
       call bfineroot(temp_cohort%dbh,ft,temp_cohort%canopy_trim,c_fnrt)
       call bsap_allom(temp_cohort%dbh,ft,temp_cohort%canopy_trim,a_sapw, c_sapw)
       call bagw_allom(temp_cohort%dbh,ft,c_agw)
       call bbgw_allom(temp_cohort%dbh,ft,c_bgw)
       call bdead_allom(c_agw,c_bgw,c_sapw,ft,c_struct)
       call bstore_allom(temp_cohort%dbh,ft,temp_cohort%canopy_trim,c_store)

       ! Default assumption is that leaves are on
       cohortstatus = leaves_on
       temp_cohort%laimemory = 0.0_r8     
       temp_cohort%sapwmemory = 0.0_r8     
       temp_cohort%structmemory = 0.0_r8     

       
       ! But if the plant is seasonally (cold) deciduous, and the site status is flagged
       ! as "cold", then set the cohort's status to leaves_off, and remember the leaf biomass
       if ((EDPftvarcon_inst%season_decid(ft) == itrue) .and. &
             (any(currentSite%cstatus == [phen_cstat_nevercold,phen_cstat_iscold]))) then
         temp_cohort%laimemory = c_leaf
         c_leaf = 0.0_r8

         ! If plant is not woody then set sapwood and structural biomass as well
         if (EDPftvarcon_inst%woody(ft).ne.itrue) then
            temp_cohort%sapwmemory = c_sapw * stem_drop_fraction
            temp_cohort%structmemory = c_struct * stem_drop_fraction
            c_sapw = (1.0_r8 - stem_drop_fraction) * c_sapw 
            c_struct = (1.0_r8 - stem_drop_fraction) * c_struct
         endif
         cohortstatus = leaves_off
       endif
       
       ! Or.. if the plant is drought deciduous, and the site status is flagged as 
       ! "in a drought", then likewise, set the cohort's status to leaves_off, and remember leaf
       ! biomass
       if ((EDPftvarcon_inst%stress_decid(ft) == itrue) .and. &
             (any(currentSite%dstatus == [phen_dstat_timeoff,phen_dstat_moistoff]))) then
         temp_cohort%laimemory = c_leaf
         c_leaf = 0.0_r8

         ! If plant is not woody then set sapwood and structural biomass as well
         if(EDPftvarcon_inst%woody(ft).ne.itrue)then
            temp_cohort%sapwmemory = c_sapw * stem_drop_fraction
            temp_cohort%structmemory = c_struct * stem_drop_fraction
            c_sapw = (1.0_r8 - stem_drop_fraction) * c_sapw 
            c_struct = (1.0_r8 - stem_drop_fraction) * c_struct
         endif
         cohortstatus = leaves_off
       endif
      

       ! Cycle through available carbon and nutrients, find the limiting element
       ! to dictate the total number of plants that can be generated

       if ( (hlm_use_ed_prescribed_phys .eq. ifalse) .or. &
            (EDPftvarcon_inst%prescribed_recruitment(ft) .lt. 0._r8) ) then

           temp_cohort%n = 1.e10_r8

           do el = 1,num_elements
               
               element_id = element_list(el)
               select case(element_id)
               case(carbon12_element)
               
                  mass_demand = (c_struct+c_leaf+c_fnrt+c_sapw+c_store)
               
               case(nitrogen_element)
               
                  mass_demand = c_struct*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,struct_organ) + &
                                 c_leaf*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,leaf_organ) + &
                                 c_fnrt*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,fnrt_organ) + & 
                                 c_sapw*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,sapw_organ) + & 
                                 c_store*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,store_organ)
               
               case(phosphorus_element)
               
                  mass_demand = c_struct*EDPftvarcon_inst%prt_phos_stoich_p1(ft,struct_organ) + &
                                 c_leaf*EDPftvarcon_inst%prt_phos_stoich_p1(ft,leaf_organ) + &
                                 c_fnrt*EDPftvarcon_inst%prt_phos_stoich_p1(ft,fnrt_organ) + & 
                                 c_sapw*EDPftvarcon_inst%prt_phos_stoich_p1(ft,sapw_organ)  + & 
                                 c_store*EDPftvarcon_inst%prt_phos_stoich_p1(ft,store_organ)
               
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
       if (temp_cohort%n > min_n_safemath )then

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

                 m_struct = c_struct*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,struct_organ)
                 m_leaf   = c_leaf*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,leaf_organ)
                 m_fnrt   = c_fnrt*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,fnrt_organ)
                 m_sapw   = c_sapw*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,sapw_organ)
                 m_store  = c_store*EDPftvarcon_inst%prt_nitr_stoich_p1(ft,store_organ)
                 m_repro  = 0._r8

              case(phosphorus_element)

                 m_struct = c_struct*EDPftvarcon_inst%prt_phos_stoich_p1(ft,struct_organ)
                 m_leaf   = c_leaf*EDPftvarcon_inst%prt_phos_stoich_p1(ft,leaf_organ)
                 m_fnrt   = c_fnrt*EDPftvarcon_inst%prt_phos_stoich_p1(ft,fnrt_organ)
                 m_sapw   = c_sapw*EDPftvarcon_inst%prt_phos_stoich_p1(ft,sapw_organ)
                 m_store  = c_store*EDPftvarcon_inst%prt_phos_stoich_p1(ft,store_organ)
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
                temp_cohort%laimemory, temp_cohort%sapwmemory, temp_cohort%structmemory, &
                cohortstatus, recruitstatus, &
                temp_cohort%canopy_trim, currentPatch%NCL_p, currentSite%spread, bc_in)
           
           ! Note that if hydraulics is on, the number of cohorts may had
           ! changed due to hydraulic constraints.
           ! This constaint is applied during "create_cohort" subroutine.
           
           ! keep track of how many individuals were recruited for passing to history
           currentSite%recruitment_rate(ft) = currentSite%recruitment_rate(ft) + temp_cohort%n
          

       endif
     enddo  !pft loop
     
     deallocate(temp_cohort) ! delete temporary cohort

  end subroutine recruitment

  ! ============================================================================

  subroutine CWDInput( currentSite, currentPatch, litt)

    !
    ! !DESCRIPTION:
    ! Generate litter fields from turnover.
    ! Note, that the when this is called, the number density of the plants
    ! has not been reduced from non-mortal turnover yet.
    ! Thus, we need to avoid double counting losses from dying trees
    ! and turnover in dying trees.
    !
    ! !USES:
    use SFParamsMod , only : SF_val_CWD_frac

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target :: currentSite
    type(ed_patch_type),intent(inout), target :: currentPatch
    type(litter_type),intent(inout),target    :: litt


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
            icontext = i_biomass_rootprof_context)

      leaf_m_turnover   = currentCohort%prt%GetTurnover(leaf_organ,element_id)
      store_m_turnover  = currentCohort%prt%GetTurnover(store_organ,element_id)
      fnrt_m_turnover   = currentCohort%prt%GetTurnover(fnrt_organ,element_id)
      sapw_m_turnover   = currentCohort%prt%GetTurnover(sapw_organ,element_id)
      struct_m_turnover = currentCohort%prt%GetTurnover(struct_organ,element_id)
      repro_m_turnover  = currentCohort%prt%GetTurnover(repro_organ,element_id)

      leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id)
      store_m         = currentCohort%prt%GetState(store_organ,element_id)
      fnrt_m          = currentCohort%prt%GetState(fnrt_organ,element_id)
      sapw_m          = currentCohort%prt%GetState(sapw_organ,element_id)
      struct_m        = currentCohort%prt%GetState(struct_organ,element_id)
      repro_m         = currentCohort%prt%GetState(repro_organ,element_id)

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
      
      do c = 1,ncwd
         litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + &
              (sapw_m_turnover + struct_m_turnover) * &
              SF_val_CWD_frac(c) * plant_dens * &
              EDPftvarcon_inst%allom_agb_frac(pft)

         flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
               (struct_m_turnover + sapw_m_turnover) * SF_val_CWD_frac(c) * &
               EDPftvarcon_inst%allom_agb_frac(pft) * currentCohort%n

         bg_cwd_tot = (sapw_m_turnover + struct_m_turnover) * &
              SF_val_CWD_frac(c) * plant_dens * & 
              (1.0_r8-EDPftvarcon_inst%allom_agb_frac(pft))

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
              SF_val_CWD_frac(c) * dead_n * &
              (1.0_r8-EDPftvarcon_inst%allom_agb_frac(pft))
         
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
                 SF_val_CWD_frac(c) * dead_n_dlogging * &
                 EDPftvarcon_inst%allom_agb_frac(pft) 
            
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
                 SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging)  * &
                 EDPftvarcon_inst%allom_agb_frac(pft)

            flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
                  SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging) * &
                  currentPatch%area * EDPftvarcon_inst%allom_agb_frac(pft)

         else

            litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + (struct_m + sapw_m) * & 
                 SF_val_CWD_frac(c) * dead_n  * &
                 EDPftvarcon_inst%allom_agb_frac(pft)

            flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
                  SF_val_CWD_frac(c) * dead_n * (struct_m + sapw_m) * &
                  currentPatch%area * EDPftvarcon_inst%allom_agb_frac(pft)
            
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
                  SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging) * & 
                  currentPatch%area
            
            currentSite%resources_management%delta_biomass_stock = &
                  currentSite%resources_management%delta_biomass_stock + &
                  (struct_m + sapw_m) * &
                  SF_val_CWD_frac(c) * dead_n * currentPatch%area
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

  subroutine SeedDecayToFines(litt)

    type(litter_type) :: litt
    !
    ! !LOCAL VARIABLES:
    integer  ::  pft

    ! Add decaying seeds to the leaf litter
    ! -----------------------------------------------------------------------------------

    do pft = 1,numpft

        litt%leaf_fines_in(ilabile) = litt%leaf_fines_in(ilabile) + & 
              (litt%seed_decay(pft) + litt%seed_germ_decay(pft)) * EDPftvarcon_inst%lf_flab(pft)
        
        litt%leaf_fines_in(icellulose) = litt%leaf_fines_in(icellulose) + & 
              (litt%seed_decay(pft) + litt%seed_germ_decay(pft)) * EDPftvarcon_inst%lf_fcel(pft)
        
        litt%leaf_fines_in(ilignin) = litt%leaf_fines_in(ilignin) + & 
              (litt%seed_decay(pft) + litt%seed_germ_decay(pft)) * EDPftvarcon_inst%lf_flig(pft)

    enddo
    
    
    return
  end subroutine SeedDecayToFines
  
  



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
    integer  :: j
    integer  :: ifp                   ! Index of a FATES Patch "ifp"
    real(r8) :: t_scalar
    real(r8) :: w_scalar
    real(r8) :: catanf                ! hyperbolic temperature function from CENTURY
    real(r8) :: catanf_30             ! hyperbolic temperature function from CENTURY
    real(r8) :: t1                    ! temperature argument
    !----------------------------------------------------------------------

    catanf(t1) = 11.75_r8 +(29.7_r8 / pi) * atan( pi * 0.031_r8  * ( t1 - 15.4_r8 ))
    catanf_30 = catanf(30._r8)
    
    ifp = currentPatch%patchno 

    if ( .not. use_century_tfunc ) then
    !calculate rate constant scalar for soil temperature,assuming that the base rate constants 
    !are assigned for non-moisture limiting conditions at 25C.
      if (bc_in%t_veg24_pa(ifp)  >=  tfrz) then
        t_scalar = q10_mr**((bc_in%t_veg24_pa(ifp)-(tfrz+25._r8))/10._r8)
                 !  Q10**((t_soisno(c,j)-(tfrz+25._r8))/10._r8)
      else
        t_scalar = (q10_mr**(-25._r8/10._r8))*(q10_froz**((bc_in%t_veg24_pa(ifp)-tfrz)/10._r8))
                  !Q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-tfrz)/10._r8)
      endif
    else
      ! original century uses an arctangent function to calculate the 
      ! temperature dependence of decomposition      
      t_scalar = max(catanf(bc_in%t_veg24_pa(ifp)-tfrz)/catanf_30,0.01_r8)
    endif    
   
    !Moisture Limitations   
    !BTRAN APPROACH - is quite simple, but max's out decomp at all unstressed 
    !soil moisture values, which is not realistic.  
    !litter decomp is proportional to water limitation on average... 
    w_scalar = sum(currentPatch%btran_ft(1:numpft))/real(numpft,r8)

    currentPatch%fragmentation_scaler =  min(1.0_r8,max(0.0_r8,t_scalar * w_scalar))
    
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
    
    real(r8),intent(in)                        :: fragmentation_scaler

    ! This is not necessarily every soil layer, this is the number
    ! of effective layers that are active and can be sent
    ! to the soil decomposition model
    integer,intent(in)                         :: nlev_eff_decomp  
    
    !
    ! !LOCAL VARIABLES:
    integer :: c
    integer :: ilyr
    integer :: dcmpy
    integer :: numlevsoil
    !----------------------------------------------------------------------

    do c = 1,ncwd  

       litt%ag_cwd_frag(c)   = litt%ag_cwd(c) * SF_val_max_decomp(c) * &
             years_per_day * fragmentation_scaler
       
       do ilyr = 1,nlev_eff_decomp
           
           litt%bg_cwd_frag(c,ilyr) = litt%bg_cwd(c,ilyr) * SF_val_max_decomp(c) * &
                years_per_day * fragmentation_scaler

       enddo
    end do

    ! this is the rate at which dropped leaves stop being part of the burnable pool 
    ! and begin to be part of the decomposing pool. This should probably be highly 
    ! sensitive to moisture, but also to the type of leaf thick leaves can dry out 
    ! before they are decomposed, for example. This section needs further scientific input. 

    do dcmpy = 1,ndcmpy

       litt%leaf_fines_frag(dcmpy) = litt%leaf_fines(dcmpy) * &
             years_per_day * SF_val_max_decomp(dl_sf) * fragmentation_scaler
       
       do ilyr = 1,nlev_eff_decomp
           litt%root_fines_frag(dcmpy,ilyr) = litt%root_fines(dcmpy,ilyr) * &
                 years_per_day *  SF_val_max_decomp(dl_sf) * fragmentation_scaler
       end do
    enddo

  end subroutine CWDOut

  ! =====================================================================================

  subroutine FluxIntoLitterPools(nsites, sites, bc_in, bc_out)
    
    ! -----------------------------------------------------------------------------------
    ! Created by Charlie Koven and Rosie Fisher, 2014-2015
    ! take the flux out of the fragmenting litter pools and port into the decomposing 
    ! litter pools. 
    ! in this implementation, decomposing pools are assumed to be humus and non-flammable, 
    ! whereas fragmenting pools are assumed to be physically fragmenting but not 
    ! respiring. This is a simplification, but allows us to 
    !
    ! a) reconcile the need to track both chemical fractions (lignin, cellulose, labile) 
    !    and size fractions (trunk, branch, etc.)
    ! b) to impose a realistic delay on the surge of nutrients into the litter pools 
    !    when large CWD is added to the system via mortality
    !
    ! Because of the different subgrid structure, this subroutine includes the functionality
    ! that in the big-leaf BGC model, is calculated in SoilBiogeochemVerticalProfileMod
    !
    ! The ED code is resolved at a daily timestep, but all of the CN-BGC fluxes are passed 
    ! in as derivatives per second, and then accumulated in the CNStateUpdate routines.  
    ! One way of doing this is to pass back the CN fluxes per second, and keep them 
    ! constant for the whole day (making sure they are not overwritten.  This means that 
    ! the carbon gets passed back and forth between the photosynthesis code 
    ! (fast timestepping) to the ED code (slow timestepping), back to the BGC code 
    ! (fast timestepping).  This means that the state update for the litter pools and 
    ! for the CWD pools occurs at different timescales. 
    ! -----------------------------------------------------------------------------------

    use EDTypesMod, only : AREA
    use FatesConstantsMod, only : sec_per_day
    use FatesInterfaceMod, only : bc_in_type, bc_out_type
    use FatesInterfaceMod, only : hlm_use_vertsoilc
    use FatesInterfaceMod, only : hlm_numlevgrnd
    use FatesConstantsMod, only : itrue
    use FatesGlobals, only : endrun => fates_endrun
    use EDParamsMod , only : ED_val_cwd_flig, ED_val_cwd_fcel
   
    

    implicit none   

    ! !ARGUMENTS    
    integer            , intent(in)            :: nsites
    type(ed_site_type) , intent(inout)         :: sites(nsites)
    type(bc_in_type)   , intent(in)            :: bc_in(:)
    type(bc_out_type)  , intent(inout), target :: bc_out(:)

    ! !LOCAL VARIABLES:
    type (ed_patch_type),  pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    real(r8), pointer              :: flux_cel_si(:)
    real(r8), pointer              :: flux_lab_si(:)
    real(r8), pointer              :: flux_lig_si(:)
    type(litter_type), pointer     :: litt
     
    real(r8) :: surface_prof(1:hlm_numlevgrnd) ! this array is used to distribute
                                               ! fragmented litter on the surface
                                               ! into the soil/decomposition
                                               ! layers. It exponentially decays
    real(r8) :: surface_prof_tot ! normalizes the surface_prof array
    integer  :: ft               ! PFT number
    integer  :: nlev_eff_soil    ! number of effective soil layers
    integer  :: nlev_eff_decomp  ! number of effective decomp layers
    real(r8) :: area_frac        ! fraction of site's area of current patch
    real(r8) :: z_decomp         ! Used for calculating depth midpoints of decomp layers
    integer  :: s                ! Site index
    integer  :: el               ! Element index (C,N,P,etc)
    integer  :: j                ! Soil layer index
    integer  :: id               ! Decomposition layer index
    integer  :: ic               ! CWD type index

    ! NOTE(rgk, 201705) this parameter was brought over from SoilBiogeochemVerticalProfile
    ! how steep profile is for surface components (1/ e_folding depth) (1/m) 
    real(r8),  parameter :: surfprof_exp  = 10.

    do s = 1,nsites

       ! This is the number of effective soil layers to transfer from
       nlev_eff_soil   = max(bc_in(s)%max_rooting_depth_index_col, 1)

       ! The decomposition layers are most likely the exact same layers
       ! as the soil layers (same depths also), unless it is a simplified
       ! single layer case, where nlevdecomp = 1

       nlev_eff_decomp = min(bc_in(s)%nlevdecomp,nlev_eff_soil)

       ! define a single shallow surface profile for surface additions 
       ! (leaves, stems, and N deposition). This sends the above ground
       ! mass into the soil pools using an exponential depth decay function.
       ! Since it is sending an absolute mass [kg] into variable layer
       ! widths, we multiply the profile by the layer width, so that
       ! wider layers get proportionally more.  After the masses
       ! are sent, each layer will normalize by depth.
       
       surface_prof(:) = 0._r8
       z_decomp = 0._r8
       do id = 1,nlev_eff_decomp
          z_decomp = z_decomp+0.5*bc_in(s)%dz_decomp_sisl(id)
          surface_prof(id) = exp(-surfprof_exp * z_decomp) *  bc_in(s)%dz_decomp_sisl(id)
          z_decomp = z_decomp+0.5*bc_in(s)%dz_decomp_sisl(id)
       end do
       surface_prof_tot = sum(surface_prof)
       do id = 1,nlev_eff_decomp
          surface_prof(id) = surface_prof(id)/surface_prof_tot
       end do

       ! Loop over the different elements. 
       do el = 1, num_elements

          ! Zero out the boundary flux arrays
          ! Make a pointer to the cellulose, labile and lignan
          ! flux partitions.

          select case (element_list(el))
          case (carbon12_element)
             bc_out(s)%litt_flux_cel_c_si(:) = 0._r8
             bc_out(s)%litt_flux_lig_c_si(:) = 0._r8
             bc_out(s)%litt_flux_lab_c_si(:) = 0._r8
             flux_cel_si => bc_out(s)%litt_flux_cel_c_si(:)
             flux_lab_si => bc_out(s)%litt_flux_lab_c_si(:)
             flux_lig_si => bc_out(s)%litt_flux_lig_c_si(:)
          case (nitrogen_element) 
             bc_out(s)%litt_flux_cel_n_si(:) = 0._r8
             bc_out(s)%litt_flux_lig_n_si(:) = 0._r8
             bc_out(s)%litt_flux_lab_n_si(:) = 0._r8
             flux_cel_si => bc_out(s)%litt_flux_cel_n_si(:)
             flux_lab_si => bc_out(s)%litt_flux_lab_n_si(:)
             flux_lig_si => bc_out(s)%litt_flux_lig_n_si(:)
          case (phosphorus_element)
             bc_out(s)%litt_flux_cel_p_si(:) = 0._r8
             bc_out(s)%litt_flux_lig_p_si(:) = 0._r8
             bc_out(s)%litt_flux_lab_p_si(:) = 0._r8
             flux_cel_si => bc_out(s)%litt_flux_cel_p_si(:)
             flux_lab_si => bc_out(s)%litt_flux_lab_p_si(:)
             flux_lig_si => bc_out(s)%litt_flux_lig_p_si(:)
          end select
          
          currentPatch => sites(s)%oldest_patch
          do while (associated(currentPatch))
             
             ! Set a pointer to the litter object
             ! for the current element on the current
             ! patch
             litt       => currentPatch%litter(el)
             area_frac  = currentPatch%area/area
             
             do ic = 1, ncwd
                
                do id = 1,nlev_eff_decomp
                   flux_cel_si(id) = flux_cel_si(id) + &
                         litt%ag_cwd_frag(ic) * ED_val_cwd_fcel * area_frac * surface_prof(id)
                   
                   flux_lig_si(id) = flux_lig_si(id) + & 
                         litt%ag_cwd_frag(ic) * ED_val_cwd_flig * area_frac * surface_prof(id)
                end do
                   
                do j = 1, nlev_eff_soil
                   
                   id = bc_in(s)%decomp_id(j)  ! Map from soil layer to decomp layer
                   
                   flux_cel_si(id) = flux_cel_si(id) + &
                         litt%bg_cwd_frag(ic,j) * ED_val_cwd_fcel * area_frac
                   
                   flux_lig_si(id) = flux_lig_si(id) + &
                         litt%bg_cwd_frag(ic,j) * ED_val_cwd_flig * area_frac
                   
                end do
             end do
             
             ! leaf and fine root fragmentation fluxes
                
             do id = 1,nlev_eff_decomp
                   
                 flux_lab_si(id) = flux_lab_si(id) + &
                       litt%leaf_fines_frag(ilabile) * area_frac* surface_prof(id)
                   
                 flux_cel_si(id) = flux_cel_si(id) + &
                       litt%leaf_fines_frag(icellulose) * area_frac* surface_prof(id)
                 
                 flux_lig_si(id) = flux_lig_si(id) + &
                       litt%leaf_fines_frag(ilignin) * area_frac* surface_prof(id)
                 
             end do

             do j = 1, nlev_eff_soil
                 
                 id = bc_in(s)%decomp_id(j)

                 flux_lab_si(id) = flux_lab_si(id) + &
                       litt%root_fines_frag(ilabile,j) * area_frac
                 flux_cel_si(id) = flux_cel_si(id) + &
                       litt%root_fines_frag(icellulose,j) * area_frac
                 flux_lig_si(id) = flux_lig_si(id) + &
                       litt%root_fines_frag(ilignin,j) * area_frac
             enddo

         
             currentPatch => currentPatch%younger
          end do
          
          ! Normalize all masses over the decomposition layer's depth
          ! Convert from kg/m2/day -> g/m3/s

          do id = 1,nlev_eff_decomp
             flux_cel_si(id) = days_per_sec * g_per_kg * &
                               flux_cel_si(id) / bc_in(s)%dz_decomp_sisl(id)
             flux_lig_si(id) = days_per_sec * g_per_kg * &
                               flux_lig_si(id) / bc_in(s)%dz_decomp_sisl(id)
             flux_lab_si(id) = days_per_sec * g_per_kg * &
                               flux_lab_si(id) / bc_in(s)%dz_decomp_sisl(id)
          end do

       end do  ! do elements
       
    end do  ! do sites(s)
    return
end subroutine FluxIntoLitterPools



end module EDPhysiologyMod
