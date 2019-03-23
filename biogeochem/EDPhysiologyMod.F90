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
  use FatesInterfaceMod, only    : bc_in_type
  use EDCohortDynamicsMod , only : zero_cohort
  use EDCohortDynamicsMod , only : create_cohort, sort_cohorts,InitPRTCohort
  use FatesAllometryMod   , only : tree_lai
  use FatesAllometryMod   , only : tree_sai
  use FatesAllometryMod   , only : decay_coeff_kn

  use EDTypesMod          , only : numlevsoil_max
  use EDTypesMod          , only : numWaterMem
  use EDTypesMod          , only : dl_sf, dinc_ed
  use EDTypesMod          , only : external_recruitment
  use EDTypesMod          , only : ncwd
  use EDTypesMod          , only : nlevleaf
  use EDTypesMod          , only : senes
  use EDTypesMod          , only : maxpft
  use EDTypesMod          , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod          , only : dump_cohort
  use EDTypesMod          , only : first_leaf_aclass
  use EDTypesMod          , only : leaves_on
  use EDTypesMod          , only : leaves_off
  use EDTypesMod          , only : min_n_safemath
  use EDTypesMod          , only : num_elements
  use EDTypesMod          , only : element_list
  use EDTypesMod          , only : element_pos
  use EDTypesMod          , only : site_fluxdiags_type

  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use FatesGlobals          , only : fates_log
  use FatesGlobals          , only : endrun => fates_endrun
  use EDParamsMod           , only : fates_mortality_disturbance_fraction
  use FatesPlantHydraulicsMod  , only : AccumulateMortalityWaterStorage
  
  use FatesConstantsMod     , only : itrue,ifalse
  use FatesConstantsMod     , only : calloc_abs_error
  use FatesConstantsMod     , only : days_per_year

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
  use FatesAllometryMod  , only : StructureResetOfDH
  
  use PRTGenericMod, only : prt_carbon_allom_hyp
  use PRTGenericMod, only : leaf_organ
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

  implicit none
  private

  public :: non_canopy_derivs
  public :: trim_canopy
  public :: phenology
  private :: phenology_leafonoff
  public :: recruitment
  private :: cwd_input
  private :: cwd_out
  private :: fragmentation_scaler
  private :: SeedsIn
  private :: SeedDecay
  private :: SeedGermination
  public :: ZeroLitterFluxes
  public :: flux_into_litter_pools
  public :: ZeroAllocationRates
  

  logical, parameter :: debug  = .false. ! local debug flag
  character(len=*), parameter, private :: sourcefile = &
        __FILE__

  integer, parameter :: i_dbh  = 1    ! Array index associated with dbh
  integer, parameter :: i_cleaf = 2   ! Array index associated with leaf carbon
  integer, parameter :: i_cfroot = 3  ! Array index associated with fine-root carbon
  integer, parameter :: i_csap   = 4  ! Array index associated with sapwood carbon
  integer, parameter :: i_cstore = 5  ! Array index associated with storage carbon
  integer, parameter :: i_cdead = 6   ! Array index associated with structural carbon
  integer, parameter :: i_crepro = 7  ! Array index associated with reproductive carbon 
  integer, parameter :: n_cplantpools = 7 ! Size of the carbon only integration framework

  ! ============================================================================

contains

  subroutine ZeroLitterFluxes( currentSite )

    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), pointer               :: currentPatch
    
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
    ! -----------------------------------------------------------------------------------
     

    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout)         :: currentPatch
    type(bc_in_type), intent(in)               :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(litt_vartype), pointer :: litt    ! Points to the litter object for 
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
       
       ! Calculate seed germination rate
       call SeedGermination(litt, currentSite%is_cold, currentSite%is_drought)
       
       ! Send fluxes from newly created litter into the litter pools
       ! This litter flux is from non-disturbance inducing mortality, as well
       ! as litter fluxes from live trees
       call CWDInput( currentSite, currentPatch, litt, bc_in)


       ! Only calculate fragmentation flux over layers that are active
       ! (RGK-Mar2019) SHOULD WE MAX THIS AT 1? DONT HAVE TO

       nlev_eff_decomp = max(bc_in(s)%max_rooting_depth_index_col, 1)
       call CWDOut(litt,currentPatch%fragmentation_scaler,nlev_eff_decomp)


       site_mass => currentSite%mass_balance(el)
       
       ! Seeds entering externally [kg/site]
       site_mass%seed_in = site_mass%seed_in + sum(litt%seed_in_extern(:))*currentPatch%area
       
       ! Fragmentation flux to soil decomposition model [kg/site]
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
    type(litt_vartype), pointer :: litt 
    integer :: el          ! Loop counter for litter element type
    integer :: pft         ! pft loop counter
    integer :: c           ! CWD loop counter
    integer :: nlevsoil    ! number of soil layers
    integer :: ilyr        ! soil layer loop counter

    do el = 1, num_elements
       
       litt => currentPatch%litter(el)
       nlevsoil = size(litt%bg_cwd,dim=2)
       
       ! Update the bank of viable seeds
       ! -----------------------------------------------------------------------------------
       
       do pft = 1,numpft
          litt%seed(pft) = litt%seed(pft) &
               + litt%seed_in_local(pft)  &
               + litt%seed_in_extern(pft) &
               - litt%seed_decay(pft)     &
               - litt%seed_germination(pft)
       enddo
       
       ! Update the Coarse Woody Debris pools (above and below)
       ! -----------------------------------------------------------------------------------
       
       do c = 1,ncwd
          litt%ag_cwd(c) = litt%ag_cwd(c)  + litt%ag_cwd_in(c) - litt%ag_cwd_frag(c)
          do ilyr=1,nlevsoil
             litt%bg_cwd(c,ilyr) =  currentPatch%bg_cwd(c,ilyr) &
                  + litt%bg_cwd_in(c,ilyr) &
                  - litt%bg_cwd_frac(c,ilyr)
          enddo
       end do
    
       ! Update the fine litter pools from leaves and fine-roots
       ! -----------------------------------------------------------------------------------
       
       do pft = 1,numpft
          litt%leaf_fines(pft) = litt%leaf_fines(pft) &
               + litt%leaf_fines_in(pft)              &
               - litt%leaf_fines_frag(pft)
          do ilyr=1,nlevsoil
             litt%root_fines(pft,ilyr) = litt%root_fines(pft,ilyr) &
                  + litt%root_fines_in(pft,ilyr)      &
                  - litt%root_fines_frag(pft,ilyr)
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
    real(r8) :: fnrt_c           ! fineroot carbon [kg]
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
                                           currentPatch%canopy_layer_tlai, currentCohort%treelai,currentCohort%vcmax25top )  

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
                if (EDPftvarcon_inst%season_decid(ipft) == 1.or. &
                     EDPftvarcon_inst%stress_decid(ipft) == 1)then 

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

    integer  :: t            ! day of year
    integer  :: ncolddays    ! no days underneath the threshold for leaf drop
    integer  :: i
    integer  :: timesincedleafon,timesincedleafoff,timesinceleafon,timesinceleafoff
    integer  :: refdate
    integer  :: curdate
    
    integer  :: yr                       ! year (0, ...)
    integer  :: mon                      ! month (1, ..., 12)
    integer  :: day                      ! day of month (1, ..., 31)
    integer  :: sec                      ! seconds of the day

    real(r8) :: leaf_c                       ! leaf carbon [kg]
    real(r8) :: fnrt_c               ! fineroot carbon [kg]
    real(r8) :: sapw_c               ! sapwood carbon [kg]
    real(r8) :: store_c              ! storage carbon [kg]
    real(r8) :: struct_c             ! structure carbon [kg]

    real(r8) :: gdd_threshold
    integer  :: ncdstart     ! beginning of counting period for chilling degree days.
    integer  :: gddstart     ! beginning of counting period for growing degree days.
    real(r8) :: temp_in_C    ! daily averaged temperature in celcius

    real(r8), parameter :: canopy_leaf_lifespan = 365.0_r8    ! Mean lifespan canopy leaves
                                                              ! FIX(RGK 07/10/17)
                                                              ! This is a band-aid on unusual code
                                                              

    ! Parameter of drought decid leaf loss in mm in top layer...FIX(RF,032414) 
    ! - this is arbitrary and poorly understood. Needs work. ED_

    !Parameters: defaults from Botta et al. 2000 GCB,6 709-725 
    !Parameters, default from from SDGVM model of senesence

    t  = hlm_day_of_year
    temp_in_C = bc_in%t_veg24_si - tfrz

    !-----------------Cold Phenology--------------------!              

    !Zero growing degree and chilling day counters
    if (currentSite%lat > 0)then
       ncdstart = 270  !Northern Hemisphere begining November
       gddstart = 1    !Northern Hemisphere begining January
    else
       ncdstart = 120  !Southern Hemisphere beginning May
       gddstart = 181  !Northern Hemisphere begining July
    endif
    
    ! FIX(SPM,032414) - this will only work for the first year, no?
    if (t == ncdstart)then
       currentSite%ncd = 0._r8
    endif

    !Accumulate growing/chilling days after start of counting period
    if (temp_in_C  <  ED_val_phen_chiltemp)then
       currentSite%ncd = currentSite%ncd + 1.0_r8
    endif

    !GDD accumulation function, which also depends on chilling days.
    gdd_threshold = ED_val_phen_a + ED_val_phen_b*exp(ED_val_phen_c*currentSite%ncd)

    !Accumulate temperature of last 10 days.
    currentSite%last_n_days(2:senes) =  currentSite%last_n_days(1:senes-1)
    currentSite%last_n_days(1) = temp_in_C                                      
    !count number of days for leaves off
    ncolddays = 0
    do i = 1,senes
       if (currentSite%last_n_days(i) < ED_val_phen_coldtemp)then
          ncolddays = ncolddays + 1
       endif
    enddo

    ! Here is where we do the GDD accumulation calculation
    !
    ! reset GDD on set dates
    if (t == gddstart)then
       currentSite%ED_GDD_site = 0._r8
    endif
    !
    ! accumulate the GDD using daily mean temperatures
    if (bc_in%t_veg24_si .gt. tfrz) then
       currentSite%ED_GDD_site = currentSite%ED_GDD_site + bc_in%t_veg24_si - tfrz
    endif
    

    timesinceleafoff = hlm_model_day - currentSite%leafoffdate
    !LEAF ON: COLD DECIDUOUS. Needs to
    !1) have exceeded the growing degree day threshold 
    !2) The leaves should not be on already
    !3) There should have been at least on chilling day in the counting period.  
    if (currentSite%ED_GDD_site > gdd_threshold)then
       if (currentSite%is_cold) then
          if (currentSite%ncd >= 1) then
             ! NOTE(bja, 2015-01) should leafondate = model_day to be consistent with leaf off?
             currentSite%is_cold = .false.      
             currentSite%leafondate = t         ! record leaf on date   
             if ( debug ) write(fates_log(),*) 'leaves on'
          endif !ncd
       endif !status
    endif !GDD

    timesinceleafon = hlm_model_day - currentSite%leafondate


    !LEAF OFF: COLD THRESHOLD
    !Needs to:
    !1) have exceeded the number of cold days threshold
    !2) have exceeded the minimum leafon time.
    !3) The leaves should not be off already
    !4) The day of the year should be larger than the counting period. 
    !   (not sure if we need this/if it will break the restarting)
    
    if (ncolddays > ED_val_phen_ncolddayslim)then
     if (timesinceleafon > ED_val_phen_mindayson)then
       if (.not.currentSite%is_cold)then
          currentSite%is_cold = .true.              ! now cold enough for leaf-off
          currentSite%leafoffdate = hlm_model_day   ! record leaf off date   
          if ( debug ) write(fates_log(),*) 'leaves off'
       endif
    endif
    endif

    !LEAF OFF: COLD LIFESPAN THRESHOLD
    if(timesinceleafoff > 400)then !remove leaves after a whole year when there is no 'off' period.  
       if(.not.currentSite%is_cold)then
          currentSite%is_cold = .true.              ! cold enough for leaf-off
          currentSite%leafoffdate = hlm_model_day   !record leaf off date   
          if ( debug ) write(fates_log(),*) 'leaves off'
       endif
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
    ! D*: We don't allow anything to happen in the first ten days to allow the water memory window to come into equlibirum. 
    ! E*: If the soil is always wet, the leaves come on at the beginning of the window, and then last for their lifespan. 
    ! ISSUES
    ! 1. It's not clear what water content we should track. Here we are tracking the top layer, 
    ! but we probably should track something like BTRAN,
    ! but BTRAN is defined for each PFT, and there could potentially be more than one stress-dec PFT.... ?
    ! 2. In the beginning, the window is set at an arbitrary time of the year, so the leaves might come on 
    ! in the dry season, using up stored reserves
    ! for the stress-dec plants, and potentially killing them. To get around this, we need to read in the 
    ! 'leaf on' date from some kind of start-up file
    ! but we would need that to happen for every resolution, etc. 
    ! 3. Will this methodology properly kill off the stress-dec trees where there is no water stress? 
    ! What about where the wet period coincides with the
    ! warm period? We would just get them overlapping with the cold-dec trees, even though that isn't appropriate.... 
    ! Why don't the drought deciduous trees grow
    ! in the North? Is cold decidousness maybe even the same as drought deciduosness there (and so does this 
    ! distinction actually matter??).... 

    !Accumulate surface water memory of last 10 days.

    do i = 1,numWaterMem-1 !shift memory along one
       currentSite%water_memory(numWaterMem+1-i) = currentSite%water_memory(numWaterMem-i)
    enddo
    currentSite%water_memory(1) = bc_in%h2o_liqvol_sl(1)   !waterstate_inst%h2osoi_vol_col(coli,1)

    ! In drought phenology, we often need to force the leaves to stay on or off as moisture fluctuates...      
    ! Here we incremend how long the leaves have been off;
    ! We set the default assumption that no time has elapsed, but if drought
    ! status is true, then we update the time
    ! If the leaves are off. How long have they been off? 
    ! leaves have come on, but last year, so at a later date than now.
    timesincedleafoff = 0
    if ( currentSite%is_drought )then
       if (currentSite%dleafoffdate > 0.and.currentSite%dleafoffdate > t)then 
          timesincedleafoff = t + (360 - currentSite%dleafoffdate)
       else
          timesincedleafoff = t - currentSite%dleafoffdate    
       endif
    endif

    timesincedleafon = 0
    !the leaves are on. How long have they been on? 
    if ( .not.currentSite%is_drought )then  
       !leaves have come on, but last year, so at a later date than now.
       if (currentSite%dleafondate > 0.and.currentSite%dleafondate > t)then 
          timesincedleafon = t + (360 - currentSite%dleafondate)
       else
          timesincedleafon = t - currentSite%dleafondate      
       endif
    endif

    !LEAF ON: DROUGHT DECIDUOUS WETNESS
    !Here, we used a window of oppurtunity to determine if we are close to the time when then leaves came on last year
    if ((t >= currentSite%dleafondate - 30.and.t <= currentSite%dleafondate + 30).or.(t > 360 - 15.and. &
         currentSite%dleafondate < 15))then ! are we in the window?

       if ( sum(currentSite%water_memory(1:numWaterMem))/real(numWaterMem,r8) &
            >= ED_val_phen_drought_threshold .and. &
            currentSite%is_drought .and. &
            (t >= 10) ) then 
          ! leave some minimum time between leaf off and leaf on to prevent 'flickering'.  
          if (timesincedleafoff > ED_val_phen_doff_time)then  
             currentSite%is_drought  = .false.     ! end the drought
             currentSite%dleafondate = t           !record leaf on date
          endif
       endif
    endif

    ! we still haven't done budburst by end of window
    if (t == currentSite%dleafondate+30 .and. currentSite%is_drought)then 
       currentSite%is_drought = .false.  ! force budburst!
       currentSite%dleafondate = t       ! record leaf on date
    endif

    !LEAF OFF: DROUGHT DECIDUOUS LIFESPAN - if the leaf gets to the end of its useful life. A*, E*
    if ( .not.currentSite%is_drought .and. (t >= 10) ) then  !D*
       !Are the leaves at the end of their lives? 
       if ( timesincedleafon > canopy_leaf_lifespan )then 
          currentSite%is_drought   = .true.   !alter status of site to 'leaves off'
          currentSite%dleafoffdate = t        !record leaf off date          
       endif
    endif

    !LEAF OFF: DROUGHT DECIDUOUS DRYNESS - if the soil gets too dry, and the leaves have already been on a while... 
    if ( .not.currentSite%is_drought .and. (t >= 10) ) then  !D*
       if (sum(currentSite%water_memory(1:10)/10._r8) <= ED_val_phen_drought_threshold)then 
          if (timesincedleafon > 100)then !B* Have the leaves been on for some reasonable length of time? To prevent flickering. 
             currentSite%is_drought   = .true.  !alter status of site to 'leaves on'
             currentSite%dleafoffdate = t       !record leaf on date           
          endif
       endif
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
    real(r8) :: store_c                ! storage carbon [kg]
    real(r8) :: store_c_transfer_frac  ! Fraction of storage carbon used to flush leaves
    integer  :: ipft
    real(r8), parameter :: leaf_drop_fraction = 1.0_r8

    !------------------------------------------------------------------------

    currentPatch => CurrentSite%oldest_patch   

    do while(associated(currentPatch))    
       currentCohort => currentPatch%tallest
       do while(associated(currentCohort))        

          ipft = currentCohort%pft

          ! Retrieve existing leaf and storage carbon

          call currentCohort%prt%CheckMassConservation(ipft,0)

          store_c = currentCohort%prt%GetState(store_organ, all_carbon_elements)
          leaf_c  = currentCohort%prt%GetState(leaf_organ, all_carbon_elements)

          ! COLD LEAF ON
          ! The site level flags signify that it is no-longer too cold
          ! for leaves. Time to signal flushing

          if (EDPftvarcon_inst%season_decid(ipft) == 1)then
             if ( .not.currentSite%is_cold )then            !we have just moved to leaves being on . 
                if (currentCohort%status_coh == leaves_off)then !Are the leaves currently off?        
                   currentCohort%status_coh = leaves_on         ! Leaves are on, so change status to 
                                                                ! stop flow of carbon out of bstore. 
                   
                   if(store_c>nearzero) then
                      store_c_transfer_frac = &
                            min(EDPftvarcon_inst%phenflush_fraction(ipft)*currentCohort%laimemory, store_c)/store_c
                   else
                      store_c_transfer_frac = 0.0_r8
                   end if

                   ! This call will request that storage carbon will be transferred to 
                   ! leaf tissues. It is specified as a fraction of the available storage
                   call PRTPhenologyFlush(currentCohort%prt, ipft, leaf_organ, store_c_transfer_frac)

                   currentCohort%laimemory = 0.0_r8

                endif !pft phenology
             endif ! growing season 

             ! COLD LEAF OFF
             ! The site level flag signifies that it is now too cold for 
             ! deciduous leaves. Time to drop if they have not already

             if (currentSite%is_cold) then 
                if (currentCohort%status_coh == leaves_on)then ! leaves have not dropped
                   
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
                   

                endif !leaf status
             endif !currentSite status
          endif  !season_decid

          ! DROUGHT LEAF ON
          ! Site level flag indicates it is no longer in drought condition
          ! deciduous plants can flush

          if (EDPftvarcon_inst%stress_decid(ipft) == 1)then
             
             if ( .not.currentSite%is_drought ) then 

                ! we have just moved to leaves being on . 
                if (currentCohort%status_coh == leaves_off)then    

                   !is it the leaf-on day? Are the leaves currently off?    

                   currentCohort%status_coh = leaves_on    ! Leaves are on, so change status to 
                                                           ! stop flow of carbon out of bstore. 

                   if(store_c>nearzero) then
                      store_c_transfer_frac = &
                            min(EDPftvarcon_inst%phenflush_fraction(ipft)*currentCohort%laimemory, store_c)/store_c
                   else
                      store_c_transfer_frac = 0.0_r8
                   end if
                   
                   ! This call will request that storage carbon will be transferred to 
                   ! leaf tissues. It is specified as a fraction of the available storage
                   call PRTPhenologyFlush(currentCohort%prt, ipft, &
                         leaf_organ, store_c_transfer_frac)

                   currentCohort%laimemory = 0.0_r8

                endif !currentCohort status again?
             endif   !currentSite status

             ! DROUGHT LEAF OFF
             ! Site level flag indicates a drought condition is in effect
             ! deciduous plants should drop leaves if have not already

             if ( currentSite%is_drought ) then        
                if (currentCohort%status_coh == leaves_on)then ! leaves have not dropped

                   ! This sets the cohort to the "leaves off" flag
                   currentCohort%status_coh      = leaves_off
                   
                   ! Remember what the lai (leaf mass actually) was for next year
                   currentCohort%laimemory   = leaf_c

                   call PRTDeciduousTurnover(currentCohort%prt,ipft, &
                         leaf_organ, leaf_drop_fraction)

                endif
             endif !status
          endif !drought dec.

          call currentCohort%prt%CheckMassConservation(ipft,1)

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
    ! -----------------------------------------------------------------------------------


    ! !USES:
    use EDTypesMod, only : area
    use EDTypesMod, only : homogenize_seed_pfts
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(bc_in_type), intent(in)               :: bc_in

    type(ed_patch_type), pointer  :: currentPatch
    type(litter_type), pointer    :: litt
    type(ed_cohort_type), pointer :: currentCohort

    integer  :: pft
    real(r8) :: store_m_to_repro       ! mass sent from storage to reproduction upon death [kg/plant]
    real(r8) :: site_seed_rain(maxpft) ! This is the sum of seed-rain for the site [kg/site/day]
    real(r8) :: mean_site_seed_rain    ! The mean site level seed rain for all PFTs
    integer  :: n_litt_types           ! number of litter element types (c,n,p, etc)
    integer  :: el                     ! loop counter for litter element types
    integer  :: element_id             ! element id consistent with parteh/PRTGenericMod.F90
    !------------------------------------------------------------------------------------

    do el = 1, num_elements
       
       site_seed_rain(:) = 0._r8
       
       site_mass => currentSite%mass_balance(el)

       ! Loop over all patches and sum up the seed input for each PFT
       currentPatch => currentSite%oldest_patch
       do while (associated(currentPatch))
          
          litt => currentPatch%litter(el)
          element_id = litt%element_id

          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             
             pft = currentCohort%pft
          
             ! a certain fraction of bstore might go to clonal reproduction when plants die
             store_m_to_repro = -currentCohort%prt%GetState(store_organ,element_id) * &
                   EDPftvarcon_inst%allom_frbstor_repro(pft)*currentCohort%dndt
             
             
             ! Transfer reproductive tissues from "on-plant" to the litter pool
             ! Transfer all reproductive tissues into seed production
             call PRTReproRelease(currentCohort%prt,repro_organ,element_id, &
                   1.0_r8, seed_prod)
             
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
          mean_site_seed_rain = sum(site_seed_rain(:))/real(numpft,r8)
          site_seed_rain(1:numpft) = mean_site_seed_rain
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
             seed_in_external =  seed_stoich*EDPftvarcon_inst%seed_rain(pft)/days_per_year
             
             litt%seed_in_extern(pft) = litt%seed_in_extern(pft) + seed_in_external

             site_mass%seed_influx = site_mass%seed_influx + seed_in_external*currentPatch%area

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
    type(litt_vartype) :: litt
    !
    ! !LOCAL VARIABLES:
    integer  ::  pft
    !----------------------------------------------------------------------

    ! default value from Liscke and Loffler 2006 ; making this a PFT-specific parameter
    ! decays the seed pool according to exponential model
    ! seed_decay_turnover is in yr-1
    ! seed_decay is kg/day
    ! Assume that decay rates are same for all chemical species

    do pft = 1,numpft 
       litt%seed_decay(pft) = litt%seed(pft) * &
             EDPftvarcon_inst%seed_decay_turnover(pft)/days_per_year
    enddo

    return
  end subroutine SeedDecay

  ! ============================================================================
  subroutine SeedGermination( litt, is_cold, is_drought )
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into sapling pool    
    !
    ! !USES:
    
    !
    ! !ARGUMENTS
    type(litt_vartype) :: litt  
    logical, intent(in) :: is_cold    ! Is the site in cold leaf-off status?
    logical, intent(in) :: is_drought ! Is the site in drought leaf-off status?
    !
    ! !LOCAL VARIABLES:
    integer :: pft

    ! Turning of this cap, because the cap will impose changes on proportionality
    ! of nutrients. (RGK 02-2019)
    !real(r8), parameter ::  max_germination 1.0_r8 ! Cap on germination rates. 
                                                    ! KgC/m2/yr Lishcke et al. 2009

    real(r8), parameter :: max_germination = 1.e6_r8  ! Force to very high number

    !----------------------------------------------------------------------

    ! germination_timescale is being pulled to PFT parameter; units are 1/yr
    ! thus the mortality rate of seed -> recruit (in units of carbon) 
    ! is seed_decay_turnover(p)/germination_timescale(p)
    ! and thus the mortlaity rate (in units of individuals) is the product of 
    ! that times the ratio of (hypothetical) seed mass to recruit biomass

    do pft = 1,numpft
       litt%seed_germ(pft) =  min(litt%seed(pft) * EDPftvarcon_inst%germination_timescale(pft), &
                                  max_germination)/days_per_year  
       
       !set the germination only under the growing season...c.xu
       if (EDPftvarcon_inst%season_decid(pft) == itrue .and. is_cold)then 
          litt%seed_germ(pft) = 0.0_r8
       endif
       if (EDPftvarcon_inst%stress_decid(pft) == itrue .and. is_drought)then
          litt%seed_germ(pft) = 0.0_r8
       endif
    enddo

  end subroutine SeedGermination

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
    integer :: ft
    type (ed_cohort_type) , pointer :: temp_cohort
    type (litter_type), pointer     :: litt          ! The litter object (carbon right now)
    integer :: cohortstatus
    integer,parameter :: recruitstatus = 1 !weather it the new created cohorts is recruited or initialized
    real(r8) :: b_leaf
    real(r8) :: b_fineroot    ! fine root biomass [kgC]
    real(r8) :: b_sapwood     ! sapwood biomass [kgC]
    real(r8) :: a_sapwood     ! sapwood cross section are [m2] (dummy)
    real(r8) :: b_agw         ! Above ground biomass [kgC]
    real(r8) :: b_bgw         ! Below ground biomass [kgC]
    real(r8) :: b_dead
    real(r8) :: b_store
    !----------------------------------------------------------------------

    allocate(temp_cohort) ! create temporary cohort
    call zero_cohort(temp_cohort)

    do ft = 1,numpft

       temp_cohort%canopy_trim = 0.8_r8  !starting with the canopy not fully expanded 
       temp_cohort%pft         = ft
       temp_cohort%hite        = EDPftvarcon_inst%hgt_min(ft)
       call h2d_allom(temp_cohort%hite,ft,temp_cohort%dbh)

       ! Initialize live pools
       call bleaf(temp_cohort%dbh,ft,temp_cohort%canopy_trim,b_leaf)
       call bfineroot(temp_cohort%dbh,ft,temp_cohort%canopy_trim,b_fineroot)
       call bsap_allom(temp_cohort%dbh,ft,temp_cohort%canopy_trim,a_sapwood, b_sapwood)
       call bagw_allom(temp_cohort%dbh,ft,b_agw)
       call bbgw_allom(temp_cohort%dbh,ft,b_bgw)
       call bdead_allom(b_agw,b_bgw,b_sapwood,ft,b_dead)
       call bstore_allom(temp_cohort%dbh,ft,temp_cohort%canopy_trim,b_store)

       ! Default assumption is that leaves are on
       cohortstatus = leaves_on
       temp_cohort%laimemory = 0.0_r8     

       ! But if the plant is seasonally (cold) deciduous, and the site status is flagged
       ! as "cold", then set the cohort's status to leaves_off, and remember the leaf biomass
       if (EDPftvarcon_inst%season_decid(ft) == itrue .and. currentSite%is_cold)then
          temp_cohort%laimemory = b_leaf
          b_leaf = 0.0_r8
          cohortstatus = leaves_off
       endif
       
       ! Or.. if the plant is drought deciduous, and the site status is flagged as 
       ! "in a drought", then likewise, set the cohort's status to leaves_off, and remember leaf
       ! biomass
       if (EDPftvarcon_inst%stress_decid(ft) == itrue .and. currentSite%is_drought )then
          temp_cohort%laimemory = b_leaf
          b_leaf = 0.0_r8
          cohortstatus = leaves_off
       endif

       ! This is somewhat (hacky), carbon12_element is index 1
       litt => currentPatch%litter(carbon12_element)

       if (hlm_use_ed_prescribed_phys .eq. ifalse .or. EDPftvarcon_inst%prescribed_recruitment(ft) .lt. 0. ) then
          temp_cohort%n           = currentPatch%area * litt%seed_germ(ft) &
               / (b_dead+b_leaf+b_fineroot+b_sapwood+b_store)
       else
          ! prescribed recruitment rates. number per sq. meter per year
          temp_cohort%n        = currentPatch%area * EDPftvarcon_inst%prescribed_recruitment(ft) * hlm_freq_day
       endif

       ! Only bother allocating a new cohort if there is a reasonable amount of it
       if (temp_cohort%n > min_n_safemath )then
          if ( debug ) write(fates_log(),*) 'EDPhysiologyMod.F90 call create_cohort '

          call create_cohort(currentSite,currentPatch, temp_cohort%pft, temp_cohort%n, temp_cohort%hite, temp_cohort%dbh, &
               b_leaf, b_fineroot, b_sapwood, b_dead, b_store, &  
               temp_cohort%laimemory, cohortstatus,recruitstatus, temp_cohort%canopy_trim, currentPatch%NCL_p, &
               currentSite%spread, first_leaf_aclass, bc_in)

          ! Note that if hydraulics is on, the number of cohorts may had changed due to hydraulic constraints.
          ! This constaint is applied during "create_cohort" subroutine.
          
          ! keep track of how many individuals were recruited for passing to history
          currentSite%recruitment_rate(ft) = currentSite%recruitment_rate(ft) + temp_cohort%n
          
          ! modify the carbon balance accumulators to take into account the different way of defining recruitment
          ! add prescribed rates as an input C flux, and the recruitment that would have otherwise occured as an output flux
          ! (since the carbon associated with them effectively vanishes)
          ! check the water for hydraulics
          if (hlm_use_ed_prescribed_phys .ne. ifalse .and. EDPftvarcon_inst%prescribed_recruitment(ft) .ge. 0. ) then
             currentSite%flux_in = currentSite%flux_in + temp_cohort%n * &
                  (b_store + b_leaf + b_fineroot + b_sapwood + b_dead)
             currentSite%flux_out = currentSite%flux_out + currentPatch%area * litt%seed_germ(ft)
          endif

       endif
    enddo  !pft loop

    deallocate(temp_cohort) ! delete temporary cohort

  end subroutine recruitment

  ! ============================================================================

  subroutine CWDInput( currentSite, currentPatch, litt, element_id, bc_in)

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
    type(dead_type),intent(inout),target      :: litt
    integer,intent(in)                        :: element_id
    type(bc_in_type), intent(in)              :: bc_in


    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: currentCohort
    type(site_fluxdiags_type), pointer :: flux_diags
    type(site_masscheck_type), pointer :: site_mass
    integer  :: c
    real(r8) :: dead_n          ! total understorey dead tree density
    real(r8) :: dead_n_dlogging ! direct logging understory dead-tree density
    real(r8) :: dead_n_ilogging ! indirect understory dead-tree density (logging)
    real(r8) :: dead_n_natural  ! understory dead density not associated
                                ! with direct logging
    real(r8) :: leaf_m          ! mass of the element of interest in the 
                                ! leaf  [kg]
    real(r8) :: fnrt_m
    real(r8) :: sapw_m
    real(r8) :: struct_m
    real(r8) :: store_m
    real(r8) :: leaf_m_turnover ! leaf turnover [kg]
    real(r8) :: fnrt_m_turnover
    real(r8) :: sapw_m_turnover
    real(r8) :: struct_m_turnover
    real(r8) :: store_m_turnover
    real(r8) :: bg_cwd_tot        ! Total below-ground coarse woody debris
                                  ! input flux
    real(r8) :: root_fines_tot    ! Total below-ground fine root coarse
                                  ! woody debris
    integer  :: element_id        ! element id consistent with parteh/PRTGenericMod.F90

    real(r8) :: trunk_product   ! carbon flux into trunk products kgC/day/site
    integer  :: pft
    integer  :: numlevsoil              ! Actual number of soil layers
    real(r8) :: rootfr(numlevsoil_max)  ! Fractional root mass profile
    !----------------------------------------------------------------------

    ! -----------------------------------------------------------------------------------
    ! Other direct litter fluxes happen in phenology and in spawn_patches. 
    ! -----------------------------------------------------------------------------------

    numlevsoil = bc_in%nlevsoil
    
    element_id = litt%element_id
    
    ! Object tracking flux diagnostics for each element
    flux_diags => currentSite%flux_diags(element_pos(element_id))
    
    ! Object tracking site level mass balance for each element
    site_mass => currentSite%mass_balance(element_pos(element_id))

    currentCohort => currentPatch%shortest
    do while(associated(currentCohort))
      pft = currentCohort%pft        

      leaf_m_turnover   = currentCohort%prt%GetTurnover(leaf_organ,element_id)
      store_m_turnover  = currentCohort%prt%GetTurnover(store_organ,element_id)
      fnrt_m_turnover   = currentCohort%prt%GetTurnover(fnrt_organ,element_id)
      sapw_m_turnover   = currentCohort%prt%GetTurnover(sapw_organ,element_id)
      struct_m_turnover = currentCohort%prt%GetTurnover(struct_organ,element_id)

      leaf_m          = currentCohort%prt%GetState(leaf_organ,element_id)
      store_m         = currentCohort%prt%GetState(store_organ,element_id)
      fnrt_m          = currentCohort%prt%GetState(fnrt_organ,element_id)
      sapw_m          = currentCohort%prt%GetState(sapw_organ,element_id)
      struct_m        = currentCohort%prt%GetState(struct_organ,element_id)

      plant_dens =  currentCohort%n/currentPatch%area

      ! ---------------------------------------------------------------------------------
      ! PART 1 Litter fluxes from non-mortal tissue turnovers  Kg/m2/day
      !        Important note:  Turnover has already been removed from the cohorts.
      !        So, in the next part of this algorithm, when we send the biomass
      !        from dying trees to the litter pools, we don't have to worry 
      !        about double counting.
      ! ---------------------------------------------------------------------------------

      litt%leaf_fines_in(pft) = litt%leaf_fines_in(pft) + & 
           leaf_m_turnover * plant_dens

      flux_diags%leaf_litter_input(pft) = &
            flux_diags%leaf_litter_input(pft) +  &
            leaf_m_turnover * currentCohort%n
      
      root_fines_tot = (fnrt_m_turnover + store_m_turnover ) * &
            plant_dens
      
      do ilyr = 1, numlevsoil
         litt%root_fines_in(ipft,ilyr) = litt%root_fines_in(ipft,ilyr) + &
              currentCohort%root_fr(ilyr) * root_fines_tot
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
                  bg_cwd_tot * currentCohort%root_fr(ilyr)
         end do
         
         flux_diags%cwd_bg_input(c)  = flux_diags%cwd_bg_input(c) + &
               bg_cwd_tot*currentPatch%area
         

      enddo


      ! ---------------------------------------------------------------------------------
      ! PART 2 Litter fluxes from non-disturbance inducing  mortality. Kg/m2/day
      ! ---------------------------------------------------------------------------------

      ! Total number of dead (n/m2/day)
      dead_n = -1.0_r8 * currentCohort%dndt/currentPatch%area / days_per_year
      
      ! Total number of dead understory from direct logging (n/m2/day)
      ! (it is possible that large harvestable trees are in the understory)
      dead_n_dlogging = currentCohort%lmort_direct*currentCohort%n/currentPatch%area
          
      ! Total number of dead understory from indirect logging
      dead_n_ilogging = ( currentCohort%lmort_collateral + currentCohort%lmort_infra) * &
           currentCohort%n/currentPatch%area
          
      dead_n_natural = dead_n - dead_n_dlogging - dead_n_ilogging

     
      litt%leaf_fines_in(pft) = litt%leaf_fines_in(pft) + leaf_m * dead_n

      flux_diags%leaf_litter_input(pft) = &
            flux_diags%leaf_litter_input(pft) +  &
            leaf_m * dead_n*currentPatch%area


      ! %n has not been updated due to mortality yet, thus
      ! the litter flux has already been counted since it captured
      ! the losses of live trees and those flagged for death
      
      root_fines_tot =  dead_n * (fnrt_m + &
           store_m*(1._r8-EDPftvarcon_inst%allom_frbstor_repro(pft)) )
      
      do ilyr = 1, numlevsoil
         litt%root_fines_in(pft,ilyr) = litt%root_fines_in(pft,ilyr) + &
              root_fines_tot * currentCohort%root_fr(ilyr)
      end do

      flux_diags%root_litter_input(pft) = &
            flux_diags%root_litter_input(pft) +  &
            root_fines_tot*currentPatch%area



      ! Track CWD inputs from mortal plants
      
      do c = 1,ncwd

         ! Below-ground 
         
         bg_cwd_tot = (struct_m + sapw_m) * & 
              SF_val_CWD_frac(c) * dead_n * &
              (1.0_r8-EDPftvarcon_inst%allom_agb_frac(pft))
         
         do ilyr = 1, numlevsoil
            litt%bg_cwd_in(c,ilyr) = litt%bg_cwd_in(c,ilyr) + &
                  currentCohort%root_fr(ilyr) * bg_cwd_tot
         end do

         flux_diags%cwd_bg_input(c)  = flux_diags%cwd_bg_input(c) + &
               SF_val_CWD_frac(c) * dead_n * currentPatch%area * &
               (1.0_r8-EDPftvarcon_inst%allom_agb_frac(pft))

         ! The bole of the plant is the last index of the cwd array. So any harvesting
         ! mortality is diverted away from above-ground CWD and sent to harvest
         ! and flux out. Send AGB component of boles from non direct-logging activities 
         ! to AGB litter pool

         if (c==ncwd) then
            
            litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + (struct_m + sapw_m) * & 
                 SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging)  * &
                 EDPftvarcon_inst%allom_agb_frac(pft)
            
            ! Send AGB component of boles from direct-logging activities to export/harvest pool
            ! Generate trunk product (kg/day/m2)
            trunk_product =  (struct_m + sapw_m) * &
                 SF_val_CWD_frac(c) * dead_n_dlogging * EDPftvarcon_inst%allom_agb_frac(pft)
            
            site_mass%wood_product = site_mass%wood_product + trunk_product
            
            flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
                  SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging)  * &
                  currentPatch%area * EDPftvarcon_inst%allom_agb_frac(pft)

         else

            litt%ag_cwd_in(c) = litt%ag_cwd_in(c) + (struct_m + sapw_m) * & 
                 SF_val_CWD_frac(c) * dead_n  * &
                 EDPftvarcon_inst%allom_agb_frac(pft)

            flux_diags%cwd_ag_input(c)  = flux_diags%cwd_ag_input(c) + &
                  SF_val_CWD_frac(c) * dead_n * &
                  currentPatch%area * EDPftvarcon_inst%allom_agb_frac(pft)
            
         end if
         
      end do


      ! Update diagnostics that track resource management

      if( (element_id .eq. all_carbon_elements) .or. &
            (element_id .eq. carbon12_element) ) then
      
         currentSite%resources_management%delta_litter_stock  = &
              currentSite%resources_management%delta_litter_stock + &
              (leaf_m + fnrt_m + store_m ) * &
              (dead_n_ilogging+dead_n_dlogging) * currentPatch%area

         currentSite%resources_management%delta_biomass_stock = &
              currentSite%resources_management%delta_biomass_stock + &
              (leaf_c + fnrt_c + store_c ) * & 
              (dead_n_ilogging+dead_n_dlogging) *currentPatch%area

         currentSite%resources_management%trunk_product_site = &
               currentSite%resources_management%trunk_product_site + &
               trunk_product * currentPatch%area

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

    type(litt_vartype) :: litt
    !
    ! !LOCAL VARIABLES:
    integer  ::  pft

    ! Add decaying seeds to the leaf litter
    ! -----------------------------------------------------------------------------------
    
    do pft = 1,numpft
       litt%leaf_fines_in(pft) = litt%leaf_fines_in(pft) + litt%seed_decay(pft)
    enddo
    
    
    return
  end subroutine SeedDecayToFines
  

  ! ============================================================================

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
    real(r8) :: Q10                   ! temperature dependence
    real(r8) :: froz_q10              ! separate q10 for frozen soil respiration rates.
                                      ! default to same as above zero rates
    !----------------------------------------------------------------------

    catanf(t1) = 11.75_r8 +(29.7_r8 / pi) * atan( pi * 0.031_r8  * ( t1 - 15.4_r8 ))
    catanf_30 = catanf(30._r8)
    
    ifp = currentPatch%patchno 
    
    ! set "froz_q10" parameter
    froz_q10  = FatesSynchronizedParamsInst%froz_q10  
    Q10       = FatesSynchronizedParamsInst%Q10

    if ( .not. use_century_tfunc ) then
    !calculate rate constant scalar for soil temperature,assuming that the base rate constants 
    !are assigned for non-moisture limiting conditions at 25C. 
      if (bc_in%t_veg24_pa(ifp)  >=  tfrz) then
        t_scalar = Q10**((bc_in%t_veg24_pa(ifp)-(tfrz+25._r8))/10._r8)
                 !  Q10**((t_soisno(c,j)-(tfrz+25._r8))/10._r8)
      else
        t_scalar = (Q10**(-25._r8/10._r8))*(froz_q10**((bc_in%t_veg24_pa(ifp)-tfrz)/10._r8))
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
    w_scalar = sum(currentPatch%btran_ft(1:numpft))/numpft

    currentPatch%fragmentation_scaler =  min(1.0_r8,max(0.0_r8,t_scalar * w_scalar))
    
  end subroutine fragmentation_scaler
  
  ! ============================================================================

  subroutine CWDOut( litt, element_id, fragmentation_scaler, nlev_eff_decomp )
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use SFParamsMod, only : SF_val_max_decomp

    !
    ! !ARGUMENTS    
    type(dead_type),intent(inout),target       :: litt
    integer,intent(in)                         :: element_id
    real(r8),intent(in)                        :: fragmentation_scaler

    ! This is not necessarily every soil layer, this is the number
    ! of effective layers that are active and can be sent
    ! to the soil decomposition model
    integer,intent(in)                         :: nlev_eff_decomp  
    
    !
    ! !LOCAL VARIABLES:
    integer :: c
    integer :: pft
    integer :: ilyr
    integer :: numlevsoil
    !----------------------------------------------------------------------

    do c = 1,ncwd  

       litt%ag_cwd_frag(c)   = litt%ag_cwd(c) * SF_val_max_decomp(c+1) * &
             fragmentation_scaler

       do ilyr = 1,nlev_eff_decomp

          litt%bg_cwd_frag(c,ilyr) = litt%bg_cwd(c,ilyr) * SF_val_max_decomp(c+1) * &
                fragmentation_scaler

       enddo
    end do

    ! this is the rate at which dropped leaves stop being part of the burnable pool 
    ! and begin to be part of the decomposing pool. This should probably be highly 
    ! sensitive to moisture, but also to the type of leaf thick leaves can dry out 
    ! before they are decomposed, for example. This section needs further scientific input. 

    do pft = 1,numpft
       
       litt%leaf_fines_frag(pft) = litt%leaf_litter(pft) * &
             SF_val_max_decomp(dl_sf) * fragmentation_scaler

       do ilyr = 1,nlev_eff_decomp
          litt%root_fines_frag(pft,ilyr) = litt%root_fines(pft,ilyr) * &
                SF_val_max_decomp(dl_sf) * fragmentation_scaler
       end do
    enddo

  end subroutine CWDOut

  ! =====================================================================================

  subroutine flux_into_litter_pools(nsites, sites, bc_in, bc_out)
    
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
    use FatesAllometryMod, only : set_root_fraction
    use FatesAllometryMod, only : i_biomass_rootprof_context 
    

    implicit none   

    ! !ARGUMENTS    
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(:)
    type(bc_out_type)       , intent(inout)         :: bc_out(:)

    ! !LOCAL VARIABLES:
    type (ed_patch_type),  pointer :: currentPatch
    type (ed_cohort_type), pointer :: currentCohort
    real(r8),dimension(:), pointer :: flux_cel_si
    real(r8),dimension(:), pointer :: flux_lab_si
    real(r8),dimension(:), pointer :: flux_lig_si
    type(litter_type),     pointer :: litt
     
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
             flux_cel_si => bc_out(s)%litt_flux_cel_c_si
             flux_lab_si => bc_out(s)%litt_flux_lab_c_si
             flux_lig_si => bc_out(s)%litt_flux_lig_c_si
          case (nitrogen_element) 
             bc_out(s)%litt_flux_cel_n_si(:) = 0._r8
             bc_out(s)%litt_flux_lig_n_si(:) = 0._r8
             bc_out(s)%litt_flux_lab_n_si(:) = 0._r8
             flux_cel_si => bc_out(s)%litt_flux_cel_n_si
             flux_lab_si => bc_out(s)%litt_flux_lab_n_si
             flux_lig_si => bc_out(s)%litt_flux_lig_n_si
          case (phosphorus_element)
             bc_out(s)%litt_flux_cel_p_si(:) = 0._r8
             bc_out(s)%litt_flux_lig_p_si(:) = 0._r8
             bc_out(s)%litt_flux_lab_p_si(:) = 0._r8
             flux_cel_si => bc_out(s)%litt_flux_cel_p_si
             flux_lab_si => bc_out(s)%litt_flux_lab_p_si
             flux_lig_si => bc_out(s)%litt_flux_lig_p_si
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
                   
                   flux_lig_si(id) = flux_lig_si(id) & 
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
             do ft = 1,numpft
                
                do id = 1,nlev_eff_decomp
                   
                   flux_lab_si(id) = flux_lab_si(id) + &
                         litt%leaf_fines_frag(ft) * EDPftvarcon_inst%lf_flab(ft) * &
                         area_frac* surface_prof(id)
                   
                   flux_cel_si(id) = flux_cel_si(id) + &
                         litt%leaf_fines_frag(ft) * EDPftvarcon_inst%lf_fcel(ft) * &
                         area_frac* surface_prof(id)
                   
                   flux_lig_si(id) = flux_lig_si(id) + &
                         litt%leaf_fines_frag(ft) * EDPftvarcon_inst%lf_flig(ft) * &
                         area_frac* surface_prof(id)

                end do

                do j = 1, nlev_eff_soil
                   
                   id = bc_in(s)%decomp_id(j)
                   
                   bc_out(s)%litt_flux_lab_c_si(id) = bc_out(s)%litt_flux_lab_c_si(id) + &
                         litt%root_fines_frag(ft,j) * EDPftvarcon_inst%fr_flab(ft) * area_frac
                   
                   bc_out(s)%litt_flux_cel_c_si(id) = bc_out(s)%litt_flux_cel_c_si(id) + &
                         litt%root_fines_frag(ft,j) * EDPftvarcon_inst%fr_fcel(ft) * area_frac
                   
                   bc_out(s)%litt_flux_lig_c_si(id) = bc_out(s)%litt_flux_lig_c_si(id) + &
                         litt%root_fines_frag(ft,j) * EDPftvarcon_inst%fr_flig(ft) * area_frac
                enddo
             end do
               
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
 end subroutine flux_into_litter_pools



end module EDPhysiologyMod
