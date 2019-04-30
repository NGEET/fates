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
  use EDPftvarcon      , only    : EDPftvarcon_inst
  use FatesInterfaceMod, only    : bc_in_type
  use EDCohortDynamicsMod , only : zero_cohort
  use EDCohortDynamicsMod , only : create_cohort, sort_cohorts,InitPRTCohort
  use FatesAllometryMod   , only : tree_lai
  use FatesAllometryMod   , only : tree_sai
  use FatesAllometryMod   , only : decay_coeff_kn

  use EDTypesMod          , only : numWaterMem
  use EDTypesMod          , only : dl_sf, dinc_ed
  use EDTypesMod          , only : ncwd
  use EDTypesMod          , only : nlevleaf
  use EDTypesMod          , only : num_vegtemp_mem
  use EDTypesMod          , only : maxpft
  use EDTypesMod          , only : ed_site_type, ed_patch_type, ed_cohort_type
  use EDTypesMod          , only : dump_cohort
  use EDTypesMod          , only : first_leaf_aclass
  use EDTypesMod          , only : leaves_on
  use EDTypesMod          , only : leaves_off
  use EDTypesMod          , only : min_n_safemath
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

  use FatesPlantHydraulicsMod  , only : AccumulateMortalityWaterStorage
  
  use FatesConstantsMod     , only : itrue,ifalse
  use FatesConstantsMod     , only : calloc_abs_error

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
  
  use PRTGenericMod, only : prt_carbon_allom_hyp
  use PRTGenericMod, only : leaf_organ
  use PRTGenericMod, only : all_carbon_elements
  use PRTGenericMod, only : carbon12_element
  use PRTGenericMod, only : nitrogen_element
  use PRTGenericMod, only : phosphorous_element
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
  private :: seeds_in
  private :: seed_decay
  private :: seed_germination
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
  integer, parameter :: dleafon_drycheck = 100 ! Drought deciduous leaves max days on check parameter 

  ! ============================================================================

contains

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

  subroutine non_canopy_derivs( currentSite, currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! Returns time differentials of the state vector
    !
    ! !USES:
    use EDTypesMod, only : AREA
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout)         :: currentPatch
    type(bc_in_type), intent(in)               :: bc_in

    !
    ! !LOCAL VARIABLES:
    integer c,p
    !----------------------------------------------------------------------

    currentPatch%leaf_litter_in(:)   = 0.0_r8
    currentPatch%root_litter_in(:)   = 0.0_r8
    currentPatch%dleaf_litter_dt(:)  = 0.0_r8
    currentPatch%droot_litter_dt(:)  = 0.0_r8
    currentPatch%leaf_litter_out(:)  = 0.0_r8
    currentPatch%root_litter_out(:)  = 0.0_r8
    currentPatch%cwd_AG_in(:)        = 0.0_r8
    currentPatch%cwd_BG_in(:)        = 0.0_r8
    currentPatch%cwd_AG_out(:)       = 0.0_r8
    currentPatch%cwd_BG_out(:)       = 0.0_r8
    currentPatch%seeds_in(:)         = 0.0_r8  
    currentPatch%seed_decay(:)       = 0.0_r8
    currentPatch%seed_germination(:) = 0.0_r8

    ! update seed fluxes 
    call seeds_in(currentSite, currentPatch)
    call seed_decay(currentSite, currentPatch)
    call seed_germination(currentSite, currentPatch)

    ! update fragmenting pool fluxes
    call cwd_input( currentSite, currentPatch)
    call cwd_out( currentSite, currentPatch, bc_in)
  
    do p = 1,numpft
       currentSite%dseed_dt(p) = currentSite%dseed_dt(p) + &
            (currentPatch%seeds_in(p) - currentPatch%seed_decay(p) - &
            currentPatch%seed_germination(p)) * currentPatch%area/AREA
    enddo   
    
    do c = 1,ncwd
       currentPatch%dcwd_AG_dt(c) = currentPatch%cwd_AG_in(c) - currentPatch%cwd_AG_out(c) 
       currentPatch%dcwd_BG_dt(c) = currentPatch%cwd_BG_in(c) - currentPatch%cwd_BG_out(c) 
    enddo

    do p = 1,numpft
       currentPatch%dleaf_litter_dt(p) = currentPatch%leaf_litter_in(p) - &
             currentPatch%leaf_litter_out(p) 
       currentPatch%droot_litter_dt(p) = currentPatch%root_litter_in(p) - &
             currentPatch%root_litter_out(p) 
    enddo

  end subroutine non_canopy_derivs

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
    if (bc_in%t_veg24_si .gt. tfrz) then
       currentSite%grow_deg_days = currentSite%grow_deg_days + bc_in%t_veg24_si - tfrz
    endif
    
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
         (currentSite%nchilldays >= 1)) then
       currentSite%cstatus = phen_cstat_notcold  ! Set to not-cold status (leaves can come on)
       currentSite%cleafondate = model_day_int  
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

          if (EDPftvarcon_inst%season_decid(ipft) == itrue)then
             if ( currentSite%cstatus == phen_cstat_notcold  )then                ! we have just moved to leaves being on . 
                if (currentCohort%status_coh == leaves_off)then ! Are the leaves currently off?        
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

             !COLD LEAF OFF
             if (currentSite%cstatus == phen_cstat_nevercold .or. &
                 currentSite%cstatus == phen_cstat_iscold) then ! past leaf drop day? Leaves still on tree?  

                if (currentCohort%status_coh == leaves_on) then ! leaves have not dropped

                   
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

                endif
             endif !status
          endif !drought dec.

          call currentCohort%prt%CheckMassConservation(ipft,1)

          currentCohort => currentCohort%shorter
       enddo !currentCohort

       currentPatch => currentPatch%younger

    enddo !currentPatch

  end subroutine phenology_leafonoff


  ! ============================================================================
  subroutine seeds_in( currentSite, cp_pnt )
    !
    ! !DESCRIPTION:
    !  Flux from plants into seed pool. 
    !
    ! !USES:
    use EDTypesMod, only : AREA
    use EDTypesMod, only : homogenize_seed_pfts
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout), target :: cp_pnt ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type),  pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort
    integer :: p
    logical :: pft_present(maxpft)
    real(r8) :: store_c_to_repro   ! carbon sent from storage to reproduction upon death [kg/plant]
    real(r8) :: npfts_present
    !----------------------------------------------------------------------

    currentPatch => cp_pnt
   
    currentPatch%seeds_in(:) = 0.0_r8

    if ( homogenize_seed_pfts ) then
       ! special mode to remove intergenerational filters on PFT existence: each PFT seeds all PFTs
       ! first loop over all patches and cohorts to see what and how many PFTs are present on this site
       pft_present(:) = .false.
       npfts_present =  0._r8
       currentPatch => currentSite%oldest_patch
       do while(associated(currentPatch))
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             p = currentCohort%pft
             if (.not. pft_present(p)) then
                pft_present(p) = .true.
                npfts_present = npfts_present + 1._r8
             endif
             currentCohort => currentCohort%shorter
          enddo !cohort loop                        
          currentPatch => currentPatch%younger
       enddo ! patch loop
       
       ! now calculate the homogenized seed flux into each PFT pool
       currentPatch => cp_pnt
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))

          ! a certain fraction of bstore goes to clonal reproduction when plants die
          store_c_to_repro = currentCohort%prt%GetState(store_organ,all_carbon_elements) * &
                EDPftvarcon_inst%allom_frbstor_repro(currentCohort%pft)
          
          do p = 1, numpft
             if (pft_present(p)) then
		  
                  currentPatch%seeds_in(p) = currentPatch%seeds_in(p) + &
                        (currentCohort%seed_prod * currentCohort%n - &
                        currentCohort%dndt*store_c_to_repro) &
                        /(currentPatch%area * npfts_present)		  
             endif
          end do
          currentCohort => currentCohort%shorter
       enddo !cohort loop                  
    else

    ! normal case: each PFT seeds its own type
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort))
       p = currentCohort%pft

       ! a certain fraction of bstore goes to clonal reproduction when plants die
       store_c_to_repro = currentCohort%prt%GetState(store_organ,all_carbon_elements) * &
             EDPftvarcon_inst%allom_frbstor_repro(p)
       
       currentPatch%seeds_in(p) = currentPatch%seeds_in(p) +  &
           (currentCohort%seed_prod * currentCohort%n - &
	   currentCohort%dndt*store_c_to_repro)/currentPatch%area

       currentCohort => currentCohort%shorter
    enddo !cohort loop

    endif

    do p = 1,numpft
       currentPatch%seeds_in(p) = currentPatch%seeds_in(p) + &
                 EDPftvarcon_inst%seed_rain(p) !KgC/m2/year
       currentSite%seed_rain_flux(p) = currentSite%seed_rain_flux(p) + &
                 EDPftvarcon_inst%seed_rain(p) * currentPatch%area/AREA !KgC/m2/year

       currentSite%flux_in = currentSite%flux_in + &
             EDPftvarcon_inst%seed_rain(p) * currentPatch%area * hlm_freq_day

    enddo


  end subroutine seeds_in
  
  ! ============================================================================
  subroutine seed_decay( currentSite, currentPatch )
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into leaf litter pool    
    !
    ! !USES:
    use EDPftvarcon       , only : EDPftvarcon_inst
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer  ::  p
    !----------------------------------------------------------------------

    ! default value from Liscke and Loffler 2006 ; making this a PFT-specific parameter
    ! decays the seed pool according to exponential model
    ! seed_decay_turnover is in yr-1
    do p = 1,numpft 
       currentPatch%seed_decay(p) =  currentSite%seed_bank(p) * EDPftvarcon_inst%seed_decay_turnover(p)
    enddo
 
  end subroutine seed_decay

  ! ============================================================================
  subroutine seed_germination( currentSite, currentPatch ) 
    !
    ! !DESCRIPTION:
    !  Flux from seed pool into sapling pool    
    !
    ! !USES:
    use EDPftvarcon       , only : EDPftvarcon_inst
    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type),intent(inout) :: currentPatch ! seeds go to these patches.
    !
    ! !LOCAL VARIABLES:
    integer :: p
    real(r8) max_germination !cap on germination rates. KgC/m2/yr Lishcke et al. 2009
    !----------------------------------------------------------------------

    max_germination = 1.0_r8 !this is arbitrary

    ! germination_timescale is being pulled to PFT parameter; units are 1/yr
    ! thus the mortality rate of seed -> recruit (in units of carbon) 
    ! is seed_decay_turnover(p)/germination_timescale(p)
    ! and thus the mortlaity rate (in units of individuals) is the product of 
    ! that times the ratio of (hypothetical) seed mass to recruit biomass

    do p = 1,numpft
       currentPatch%seed_germination(p) =  min(currentSite%seed_bank(p) * &
             EDPftvarcon_inst%germination_timescale(p),max_germination)     
       !set the germination only under the growing season...c.xu
       if ( (EDPftvarcon_inst%season_decid(p) == itrue) .and. &
            (any(currentSite%cstatus == [phen_cstat_nevercold,phen_cstat_iscold]))) then
          currentPatch%seed_germination(p) = 0.0_r8
       endif
       if ( (EDPftvarcon_inst%stress_decid(p) == itrue) .and. & 
            (any(currentSite%dstatus == [phen_dstat_timeoff,phen_dstat_moistoff]))) then
          currentPatch%seed_germination(p) = 0.0_r8
       endif
    enddo

  end subroutine seed_germination

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

       if ( (EDPftvarcon_inst%season_decid(temp_cohort%pft) == itrue) .and. &
            (any(currentSite%cstatus == [phen_cstat_nevercold,phen_cstat_iscold]))) then
          temp_cohort%laimemory = b_leaf
          b_leaf = 0.0_r8
          cohortstatus = leaves_off
       endif

       if ( (EDPftvarcon_inst%stress_decid(temp_cohort%pft) == itrue) .and. &
            (any(currentSite%dstatus == [phen_dstat_timeoff,phen_dstat_moistoff]))) then
          temp_cohort%laimemory = b_leaf
          b_leaf = 0.0_r8
          cohortstatus = leaves_off
       endif


       if (hlm_use_ed_prescribed_phys .eq. ifalse .or. EDPftvarcon_inst%prescribed_recruitment(ft) .lt. 0. ) then
          temp_cohort%n           = currentPatch%area * currentPatch%seed_germination(ft)*hlm_freq_day &
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
             currentSite%flux_out = currentSite%flux_out + currentPatch%area * currentPatch%seed_germination(ft)*hlm_freq_day
          endif


       endif
    enddo  !pft loop

    deallocate(temp_cohort) ! delete temporary cohort

  end subroutine recruitment

  ! ============================================================================
  subroutine CWD_Input( currentSite, currentPatch)
    !
    ! !DESCRIPTION:
    ! Generate litter fields from turnover.  
    !
    ! !USES:
    use SFParamsMod , only : SF_val_CWD_frac

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target :: currentSite
    type(ed_patch_type),intent(inout), target :: currentPatch
    !
    ! !LOCAL VARIABLES:
    type(ed_cohort_type), pointer :: currentCohort
    integer  :: c,p
    real(r8) :: dead_n          ! total understorey dead tree density
    real(r8) :: dead_n_dlogging ! direct logging understory dead-tree density
    real(r8) :: dead_n_ilogging ! indirect understory dead-tree density (logging)
    real(r8) :: dead_n_natural  ! understory dead density not associated
                                ! with direct logging
    real(r8) :: leaf_c          ! leaf carbon [kg]
    real(r8) :: fnrt_c
    real(r8) :: sapw_c
    real(r8) :: struct_c
    real(r8) :: store_c
    real(r8) :: leaf_c_turnover ! leaf turnover [kg]
    real(r8) :: fnrt_c_turnover
    real(r8) :: sapw_c_turnover
    real(r8) :: struct_c_turnover
    real(r8) :: store_c_turnover

    real(r8) :: trunk_product   ! carbon flux into trunk products kgC/day/site
    integer  :: pft
    !----------------------------------------------------------------------

    ! ================================================        
    ! Other direct litter fluxes happen in phenology and in spawn_patches. 
    ! ================================================   

    currentCohort => currentPatch%shortest

    do while(associated(currentCohort))
      pft = currentCohort%pft        

      leaf_c_turnover   = currentCohort%prt%GetTurnover(leaf_organ,all_carbon_elements)
      store_c_turnover  = currentCohort%prt%GetTurnover(store_organ,all_carbon_elements)
      fnrt_c_turnover   = currentCohort%prt%GetTurnover(fnrt_organ,all_carbon_elements)
      sapw_c_turnover   = currentCohort%prt%GetTurnover(sapw_organ,all_carbon_elements)
      struct_c_turnover = currentCohort%prt%GetTurnover(struct_organ,all_carbon_elements)

      leaf_c          = currentCohort%prt%GetState(leaf_organ,all_carbon_elements)
      store_c         = currentCohort%prt%GetState(store_organ,all_carbon_elements)
      fnrt_c          = currentCohort%prt%GetState(fnrt_organ,all_carbon_elements)
      sapw_c          = currentCohort%prt%GetState(sapw_organ,all_carbon_elements)
      struct_c        = currentCohort%prt%GetState(struct_organ,all_carbon_elements)

      ! ================================================        
      ! Litter from tissue turnover. KgC/m2/year
      ! ================================================   

      currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
           leaf_c_turnover * currentCohort%n/currentPatch%area/hlm_freq_day
      
      currentPatch%root_litter_in(pft) = currentPatch%root_litter_in(pft) + &
           (fnrt_c_turnover + store_c_turnover ) * &
           currentCohort%n/currentPatch%area/hlm_freq_day
      
      
      !daily leaf loss needs to be scaled up to the annual scale here. 

      ! ---------------------------------------------------------------------------------
      ! Assumption: turnover from deadwood and sapwood are lumped together in CWD pool
      ! ---------------------------------------------------------------------------------
      
      do c = 1,ncwd
         currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + &
              (sapw_c_turnover + struct_c_turnover)/hlm_freq_day   * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area * &
              EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)

         currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + &
              (sapw_c_turnover + struct_c_turnover)/hlm_freq_day * &
              SF_val_CWD_frac(c) * currentCohort%n/currentPatch%area * &
              (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft))
      enddo

      !if (currentCohort%canopy_layer > 1)then   

          ! ================================================        
          ! Litter fluxes for understorey  mortality. KgC/m2/year
          ! ================================================

          ! Total number of dead understory (n/m2)
          dead_n = -1.0_r8 * currentCohort%dndt / currentPatch%area

          ! Total number of dead understory from direct logging
          ! (it is possible that large harvestable trees are in the understory)
          dead_n_dlogging = ( currentCohort%lmort_direct) * &
               currentCohort%n/hlm_freq_day/currentPatch%area
          
          ! Total number of dead understory from indirect logging
          dead_n_ilogging = ( currentCohort%lmort_collateral + currentCohort%lmort_infra) * &
                currentCohort%n/hlm_freq_day/currentPatch%area
          
          dead_n_natural = dead_n - dead_n_dlogging - dead_n_ilogging

          
          currentPatch%leaf_litter_in(pft) = currentPatch%leaf_litter_in(pft) + &
               (leaf_c)* dead_n
                ! %n has not been updated due to mortality yet, thus
                ! the litter flux has already been counted since it captured
                ! the losses of live trees and those flagged for death

          currentPatch%root_litter_in(pft) = currentPatch%root_litter_in(pft) + &
               (fnrt_c + store_c*(1._r8-EDPftvarcon_inst%allom_frbstor_repro(pft)) ) * dead_n

          ! Update diagnostics that track resource management
          currentSite%resources_management%delta_litter_stock  = &
                currentSite%resources_management%delta_litter_stock + &
                (leaf_c + fnrt_c + store_c ) * &
                (dead_n_ilogging+dead_n_dlogging) * & 
                hlm_freq_day * currentPatch%area

          ! Update diagnostics that track resource management
          currentSite%resources_management%delta_biomass_stock = &
                currentSite%resources_management%delta_biomass_stock + &
                (leaf_c + fnrt_c + store_c ) * & 
                (dead_n_ilogging+dead_n_dlogging) * & 
                hlm_freq_day * currentPatch%area

          if( hlm_use_planthydro == itrue ) then
             !call AccumulateMortalityWaterStorage(currentSite,currentCohort,dead_n)
             call AccumulateMortalityWaterStorage(currentSite,currentCohort,&
                                                  -1.0_r8 * currentCohort%dndt * hlm_freq_day)
          end if
          

          do c = 1,ncwd
             
             currentPatch%cwd_BG_in(c) = currentPatch%cwd_BG_in(c) + (struct_c + sapw_c) * & 
                   SF_val_CWD_frac(c) * dead_n * (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft))

             ! Send AGB component of boles from non direct-logging activities to AGB litter pool
             if (c==ncwd) then
                
                currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + (struct_c + sapw_c) * & 
                     SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging)  * &
                     EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)
                
             else

                currentPatch%cwd_AG_in(c) = currentPatch%cwd_AG_in(c) + (struct_c + sapw_c) * & 
                     SF_val_CWD_frac(c) * dead_n  * &
                     EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)

                ! Send AGB component of boles from direct-logging activities to export/harvest pool
                ! Generate trunk product (kgC/day/site)
                trunk_product =  (struct_c + sapw_c) * &
                      SF_val_CWD_frac(c) * dead_n_dlogging * EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * &
                      hlm_freq_day * currentPatch%area
                
                currentSite%flux_out = currentSite%flux_out + trunk_product

                ! Update diagnostics that track resource management
                currentSite%resources_management%trunk_product_site  = &
                      currentSite%resources_management%trunk_product_site + &
                      trunk_product
                ! Update diagnostics that track resource management
                currentSite%resources_management%trunk_product_site  = &
                      currentSite%resources_management%trunk_product_site + &
                      trunk_product
             end if

             ! Update diagnostics that track resource management
             currentSite%resources_management%delta_litter_stock  = &
                   currentSite%resources_management%delta_litter_stock + &
                   (struct_c + sapw_c) * &
                   SF_val_CWD_frac(c) * (dead_n_natural+dead_n_ilogging) * & 
                   hlm_freq_day * currentPatch%area
             ! Update diagnostics that track resource management
             currentSite%resources_management%delta_biomass_stock = &
                   currentSite%resources_management%delta_biomass_stock + &
                   (struct_c + sapw_c) * &
                   SF_val_CWD_frac(c) * dead_n * & 
                   hlm_freq_day * currentPatch%area
             
             if (currentPatch%cwd_AG_in(c) < 0.0_r8)then
                write(fates_log(),*) 'negative CWD in flux',currentPatch%cwd_AG_in(c), &
                      (struct_c + sapw_c), dead_n
             endif

          end do
          ! Update diagnostics that track resource management
          currentSite%resources_management%delta_individual    = &
                currentSite%resources_management%delta_individual + &
                (dead_n_dlogging+dead_n_ilogging) * hlm_freq_day * currentPatch%area
          
       !endif !canopy layer
       
       currentCohort => currentCohort%taller
    enddo  ! end loop over cohorts 

    do p = 1,numpft
       currentPatch%leaf_litter_in(p) = currentPatch%leaf_litter_in(p) + currentPatch%seed_decay(p) !KgC/m2/yr
    enddo

  end subroutine CWD_Input

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
  subroutine cwd_out( currentSite, currentPatch, bc_in )
    !
    ! !DESCRIPTION:
    ! Simple CWD fragmentation Model
    ! spawn new cohorts of juveniles of each PFT             
    !
    ! !USES:
    use SFParamsMod, only : SF_val_max_decomp

    !
    ! !ARGUMENTS    
    type(ed_site_type), intent(inout), target  :: currentSite
    type(ed_patch_type), intent(inout), target :: currentPatch
    type(bc_in_type), intent(in)               :: bc_in
    
    !
    ! !LOCAL VARIABLES:
    integer :: c,ft
    !----------------------------------------------------------------------

    currentPatch%root_litter_out(:) = 0.0_r8
    currentPatch%leaf_litter_out(:) = 0.0_r8
    
    call fragmentation_scaler(currentPatch, bc_in)

    !Flux of coarse woody debris into decomposing litter pool. 

    currentPatch%cwd_ag_out(1:ncwd) = 0.0_r8
    currentPatch%cwd_bg_out(1:ncwd) = 0.0_r8
    currentPatch%leaf_litter_out(:) = 0.0_r8
    currentPatch%root_litter_out(:) = 0.0_r8
    
    do c = 1,ncwd  
       currentPatch%cwd_ag_out(c)      = max(0.0_r8,   currentPatch%cwd_ag(c) * &
            SF_val_max_decomp(c+1) * currentPatch%fragmentation_scaler )  
       currentPatch%cwd_bg_out(c)      = max(0.0_r8,   currentPatch%cwd_bg(c) * &
            SF_val_max_decomp(c+1) * currentPatch%fragmentation_scaler )
    enddo

    ! this is the rate at which dropped leaves stop being part of the burnable pool and begin to be part of the 
    ! decomposing pool. This should probably be highly sensitive to moisture, but also to the type of leaf 
    ! thick leaves can dry out before they are decomposed, for example. 
    ! this section needs further scientific input. 

    do ft = 1,numpft
       currentPatch%leaf_litter_out(ft) = max(0.0_r8,currentPatch%leaf_litter(ft)* SF_val_max_decomp(dl_sf) * &
            currentPatch%fragmentation_scaler )
       currentPatch%root_litter_out(ft) = max(0.0_r8,currentPatch%root_litter(ft)* SF_val_max_decomp(dl_sf) * &
            currentPatch%fragmentation_scaler )
       if ( currentPatch%leaf_litter_out(ft)<0.0_r8.or.currentPatch%root_litter_out(ft)<0.0_r8)then
         write(fates_log(),*) 'root or leaf out is negative?',SF_val_max_decomp(dl_sf),currentPatch%fragmentation_scaler
       endif
    enddo

    !add up carbon going into fragmenting pools
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%leaf_litter_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%root_litter_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_ag_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day
    currentSite%flux_out = currentSite%flux_out + sum(currentPatch%cwd_bg_out) * &
         currentPatch%area *hlm_freq_day!kgC/site/day

  end subroutine cwd_out

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
    use EDPftvarcon, only : EDPftvarcon_inst
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
    type(bc_out_type)       , intent(inout)           :: bc_out(:)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort
    type(ed_site_type), pointer :: cs
    integer p,ci,j,s
    real(r8) time_convert    ! from year to seconds
    real(r8) mass_convert    ! ED uses kg, CLM uses g
    integer           :: begp,endp
    integer           :: begc,endc                                    !bounds 

    !------------------------------------------------------------------------
    ! The following scratch arrays are allocated for maximum possible
    ! pft and layer usage
    
    real(r8) :: cinput_rootfr(1:maxpft, 1:hlm_numlevgrnd)
    real(r8) :: croot_prof_perpatch(1:hlm_numlevgrnd)
    real(r8) :: surface_prof(1:hlm_numlevgrnd)
    integer  :: ft
    integer  :: nlev_eff_decomp
    real(r8) :: rootfr_tot(1:maxpft)
    real(r8) :: biomass_bg_ft(1:maxpft)
    real(r8) :: surface_prof_tot, leaf_prof_sum, stem_prof_sum, froot_prof_sum, biomass_bg_tot
    real(r8) :: delta
    real(r8) :: leaf_c
    real(r8) :: store_c
    real(r8) :: fnrt_c
    real(r8) :: sapw_c
    real(r8) :: struct_c

    ! NOTE(rgk, 201705) this parameter was brought over from SoilBiogeochemVerticalProfile
    ! how steep profile is for surface components (1/ e_folding depth) (1/m) 
    real(r8),  parameter :: surfprof_exp  = 10.

    real(r8) :: leaf_prof(1:nsites, 1:hlm_numlevgrnd)
    real(r8) :: froot_prof(1:nsites,  1:maxpft, 1:hlm_numlevgrnd)
    real(r8) :: croot_prof(1:nsites, 1:hlm_numlevgrnd)
    real(r8) :: stem_prof(1:nsites, 1:hlm_numlevgrnd)

    

    delta = 0.001_r8    
    !no of seconds in a year. 
    time_convert =  365.0_r8*sec_per_day

    ! number of grams in a kilogram
    mass_convert = 1000._r8
    
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! first calculate vertical profiles
    ! define two types of profiles: 
    ! (1) a surface profile, for leaves and stem inputs, which is the same for each
    ! pft but differs from one site to the next to avoid inputting any C into permafrost or bedrock
    ! (2) a fine root profile, which is indexed by both site and pft, differs for 
    ! each pft and also from one site to the next to avoid inputting any C into permafrost or bedrock
    ! (3) a coarse root profile, which is the root-biomass=weighted average of the fine root profiles
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    if (hlm_use_vertsoilc == itrue) then

       ! initialize profiles to zero
       leaf_prof(1:nsites, 1:hlm_numlevgrnd )               = 0._r8
       froot_prof(1:nsites, 1:maxpft, 1:hlm_numlevgrnd)     = 0._r8
       stem_prof(1:nsites, 1:hlm_numlevgrnd)                = 0._r8
       
       do s = 1,nsites

          ! Calculate the number of effective decomposition layers
          ! This takes into account if vertical soil biogeochem is on, how deep the soil column
          ! is, and also which layers may be frozen
          nlev_eff_decomp = min(max(bc_in(s)%max_rooting_depth_index_col, 1), bc_in(s)%nlevdecomp)

          ! define a single shallow surface profile for surface additions (leaves, stems, and N deposition)
          surface_prof(:) = 0._r8
          do j = 1,  bc_in(s)%nlevdecomp
             surface_prof(j) = exp(-surfprof_exp * bc_in(s)%z_sisl(j)) / bc_in(s)%dz_decomp_sisl(j)
          end do
          
          ! -----------------------------------------------------------------------------
          ! This is the rooting profile.  cinput_rootfr
          ! This array will calculate root
          ! mass as far down as the soil column goes.  It is possible
          ! that the active layers are not as deep as the roots go.
          ! That is ok, the roots in the active layers will be talied up and
          ! normalized.
          ! -----------------------------------------------------------------------------

          cinput_rootfr(:,:)     = 0._r8
          do ft = 1, numpft
             
             ! This generates a rooting profile over the whole soil column for each pft
             ! Note that we are calling for the root fractions in the biomass
             ! for litter context, and not the hydrologic uptake context.

             call set_root_fraction(cinput_rootfr(ft,1:bc_in(s)%nlevsoil), ft, &
                  bc_in(s)%zi_sisl, lowerb=lbound(bc_in(s)%zi_sisl,1), &
                  icontext=i_biomass_rootprof_context)
             
             do j=1,nlev_eff_decomp
                cinput_rootfr(ft,j) = cinput_rootfr(ft,j)/bc_in(s)%dz_decomp_sisl(j)
             end do

          end do

          !
          ! now add permafrost constraint: integrate rootfr over active layer of soil site,
          ! truncate below permafrost or bedrock table where present, and rescale so that integral = 1
          rootfr_tot(:) = 0._r8
          
          surface_prof_tot = 0._r8
          !
          do j = 1, nlev_eff_decomp
             surface_prof_tot = surface_prof_tot + surface_prof(j)  * bc_in(s)%dz_decomp_sisl(j)
          end do

          do ft = 1,numpft
             do j = 1, nlev_eff_decomp
                rootfr_tot(ft) = rootfr_tot(ft) + cinput_rootfr(ft,j)*bc_in(s)%dz_decomp_sisl(j)
             end do
          end do
          !
          ! rescale the fine root profile
          do ft = 1,numpft
             if ( (bc_in(s)%max_rooting_depth_index_col > 0) .and. (rootfr_tot(ft) > 0._r8) ) then
                ! where there is not permafrost extending to the surface, integrate the profiles 
                ! over the active layer this is equivalent to integrating over all soil layers 
                ! outside of permafrost regions
                do j = 1, nlev_eff_decomp
                   froot_prof(s,ft,j) = cinput_rootfr(ft,j) / rootfr_tot(ft)
                end do
             else
                ! if fully frozen, or no roots, put everything in the top layer
                froot_prof(s,ft,1) = 1._r8/bc_in(s)%dz_decomp_sisl(1)
             endif
          end do

          !
          ! rescale the shallow profiles
          if ( (bc_in(s)%max_rooting_depth_index_col > 0) .and. (surface_prof_tot > 0._r8) ) then
             ! where there is not permafrost extending to the surface, integrate the profiles over 
             ! the active layer this is equivalent to integrating over all soil layers outside of 
             ! permafrost regions
             do j = 1, nlev_eff_decomp
                ! set all surface processes to shallower profile
                leaf_prof(s,j) = surface_prof(j)/ surface_prof_tot
                stem_prof(s,j) = surface_prof(j)/ surface_prof_tot
             end do
          else
             ! if fully frozen, or no roots, put everything in the top layer
             leaf_prof(s,1) = 1._r8/bc_in(s)%dz_decomp_sisl(1)
             stem_prof(s,1) = 1._r8/bc_in(s)%dz_decomp_sisl(1)
             do j = 2, bc_in(s)%nlevdecomp
                leaf_prof(s,j) = 0._r8
                stem_prof(s,j) = 0._r8
             end do
          endif
       end do
       
    else
       
       ! for one layer decomposition model, set profiles to unity
       leaf_prof(1:nsites, :) = 1._r8
       froot_prof(1:nsites, 1:numpft, :) = 1._r8
       stem_prof(1:nsites, :) = 1._r8
       
    end if
    
    ! sanity check to ensure they integrate to 1
    do s = 1, nsites
       ! check the leaf and stem profiles
       leaf_prof_sum = 0._r8
       stem_prof_sum = 0._r8
       do j = 1, bc_in(s)%nlevdecomp
          leaf_prof_sum = leaf_prof_sum + leaf_prof(s,j) *  bc_in(s)%dz_decomp_sisl(j)
          stem_prof_sum = stem_prof_sum + stem_prof(s,j) *  bc_in(s)%dz_decomp_sisl(j)
       end do
       if ( ( abs(stem_prof_sum - 1._r8) > delta ) .or.  ( abs(leaf_prof_sum - 1._r8) > delta ) ) then
          write(fates_log(), *) 'profile sums: ',  leaf_prof_sum, stem_prof_sum
          write(fates_log(), *) 'surface_prof: ', surface_prof
          write(fates_log(), *) 'surface_prof_tot: ', surface_prof_tot
          write(fates_log(), *) 'leaf_prof: ',  leaf_prof(s,:)
          write(fates_log(), *) 'stem_prof: ',  stem_prof(s,:)
          write(fates_log(), *) 'max_rooting_depth_index_col: ', bc_in(s)%max_rooting_depth_index_col
          write(fates_log(), *) 'bc_in(s)%dz_decomp_sisl: ',  bc_in(s)%dz_decomp_sisl            
          call endrun(msg=errMsg(sourcefile, __LINE__))
       endif
       ! now check each fine root profile
       do ft = 1,numpft 
          froot_prof_sum = 0._r8
          do j = 1, bc_in(s)%nlevdecomp
             froot_prof_sum = froot_prof_sum + froot_prof(s,ft,j) *  bc_in(s)%dz_decomp_sisl(j)
          end do
          if ( ( abs(froot_prof_sum - 1._r8) > delta ) ) then
             write(fates_log(), *) 'profile sums: ', froot_prof_sum
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
       end do
    end do
    
    ! zero the site-level C input variables
    do s = 1, nsites
       do j = 1, bc_in(s)%nlevdecomp
          bc_out(s)%FATES_c_to_litr_lab_c_col(j) = 0._r8
          bc_out(s)%FATES_c_to_litr_cel_c_col(j) = 0._r8
          bc_out(s)%FATES_c_to_litr_lig_c_col(j) = 0._r8
          croot_prof(s,j)         = 0._r8
       end do
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! now disaggregate the inputs vertically, using the vertical profiles
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do s = 1,nsites
         

         currentPatch => sites(s)%oldest_patch
         do while(associated(currentPatch))
            
            ! the CWD pools lose information about which PFT they came from; 
            ! for the stems this doesn't matter as they all have the same profile, 
            ! however for the coarse roots they may have different profiles.  
            ! to approximately recover this information, loop over all cohorts in patch 
            ! to calculate the total root biomass in that patch of each pft, and then 
            ! rescale the croot_prof as the weighted average of the froot_prof
            biomass_bg_ft(1:numpft) = 0._r8
            currentCohort => currentPatch%tallest
            do while(associated(currentCohort))      
               
               leaf_c          = currentCohort%prt%GetState(leaf_organ,all_carbon_elements)
               store_c         = currentCohort%prt%GetState(store_organ,all_carbon_elements)
               fnrt_c          = currentCohort%prt%GetState(fnrt_organ,all_carbon_elements)
               sapw_c          = currentCohort%prt%GetState(sapw_organ,all_carbon_elements)
               struct_c        = currentCohort%prt%GetState(struct_organ,all_carbon_elements)

               biomass_bg_ft(currentCohort%pft) = biomass_bg_ft(currentCohort%pft) + &
                    ( (struct_c + sapw_c) * & 
                    (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) + &
                    (fnrt_c + store_c ) ) * & 
                    (currentCohort%n / currentPatch%area)

               currentCohort => currentCohort%shorter
            enddo !currentCohort
            ! 
            biomass_bg_tot = 0._r8
            do ft = 1,numpft 
               biomass_bg_tot = biomass_bg_tot + biomass_bg_ft(ft)
            end do
            !         

            ! zero this for each patch
            croot_prof_perpatch(1:bc_in(s)%nlevdecomp) = 0._r8

            !
            if ( biomass_bg_tot .gt. 0._r8) then
               do ft = 1,numpft 
                  do j = 1, bc_in(s)%nlevdecomp
                     croot_prof_perpatch(j) = croot_prof_perpatch(j) + &
                          froot_prof(s,ft,j) * biomass_bg_ft(ft) / biomass_bg_tot
                  end do
               end do
            else ! no biomass
               croot_prof_perpatch(1) = 1./bc_in(s)%dz_decomp_sisl(1)
            end if

            !
            ! add croot_prof as weighted average (weighted by patch area) of croot_prof_perpatch
            do j = 1, bc_in(s)%nlevdecomp
               croot_prof(s, j) = croot_prof(s, j) + croot_prof_perpatch(j) * currentPatch%area / AREA
            end do
            !
            ! now disaggregate, vertically and by decomposition substrate type, the 
            ! actual fluxes from CWD and litter pools
            !
            ! do c = 1, ncwd
            !    write(fates_log(),*)'cdk CWD_AG_out', c, currentpatch%CWD_AG_out(c), 
            !                         ED_val_cwd_fcel, currentpatch%area/AREA
            !    write(fates_log(),*)'cdk CWD_BG_out', c, currentpatch%CWD_BG_out(c), 
            !                         ED_val_cwd_fcel, currentpatch%area/AREA
            ! end do
            ! do ft = 1,numpft
            !    write(fates_log(),*)'cdk leaf_litter_out', ft, currentpatch%leaf_litter_out(ft), 
            !                         ED_val_cwd_fcel, currentpatch%area/AREA
            !    write(fates_log(),*)'cdk root_litter_out', ft, currentpatch%root_litter_out(ft), 
            !                         ED_val_cwd_fcel, currentpatch%area/AREA
            ! end do
            ! !
            ! CWD pools fragmenting into decomposing litter pools. 
            do ci = 1, ncwd
               do j = 1, bc_in(s)%nlevdecomp
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%CWD_AG_out(ci) * ED_val_cwd_fcel * &
                       currentpatch%area/AREA * stem_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%CWD_AG_out(ci) * ED_val_cwd_flig * &
                       currentpatch%area/AREA * stem_prof(s,j)
                  !
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%CWD_BG_out(ci) * ED_val_cwd_fcel * &
                       currentpatch%area/AREA * croot_prof_perpatch(j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%CWD_BG_out(ci) * ED_val_cwd_flig * &
                       currentpatch%area/AREA * croot_prof_perpatch(j)
               end do
            end do
            
            ! leaf and fine root pools. 
            do ft = 1,numpft
               do j = 1, bc_in(s)%nlevdecomp
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * EDPftvarcon_inst%lf_flab(ft) * &
                       currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * EDPftvarcon_inst%lf_fcel(ft) * &
                       currentpatch%area/AREA * leaf_prof(s,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%leaf_litter_out(ft) * EDPftvarcon_inst%lf_flig(ft) * &
                       currentpatch%area/AREA * leaf_prof(s,j)
                  !
                  bc_out(s)%FATES_c_to_litr_lab_c_col(j) = bc_out(s)%FATES_c_to_litr_lab_c_col(j) + &
                       currentpatch%root_litter_out(ft) * EDPftvarcon_inst%fr_flab(ft) * &
                       currentpatch%area/AREA * froot_prof(s,ft,j)
                  bc_out(s)%FATES_c_to_litr_cel_c_col(j) = bc_out(s)%FATES_c_to_litr_cel_c_col(j) + &
                       currentpatch%root_litter_out(ft) * EDPftvarcon_inst%fr_fcel(ft) * &
                       currentpatch%area/AREA * froot_prof(s,ft,j)
                  bc_out(s)%FATES_c_to_litr_lig_c_col(j) = bc_out(s)%FATES_c_to_litr_lig_c_col(j) + &
                       currentpatch%root_litter_out(ft) * EDPftvarcon_inst%fr_flig(ft) * &
                       currentpatch%area/AREA * froot_prof(s,ft,j)
               enddo
            end do
              
              currentPatch => currentPatch%younger
           end do !currentPatch

        end do  ! do sites(s)
     
        do s = 1, nsites
           do j = 1, bc_in(s)%nlevdecomp
              ! time unit conversion
              bc_out(s)%FATES_c_to_litr_lab_c_col(j)=bc_out(s)%FATES_c_to_litr_lab_c_col(j) * &
                   mass_convert / time_convert
              bc_out(s)%FATES_c_to_litr_cel_c_col(j)=bc_out(s)%FATES_c_to_litr_cel_c_col(j) * &
                   mass_convert / time_convert
              bc_out(s)%FATES_c_to_litr_lig_c_col(j)=bc_out(s)%FATES_c_to_litr_lig_c_col(j) * &
                   mass_convert / time_convert
           end do
        end do
        
        ! write(fates_log(),*)'cdk FATES_c_to_litr_lab_c: ', FATES_c_to_litr_lab_c
        ! write_col(fates_log(),*)'cdk FATES_c_to_litr_cel_c: ', FATES_c_to_litr_cel_c    
        ! write_col(fates_log(),*)'cdk FATES_c_to_litr_lig_c: ', FATES_c_to_litr_lig_c
        ! write_col(fates_log(),*)'cdk bounds%begc, bounds%endc: ', bounds%begc, bounds%endc
        ! write(fates_log(),*)'cdk leaf_prof: ', leaf_prof
        ! write(fates_log(),*)'cdk stem_prof: ', stem_prof    
        ! write(fates_log(),*)'cdk froot_prof: ', froot_prof
        ! write(fates_log(),*)'cdk croot_prof_perpatch: ', croot_prof_perpatch
        ! write(fates_log(),*)'cdk croot_prof: ', croot_prof

    end subroutine flux_into_litter_pools



end module EDPhysiologyMod
