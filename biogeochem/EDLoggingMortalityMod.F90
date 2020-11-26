module EDLoggingMortalityMod

   ! ====================================================================================
   !  Purpose: 1. create logging mortalities: 
   !           (a) direct logging mortality (cohort level)
   !           (b) collateral mortality (cohort level)
   !           (c) infrastructure mortality (cohort level)
   !           2. move the logged trunk fluxes from live into product pool 
   !           3. move logging-associated mortality fluxes from live to CWD
   !           4. keep carbon balance (in ed_total_balance_check)
   !
   !  Yi Xu & M.Huang
   !  Date: 09/2017
   ! ====================================================================================

   use FatesConstantsMod , only : r8 => fates_r8
   use EDTypesMod        , only : ed_cohort_type
   use EDTypesMod        , only : ed_patch_type
   use EDTypesMod        , only : site_massbal_type
   use EDTypesMod        , only : site_fluxdiags_type
   use FatesLitterMod    , only : ncwd
   use FatesLitterMod    , only : ndcmpy
   use FatesLitterMod    , only : litter_type
   use EDTypesMod        , only : ed_site_type
   use EDTypesMod        , only : ed_resources_management_type
   use EDTypesMod        , only : dtype_ilog
   use EDTypesMod        , only : dtype_ifall
   use EDTypesMod        , only : dtype_ifire
   use EDPftvarcon       , only : EDPftvarcon_inst
   use EDPftvarcon       , only : GetDecompyFrac
   use PRTParametersMod  , only : prt_params
   use PRTGenericMod     , only : num_elements
   use PRTGenericMod     , only : element_list
   use EDParamsMod       , only : logging_export_frac
   use EDParamsMod       , only : logging_event_code
   use EDParamsMod       , only : logging_dbhmin
   use EDParamsMod       , only : logging_dbhmax
   use EDParamsMod       , only : logging_collateral_frac 
   use EDParamsMod       , only : logging_direct_frac
   use EDParamsMod       , only : logging_mechanical_frac 
   use EDParamsMod       , only : logging_coll_under_frac 
   use EDParamsMod       , only : logging_dbhmax_infra
   use FatesInterfaceTypesMod , only : hlm_current_year
   use FatesInterfaceTypesMod , only : hlm_current_month
   use FatesInterfaceTypesMod , only : hlm_current_day
   use FatesInterfaceTypesMod , only : hlm_model_day
   use FatesInterfaceTypesMod , only : hlm_day_of_year 
   use FatesInterfaceTypesMod , only : hlm_days_per_year
   use FatesInterfaceTypesMod , only : hlm_use_lu_harvest
   use FatesInterfaceTypesMod , only : hlm_num_lu_harvest_cats
   use FatesInterfaceTypesMod , only : hlm_use_logging 
   use FatesInterfaceTypesMod , only : hlm_use_planthydro
   use FatesConstantsMod , only : itrue,ifalse
   use FatesGlobals      , only : endrun => fates_endrun 
   use FatesGlobals      , only : fates_log
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use FatesPlantHydraulicsMod, only : AccumulateMortalityWaterStorage
   use PRTGenericMod     , only : all_carbon_elements,carbon12_element
   use PRTGenericMod     , only : sapw_organ, struct_organ, leaf_organ
   use PRTGenericMod     , only : fnrt_organ, store_organ, repro_organ
   use FatesAllometryMod , only : set_root_fraction
   use FatesConstantsMod , only : primaryforest, secondaryforest, secondary_age_threshold
   use FatesConstantsMod , only : fates_tiny
   use FatesConstantsMod , only : months_per_year
   use FatesConstantsMod , only : hlm_harvest_area_fraction
   use FatesConstantsMod , only : hlm_harvest_carbon
   use FatesConstantsMod, only : fates_check_param_set

   implicit none
   private

   logical, protected :: logging_time   ! If true, logging should be 
                                        ! performed during the current time-step


   ! harvest litter localization specifies how much of the litter from a falling
   ! tree lands within the newly generated patch, and how much lands outside of 
   ! the new patch, and thus in the original patch.  By setting this to zero,
   ! it is assumed that there is no preference, and thus the mass is distributed
   ! equally.  If this is set to 1, then all of the mass lands in the new
   ! patch, and is thus "completely local".


   real(r8), parameter :: harvest_litter_localization = 0.0_r8

   character(len=*), parameter, private :: sourcefile = &
         __FILE__
   
   public :: LoggingMortality_frac
   public :: logging_litter_fluxes
   public :: logging_time
   public :: IsItLoggingTime
   public :: get_harvest_rate_area

contains

   subroutine IsItLoggingTime(is_master,currentSite)

      ! -------------------------------------------------------------------------------
      ! This subroutine determines if the current dynamics step should enact
      ! the logging module.
      ! This is done by comparing the current model time to the logging event
      ! ids.  If there is a match, it is logging time.
      ! -------------------------------------------------------------------------------
     
      integer, intent(in) :: is_master
      type(ed_site_type), intent(inout), target :: currentSite     ! site structure

      integer :: icode   ! Integer equivalent of the event code (parameter file only allows reals)
      integer :: log_date  ! Day of month for logging exctracted from event code
      integer :: log_month ! Month of year for logging extraced from event code
      integer :: log_year  ! Year for logging extracted from event code
      character(len=64) :: fmt = '(a,i2.2,a,i2.2,a,i4.4)'

      logging_time = .false.
      icode = int(logging_event_code)

      ! this is true for either hlm harvest or fates logging
      if(hlm_use_logging.eq.ifalse) return

      ! all of these are valid for hlm harvest
      ! so adjust annual harvest inputs accordingly in LoggingMortality_frac
      ! note that the specific event will allow only one hlm harvest, regardless of input
      ! code 3, every day, may not work properly because of the exclusive mortality selection
      ! even less fequent low harvest rates may be excluded - may need to give harvest priority

      if(icode .eq. 1) then
         ! Logging is turned off - not sure why we need another switch
         logging_time = .false.

      else if(icode .eq. 2) then
         ! Logging event on the first step
         if( hlm_model_day.eq.1 ) then
            logging_time = .true.
         end if

      else if(icode .eq. 3) then
         ! Logging event every day
         logging_time = .true.

      else if(icode .eq. 4) then
         ! logging event once a month
         if(hlm_current_day.eq.1  ) then
            logging_time = .true.
         end if

      else if(icode < 0 .and. icode > -366) then
         ! Logging event every year on specific day of year
         if(hlm_day_of_year .eq. abs(icode)  ) then
            logging_time = .true.
         end if

      else if(icode > 10000 ) then
         ! Specific Event: YYYYMMDD
         log_date  = icode - int(100* floor(real(icode)/100))
         log_year  = floor(real(icode)/10000)
         log_month = floor(real(icode)/100) - log_year*100

         if( hlm_current_day.eq.log_date    .and. &
               hlm_current_month.eq.log_month .and. &
               hlm_current_year.eq.log_year ) then
            logging_time = .true.
         end if
      else 
         ! Bad logging event flag
         write(fates_log(),*) 'An invalid logging code was specified in fates_params'
         write(fates_log(),*) 'Check EDLoggingMortalityMod.F90:IsItLoggingTime()'
         write(fates_log(),*) 'for a breakdown of the valid codes and change'
         write(fates_log(),*) 'fates_logging_event_code in the file accordingly.'
         write(fates_log(),*) 'exiting'
         call endrun(msg=errMsg(sourcefile, __LINE__))
      end if

      ! Initialize some site level diagnostics that are calculated for each event
      currentSite%resources_management%delta_litter_stock  = 0.0_r8
      currentSite%resources_management%delta_biomass_stock = 0.0_r8
      currentSite%resources_management%delta_individual    = 0.0_r8

      if(logging_time .and. (is_master.eq.itrue) ) then
         write(fates_log(),fmt) 'Logging Event Enacted on date: ', &
               hlm_current_month,'-',hlm_current_day,'-',hlm_current_year
      end if
      return
   end subroutine IsItLoggingTime

   ! ======================================================================================

   subroutine LoggingMortality_frac( pft_i, dbh, canopy_layer, lmort_direct, &
                                     lmort_collateral,lmort_infra, l_degrad, &
                                     hlm_harvest_rates, hlm_harvest_catnames, &
                                     hlm_harvest_units, &
                                     patch_anthro_disturbance_label, secondary_age, &
                                     frac_site_primary)

      ! Arguments
      integer,  intent(in)  :: pft_i            ! pft index 
      real(r8), intent(in)  :: dbh              ! diameter at breast height (cm)
      integer,  intent(in)  :: canopy_layer     ! canopy layer of this cohort
      real(r8), intent(in) :: hlm_harvest_rates(:) ! annual harvest rate per hlm category
      character(len=64), intent(in) :: hlm_harvest_catnames(:) ! names of hlm harvest categories
      integer, intent(in) :: hlm_harvest_units     ! unit type of hlm harvest rates: [area vs. mass]
      integer, intent(in) :: patch_anthro_disturbance_label    ! patch level anthro_disturbance_label
      real(r8), intent(in) :: secondary_age     ! patch level age_since_anthro_disturbance
      real(r8), intent(out) :: lmort_direct     ! direct (harvestable) mortality fraction
      real(r8), intent(out) :: lmort_collateral ! collateral damage mortality fraction
      real(r8), intent(out) :: lmort_infra      ! infrastructure mortality fraction
      real(r8), intent(out) :: l_degrad         ! fraction of trees that are not killed
                                                ! but suffer from forest degradation (i.e. they
                                                ! are moved to newly-anthro-disturbed secondary
                                                ! forest patch)
      real(r8), intent(in) :: frac_site_primary

      ! Local variables
      real(r8) :: harvest_rate ! the final harvest rate to apply to this cohort today

      ! todo: probably lower the dbhmin default value to 30 cm
      ! todo: change the default logging_event_code to 1 september (-244)
      ! todo: change the default logging_direct_frac to 1.0 for cmip inputs
      ! todo: check outputs against the LUH2 carbon data
      ! todo: eventually set up distinct harvest practices, each with a set of input paramaeters
      ! todo: implement harvested carbon inputs
      
      if (logging_time) then 

         ! Pass logging rates to cohort level 
         
         if (hlm_use_lu_harvest == ifalse) then
            ! 0=use fates logging parameters directly when logging_time == .true.
            ! this means harvest the whole cohort area
            harvest_rate = 1._r8
            
         else if (hlm_use_lu_harvest == itrue .and. hlm_harvest_units == hlm_harvest_area_fraction) then
            ! We are harvesting based on areal fraction, not carbon/biomass terms. 
            ! 1=use area fraction from hlm
            ! combine forest and non-forest fracs and then apply:
            ! primary and secondary area fractions to the logging rates, which are fates parameters
            
            ! Definitions of the underlying harvest land category variables
            ! these are hardcoded to match the LUH input data via landuse.timseries file (see dynHarvestMod)
            ! these are fractions of vegetated area harvested, split into five land category variables
            ! HARVEST_VH1 = harvest from primary forest
            ! HARVEST_VH2 = harvest from primary non-forest
            ! HARVEST_SH1 = harvest from secondary mature forest
            ! HARVEST_SH2 = harvest from secondary young forest
            ! HARVEST_SH3 = harvest from secondary non-forest (assume this is young for biomass)

            ! Get the area-based harvest rates based on info passed to FATES from the boundary condition
            call get_harvest_rate_area (patch_anthro_disturbance_label, hlm_harvest_catnames, &
                 hlm_harvest_rates, frac_site_primary, secondary_age, harvest_rate)

         else if (hlm_use_lu_harvest == itrue .and. hlm_harvest_units == hlm_harvest_carbon) then
            ! 2=use carbon from hlm
            ! not implemented yet
            write(fates_log(),*) 'HLM harvest carbon data not implemented yet. Exiting.'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif

         ! transfer of area to secondary land is based on overall area affected, not just logged crown area
         ! l_degrad accounts for the affected area between logged crowns
         if(int(prt_params%woody(pft_i)) == 1)then ! only set logging rates for trees
            
            ! direct logging rates, based on dbh min and max criteria
            if (dbh >= logging_dbhmin .and. .not. &
                 ((logging_dbhmax < fates_check_param_set) .and. (dbh >= logging_dbhmax )) ) then
               ! the logic of the above line is a bit unintuitive but allows turning off the dbhmax comparison entirely.
               ! since there is an .and. .not. after the first conditional, the dbh:dbhmax comparison needs to be 
               ! the opposite of what would otherwise be expected...
               lmort_direct = harvest_rate * logging_direct_frac

            else
               lmort_direct = 0.0_r8
            end if

            ! infrastructure (roads, skid trails, etc) mortality rates
            if (dbh >= logging_dbhmax_infra) then
               lmort_infra      = 0.0_r8
            else
               lmort_infra      = harvest_rate * logging_mechanical_frac
            end if

            ! Collateral damage to smaller plants below the direct logging size threshold
            ! will be applied via "understory_death" via the disturbance algorithm
            if (canopy_layer .eq. 1) then
               lmort_collateral = harvest_rate * logging_collateral_frac
            else
               lmort_collateral = 0._r8
            endif

         else  ! non-woody plants still killed by infrastructure
            lmort_direct    = 0.0_r8
            lmort_collateral = 0.0_r8
            lmort_infra      = harvest_rate * logging_mechanical_frac
         end if

         ! the area occupied by all plants in the canopy that aren't killed is still disturbed at the harvest rate
         if (canopy_layer .eq. 1) then
            l_degrad = harvest_rate - (lmort_direct + lmort_infra + lmort_collateral) ! fraction passed to 'degraded' forest.
         else
            l_degrad = 0._r8
         endif
         
      else 
         lmort_direct    = 0.0_r8
         lmort_collateral = 0.0_r8
         lmort_infra      = 0.0_r8
         l_degrad         = 0.0_r8
      end if

   end subroutine LoggingMortality_frac

   ! ============================================================================

   subroutine get_harvest_rate_area (patch_anthro_disturbance_label, hlm_harvest_catnames, hlm_harvest_rates, &
                 frac_site_primary, secondary_age, harvest_rate)


     ! -------------------------------------------------------------------------------------------
     !
     !  DESCRIPTION:
     !  get the area-based harvest rates based on info passed to FATES from the bioundary conditions in.
     !  assumes logging_time == true

      ! Arguments
      real(r8), intent(in) :: hlm_harvest_rates(:) ! annual harvest rate per hlm category
      character(len=64), intent(in) :: hlm_harvest_catnames(:) ! names of hlm harvest categories
      integer, intent(in) :: patch_anthro_disturbance_label    ! patch level anthro_disturbance_label
      real(r8), intent(in) :: secondary_age     ! patch level age_since_anthro_disturbance
      real(r8), intent(in) :: frac_site_primary
      real(r8), intent(out) :: harvest_rate

      ! Local Variables
      integer :: h_index   ! for looping over harvest categories
      integer :: icode   ! Integer equivalent of the event code (parameter file only allows reals)

     !  Loop around harvest categories to determine the annual hlm harvest rate for the current cohort based on patch history info
     harvest_rate = 0._r8
     do h_index = 1,hlm_num_lu_harvest_cats
        if (patch_anthro_disturbance_label .eq. primaryforest) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_VH1" .or. &
                hlm_harvest_catnames(h_index) .eq. "HARVEST_VH2") then
              harvest_rate = harvest_rate + hlm_harvest_rates(h_index)
           endif
        else if (patch_anthro_disturbance_label .eq. secondaryforest .and. &
             secondary_age >= secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH1") then
              harvest_rate = harvest_rate + hlm_harvest_rates(h_index)
           endif
        else if (patch_anthro_disturbance_label .eq. secondaryforest .and. &
             secondary_age < secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH2" .or. &
                hlm_harvest_catnames(h_index) .eq. "HARVEST_SH3") then
              harvest_rate = harvest_rate + hlm_harvest_rates(h_index)
           endif
        endif
     end do

     !  Normalize by site-level primary or secondary forest fraction
     !  since harvest_rate is specified as a fraction of the gridcell
     ! also need to put a cap so as not to harvest more primary or secondary area than there is in a gridcell
     if (patch_anthro_disturbance_label .eq. primaryforest) then
        if (frac_site_primary .gt. fates_tiny) then
           harvest_rate = min((harvest_rate / frac_site_primary),frac_site_primary)
        else
           harvest_rate = 0._r8
        endif
     else
        if ((1._r8-frac_site_primary) .gt. fates_tiny) then
           harvest_rate = min((harvest_rate / (1._r8-frac_site_primary)),&
                (1._r8-frac_site_primary))
        else
           harvest_rate = 0._r8
        endif
     endif

     ! calculate today's harvest rate
     ! whether to harvest today has already been determined by IsItLoggingTime
     ! for icode == 2, icode < 0, and icode > 10000 apply the annual rate one time (no calc)
     ! Bad logging event flag is caught in IsItLoggingTime, so don't check it here
     icode = int(logging_event_code)
     if(icode .eq. 1) then
        ! Logging is turned off - not sure why we need another switch
        harvest_rate = 0._r8
     else if(icode .eq. 3) then
        ! Logging event every day - this may not work due to the mortality exclusivity
        harvest_rate = harvest_rate / hlm_days_per_year
     else if(icode .eq. 4) then
        ! logging event once a month
        if(hlm_current_day.eq.1  ) then
           harvest_rate = harvest_rate / months_per_year
        end if
     end if

   end subroutine get_harvest_rate_area

   ! ============================================================================

   subroutine logging_litter_fluxes(currentSite, currentPatch, newPatch, patch_site_areadis)

      ! -------------------------------------------------------------------------------------------
      !
      !  DESCRIPTION:
      !  Carbon going from ongoing mortality into CWD pools. 
      !  This module includes only those fluxes associated with a disturbance generated by logging.
      !  Purpose: 
      !	  1) move logging-associated carbon to CWD and litter pool
      !   2) move the logging trunk from live into product pool 
      !   3) generate fluxes used in carbon balance checking
      !  E.g,:
      !  Remove trunk of logged trees from litter/CWD
      !  Add other parts of logged trees and all parts of collaterally and mechanically 
      !  damaged trees into CWD/litter  
      !
      !  This routine is only called if logging disturbance is the dominant disturbance.
      !
      !
      !  Note: The litter losses due to disturbance in the logging case is almost
      !        exactly like the natural tree-fall case.  The big differences are that
      !        the mortality rates governing the fluxes, follow a different rule set.
      !        We also compute an export flux (product) that does not go to litter.  
      !
      !  Trunk Product Flux: Only usable wood is exported from a site, substracted by a 
      !        transportation loss fraction. This is the above-ground portion of the bole,
      !        and only boles associated with direct-logging, not inftrastructure or
      !        collateral damage mortality.
      !        
      ! -------------------------------------------------------------------------------------------


      !USES:
      use SFParamsMod,  only : SF_val_cwd_frac
      use EDtypesMod,   only : area
      use EDtypesMod,   only : ed_site_type
      use EDtypesMod,   only : ed_patch_type
      use EDtypesMod,   only : ed_cohort_type
      use FatesAllometryMod , only : carea_allom


      ! !ARGUMENTS:
      type(ed_site_type)  , intent(inout), target  :: currentSite 
      type(ed_patch_type) , intent(inout), target  :: currentPatch
      type(ed_patch_type) , intent(inout), target  :: newPatch
      real(r8)            , intent(in)             :: patch_site_areadis

      !LOCAL VARIABLES:
      type(ed_cohort_type), pointer      :: currentCohort
      type(site_massbal_type), pointer   :: site_mass
      type(site_fluxdiags_type), pointer :: flux_diags
      type(litter_type),pointer          :: new_litt
      type(litter_type),pointer          :: cur_litt

      real(r8) :: direct_dead         ! Mortality count through direct logging
      real(r8) :: indirect_dead       ! Mortality count through: impacts, infrastructure and collateral damage
      real(r8) :: trunk_product_site  ! flux of carbon in trunk products exported off site      [ kgC/site ] 
                                      ! (note we are accumulating over the patch, but scale is site level)
      real(r8) :: delta_litter_stock  ! flux of carbon in total litter flux                     [ kgC/site ]
      real(r8) :: delta_biomass_stock ! total flux of carbon through mortality (litter+product) [ kgC/site ]
      real(r8) :: delta_individual    ! change in plant number through mortality [ plants/site ]
      real(r8) :: leaf_litter         ! Leafy biomass transferred through mortality [kgC/site]
      real(r8) :: root_litter         ! Rooty + storage biomass transferred through mort [kgC/site]
      real(r8) :: ag_wood             ! above ground wood mass [kg]
      real(r8) :: bg_wood             ! below ground wood mass [kg]
      real(r8) :: remainder_area      ! current patch's remaining area after donation [m2]
      real(r8) :: leaf_m              ! leaf element mass [kg]
      real(r8) :: fnrt_m              ! fineroot element mass [kg]
      real(r8) :: sapw_m              ! sapwood element mass [kg]
      real(r8) :: store_m             ! storage element mass [kg]
      real(r8) :: struct_m            ! structure element mass [kg]
      real(r8) :: repro_m             ! reproductive mass [kg]
      real(r8) :: retain_frac         ! fraction of litter retained in the donor patch
      real(r8) :: donate_frac         ! fraction of litter sent to newly formed patch
      real(r8) :: dcmpy_frac          ! fraction going into each decomposability pool
      integer  :: dcmpy               ! index for decomposability pools
      integer  :: element_id          ! parteh global element index
      integer  :: pft                 ! pft index
      integer  :: c                   ! cwd index
      integer  :: nlevsoil            ! number of soil layers
      integer  :: ilyr                ! soil layer loop index
      integer  :: el                  ! elemend loop index
      

      nlevsoil = currentSite%nlevsoil

      ! If/when sending litter fluxes to the old patch, we divide the total 
      ! mass sent to that patch, by the area it will have remaining
      ! after it donates area.
      ! i.e. subtract the area it is donating.
      
      remainder_area = currentPatch%area - patch_site_areadis


      ! Calculate the fraction of litter to be retained versus donated
      ! vis-a-vis the new and donor patch
      
      retain_frac = (1.0_r8-harvest_litter_localization) * &
            remainder_area/(newPatch%area+remainder_area)
      donate_frac = 1.0_r8-retain_frac

      do el = 1,num_elements
         
         element_id = element_list(el)
         site_mass => currentSite%mass_balance(el)
         flux_diags=> currentSite%flux_diags(el)
         cur_litt  => currentPatch%litter(el)   ! Litter pool of "current" patch
         new_litt  => newPatch%litter(el)       ! Litter pool of "new" patch
         

         ! Zero some site level accumulator diagnsotics
         trunk_product_site  = 0.0_r8
         delta_litter_stock  = 0.0_r8
         delta_biomass_stock = 0.0_r8
         delta_individual    = 0.0_r8

         ! -----------------------------------------------------------------------------
         ! Part 1: Send parts of dying plants to the litter pool.
         ! -----------------------------------------------------------------------------
         
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
               direct_dead   = currentCohort%n * currentCohort%lmort_direct
               indirect_dead = currentCohort%n * &
                     (currentCohort%lmort_collateral + currentCohort%lmort_infra)
               
            else

               ! This routine is only called during disturbance.  The litter
               ! fluxes from non-disturbance generating mortality are 
               ! handled in EDPhysiology.  Disturbance generating mortality
               ! are those cohorts in the top canopy layer, or those
               ! plants that were impacted. Thus, no direct dead can occur
               ! here, and indirect are impacts.

               if(int(prt_params%woody(pft)) == itrue) then
                  direct_dead   = 0.0_r8
                  indirect_dead = logging_coll_under_frac * &
                       (1._r8-currentPatch%fract_ldist_not_harvested) * currentCohort%n * &
                       (patch_site_areadis/currentPatch%area)   !kgC/site/day
               else
                  ! If the cohort of interest is grass, it will not experience
                  ! any mortality associated with the logging disturbance
                  direct_dead   = 0.0_r8
                  indirect_dead = 0.0_r8
               end if
            end if
            
            if( (element_id .eq. carbon12_element) .and. &
               hlm_use_planthydro == itrue ) then
               call AccumulateMortalityWaterStorage(currentSite, &
                     currentCohort,(direct_dead+indirect_dead))
            end if

            ! ----------------------------------------------------------------------------------------
            ! Handle woody litter flux for non-bole components of biomass
            ! This litter is distributed between the current and new patches, &
            ! not to any other patches. This is really the eventually area of the current patch &
            ! (currentPatch%area-patch_site_areadis) +patch_site_areadis...
            ! For the new patch, only some fraction of its land area (patch_areadis/np%area) is 
            ! derived from the current patch, so we need to multiply by patch_areadis/np%area
            ! ----------------------------------------------------------------------------------------

            call set_root_fraction(currentSite%rootfrac_scr, pft, currentSite%zi_soil)
         
            ag_wood = (direct_dead+indirect_dead) * (struct_m + sapw_m ) * &
                  prt_params%allom_agb_frac(currentCohort%pft)
            bg_wood = (direct_dead+indirect_dead) * (struct_m + sapw_m ) * &
                  (1._r8 - prt_params%allom_agb_frac(currentCohort%pft))
         
            do c = 1,ncwd-1
               
               new_litt%ag_cwd(c)     = new_litt%ag_cwd(c) + &
                     ag_wood * SF_val_CWD_frac(c) * donate_frac/newPatch%area
               cur_litt%ag_cwd(c)     = cur_litt%ag_cwd(c) + &
                     ag_wood * SF_val_CWD_frac(c) * retain_frac/remainder_area

               do ilyr = 1,nlevsoil
                  
                  new_litt%bg_cwd(c,ilyr) = new_litt%bg_cwd(c,ilyr) + &
                        bg_wood * currentSite%rootfrac_scr(ilyr) * &
                        SF_val_CWD_frac(c) * donate_frac/newPatch%area
                  
                  cur_litt%bg_cwd(c,ilyr) = cur_litt%bg_cwd(c,ilyr) + &
                        bg_wood * currentSite%rootfrac_scr(ilyr) * &
                        SF_val_CWD_frac(c) * retain_frac/remainder_area
               end do

               
               ! Diagnostics on fluxes into the AG and BG CWD pools
               flux_diags%cwd_ag_input(c) = flux_diags%cwd_ag_input(c) + & 
                    SF_val_CWD_frac(c) * ag_wood
               
               flux_diags%cwd_bg_input(c) = flux_diags%cwd_bg_input(c) + & 
                    SF_val_CWD_frac(c) * bg_wood
            
               ! Diagnostic specific to resource management code
               if( element_id .eq. carbon12_element) then
                   delta_litter_stock  = delta_litter_stock  + &
                         (ag_wood + bg_wood) * SF_val_CWD_frac(c)
               end if

            enddo
            
            ! ----------------------------------------------------------------------------------------
            ! Handle litter flux for the trunk wood of infrastucture and collateral damage mort
            ! ----------------------------------------------------------------------------------------
            
            ag_wood = indirect_dead * (struct_m + sapw_m ) * &
                  prt_params%allom_agb_frac(currentCohort%pft)
            bg_wood = indirect_dead * (struct_m + sapw_m ) * &
                  (1._r8 - prt_params%allom_agb_frac(currentCohort%pft))

            new_litt%ag_cwd(ncwd) = new_litt%ag_cwd(ncwd) + ag_wood * &
                  SF_val_CWD_frac(ncwd) * donate_frac/newPatch%area

            cur_litt%ag_cwd(ncwd) = cur_litt%ag_cwd(ncwd) + ag_wood * &
                  SF_val_CWD_frac(ncwd) * retain_frac/remainder_area
            
            do ilyr = 1,nlevsoil
               
               new_litt%bg_cwd(ncwd,ilyr) = new_litt%bg_cwd(ncwd,ilyr) + &
                     bg_wood * currentSite%rootfrac_scr(ilyr) * &
                     SF_val_CWD_frac(ncwd) * donate_frac/newPatch%area
               
               cur_litt%bg_cwd(ncwd,ilyr) = cur_litt%bg_cwd(ncwd,ilyr) + &
                     bg_wood * currentSite%rootfrac_scr(ilyr) * &
                     SF_val_CWD_frac(ncwd) * retain_frac/remainder_area

            end do

            flux_diags%cwd_ag_input(ncwd) = flux_diags%cwd_ag_input(ncwd) + & 
                 SF_val_CWD_frac(ncwd) * ag_wood
            
            flux_diags%cwd_bg_input(ncwd) = flux_diags%cwd_bg_input(ncwd) + & 
                 SF_val_CWD_frac(ncwd) * bg_wood

            if( element_id .eq. carbon12_element) then
                delta_litter_stock  = delta_litter_stock + &
                      (ag_wood+bg_wood) * SF_val_CWD_frac(ncwd)
            end if

            ! ---------------------------------------------------------------------------------------
            ! Handle below-ground trunk flux for directly logged trees (c = ncwd)
            ! ----------------------------------------------------------------------------------------
            
            bg_wood = direct_dead * (struct_m + sapw_m ) * SF_val_CWD_frac(ncwd) * &
                  (1._r8 - prt_params%allom_agb_frac(currentCohort%pft))

            do ilyr = 1,nlevsoil
                new_litt%bg_cwd(ncwd,ilyr) = new_litt%bg_cwd(ncwd,ilyr) + &
                      bg_wood * currentSite%rootfrac_scr(ilyr) * &
                      donate_frac/newPatch%area
                
                cur_litt%bg_cwd(ncwd,ilyr) = cur_litt%bg_cwd(ncwd,ilyr) + &
                      bg_wood * currentSite%rootfrac_scr(ilyr) * &
                      retain_frac/remainder_area
            end do
            
            flux_diags%cwd_bg_input(ncwd) = flux_diags%cwd_bg_input(ncwd) + &
                  bg_wood
            
            ! ----------------------------------------------------------------------------------------
            ! Handle harvest (export, flux-out) flux for the above ground boles 
            ! In this case a fraction (export_frac) of the boles from direct logging are 
            ! exported off-site, while the remainder (1-export_frac) is added to the litter pools.
            ! 
            ! Losses to the system as a whole, for C-balancing (kGC/site/day)
            ! Site level product, (kgC/site, accumulated over simulation)
            ! ----------------------------------------------------------------------------------------

            ag_wood = direct_dead * (struct_m + sapw_m ) * &
                  prt_params%allom_agb_frac(currentCohort%pft) * &
                  SF_val_CWD_frac(ncwd)

            trunk_product_site = trunk_product_site + &
                  ag_wood * logging_export_frac

            ! This is for checking the total mass balance [kg/site/day]
            site_mass%wood_product = site_mass%wood_product + &
                  ag_wood * logging_export_frac
            
            new_litt%ag_cwd(ncwd) = new_litt%ag_cwd(ncwd) + ag_wood * &
                  (1._r8-logging_export_frac)*donate_frac/newPatch%area
            
            cur_litt%ag_cwd(ncwd) = cur_litt%ag_cwd(ncwd) + ag_wood * &
                  (1._r8-logging_export_frac)*retain_frac/remainder_area

            ! ---------------------------------------------------------------------------
            ! Handle fluxes of leaf, root and storage carbon into litter pools. 
            !  (none of these are exported)
            ! ---------------------------------------------------------------------------
            
            leaf_litter = (direct_dead+indirect_dead)*(leaf_m + repro_m)
            root_litter = (direct_dead+indirect_dead)*(fnrt_m + store_m)


            do dcmpy=1,ndcmpy

               dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)

               new_litt%leaf_fines(dcmpy) = new_litt%leaf_fines(dcmpy) + &
                    leaf_litter * donate_frac/newPatch%area * dcmpy_frac
               
               cur_litt%leaf_fines(dcmpy) = cur_litt%leaf_fines(dcmpy) + &
                    leaf_litter * retain_frac/remainder_area * dcmpy_frac

               dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
               do ilyr = 1,nlevsoil
                  new_litt%root_fines(dcmpy,ilyr) = new_litt%root_fines(dcmpy,ilyr) + &
                       root_litter * currentSite%rootfrac_scr(ilyr) * dcmpy_frac * &
                       donate_frac/newPatch%area
                  
                  cur_litt%root_fines(dcmpy,ilyr) = cur_litt%root_fines(dcmpy,ilyr) + &
                       root_litter * currentSite%rootfrac_scr(ilyr) * dcmpy_frac * &
                       retain_frac/remainder_area
               end do
            end do
               
            ! track as diagnostic fluxes
            flux_diags%leaf_litter_input(pft) = flux_diags%leaf_litter_input(pft) + & 
                 leaf_litter
            
            flux_diags%root_litter_input(pft) = flux_diags%root_litter_input(pft) + & 
                 root_litter
            
            ! Logging specific diagnostics
            ! ----------------------------------------------------------------------------------------
            
            ! Note that litter stock also has terms above in the CWD loop
            if( element_id .eq. carbon12_element) then
                delta_litter_stock  = delta_litter_stock  + &
                      leaf_litter         + &
                      root_litter
                
                delta_biomass_stock = delta_biomass_stock + &
                      leaf_litter         + &
                      root_litter         + &
                      (direct_dead+indirect_dead) * (struct_m + sapw_m)
                
                delta_individual    = delta_individual    + &
                      direct_dead         + &
                      indirect_dead
            end if

            currentCohort => currentCohort%taller
         end do

         ! Update the amount of carbon exported from the site through logging
         ! operations.  Currently we assume only above-ground portion
         ! of the tree bole that experienced "direct" logging is exported
         ! This portion is known as "trunk_product_site
         
         if(element_id .eq. carbon12_element) then
            currentSite%resources_management%trunk_product_site  = &
                  currentSite%resources_management%trunk_product_site + &
                  trunk_product_site
            
            currentSite%resources_management%delta_litter_stock  = &
                  currentSite%resources_management%delta_litter_stock + &
                  delta_litter_stock
            
            currentSite%resources_management%delta_biomass_stock = &
                  currentSite%resources_management%delta_biomass_stock + &
                  delta_biomass_stock
            
            currentSite%resources_management%delta_individual    = &
                  currentSite%resources_management%delta_individual + &
                  delta_individual

         end if

      end do

      ! Not sure why this is called here, but I suppose it can't hurt
      ! (rgk 06-2019)

      currentCohort => newPatch%shortest
      do while(associated(currentCohort))
         call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread, &
               currentCohort%pft,currentCohort%c_area)
         currentCohort => currentCohort%taller
      enddo
      
      return
   end subroutine logging_litter_fluxes

end module EDLoggingMortalityMod
