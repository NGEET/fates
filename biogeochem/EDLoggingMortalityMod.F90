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
   use FatesConstantsMod , only : rsnbl_math_prec
   use FatesCohortMod    , only : fates_cohort_type
   use FatesPatchMod     , only : fates_patch_type
   use EDTypesMod        , only : site_massbal_type
   use EDTypesMod        , only : site_fluxdiags_type
   use EDTypesMod        , only : elem_diag_type
   use FatesLitterMod    , only : ncwd
   use FatesLitterMod    , only : ndcmpy
   use FatesLitterMod    , only : litter_type
   use FatesLitterMod    , only : adjust_SF_CWD_frac
   use EDTypesMod        , only : ed_site_type
   use EDTypesMod        , only : ed_resources_management_type
   use FatesConstantsMod , only : dtype_ilog
   use FatesConstantsMod , only : dtype_ifall
   use FatesConstantsMod , only : dtype_ifire
   use EDTypesMod        , only : area_inv
   use FatesConstantsMod , only : n_landuse_cats
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
   use FatesInterfaceTypesMod , only : bc_in_type
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
   use FatesInterfaceTypesMod , only : hlm_use_luh
   use FatesConstantsMod , only : itrue,ifalse
   use FatesGlobals      , only : endrun => fates_endrun 
   use FatesGlobals      , only : fates_log
   use FatesGlobals      , only : fates_global_verbose
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use FatesPlantHydraulicsMod, only : AccumulateMortalityWaterStorage
   use PRTGenericMod     , only : carbon12_element
   use PRTGenericMod     , only : sapw_organ, struct_organ, leaf_organ
   use PRTGenericMod     , only : fnrt_organ, store_organ, repro_organ
   use FatesAllometryMod , only : set_root_fraction
   use FatesConstantsMod , only : primaryland, secondaryland, secondary_age_threshold, nocomp_bareground_land
   use FatesConstantsMod , only : fates_tiny
   use FatesConstantsMod , only : months_per_year, days_per_sec, years_per_day, g_per_kg
   use FatesConstantsMod , only : hlm_harvest_area_fraction
   use FatesConstantsMod , only : hlm_harvest_carbon
   use FatesConstantsMod, only : fates_check_param_set
   use FatesConstantsMod, only : fates_no_harvest_debt, fates_with_harvest_debt, fates_bypass_harvest_debt

   use FatesInterfaceTypesMod , only : numpft
   use FatesLandUseChangeMod, only : GetInitLanduseHarvestRate
   use FatesLandUseChangeMod, only : GetLUHStatedata
     
   implicit none
   private

   logical, protected :: logging_time   ! If true, logging should be 
                                        ! performed during the current time-step

   logical, parameter :: debug = .false.
   
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
   public :: get_harvestable_carbon
   public :: get_harvest_rate_carbon
   public :: get_harvest_debt
   public :: UpdateHarvestC

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

   subroutine LoggingMortality_frac( currentSite, bc_in, pft_i, dbh, canopy_layer, lmort_direct, &
                                     lmort_collateral,lmort_infra, l_degrad, &
                                     hlm_harvest_rates, hlm_harvest_catnames, &
                                     hlm_harvest_units, &
                                     patch_land_use_label, secondary_age, &
                                     current_fates_landuse_state_vector, harvestable_forest_c, &
                                     harvest_tag)

     ! Arguments
      type(ed_site_type), intent(inout), target :: currentSite     ! site structure
      type(bc_in_type), intent(in) :: bc_in
      integer,  intent(in)  :: pft_i            ! pft index 
      real(r8), intent(in)  :: dbh              ! diameter at breast height (cm)
      integer,  intent(in)  :: canopy_layer     ! canopy layer of this cohort
      real(r8), intent(in) :: hlm_harvest_rates(:) ! annual harvest rate per hlm category
      character(len=64), intent(in) :: hlm_harvest_catnames(:) ! names of hlm harvest categories
      integer, intent(in) :: hlm_harvest_units     ! unit type of hlm harvest rates: [area vs. mass]
      integer, intent(in) :: patch_land_use_label    ! patch level land_use_label
      real(r8), intent(in) :: secondary_age     ! patch level age_since_anthro_disturbance
      real(r8), intent(in) :: harvestable_forest_c(:)  ! total harvestable forest carbon 
                                                       ! of all hlm harvest categories
      real(r8), intent(in) :: current_fates_landuse_state_vector(n_landuse_cats)  ! [m2/m2]
      real(r8), intent(out) :: lmort_direct     ! direct (harvestable) mortality fraction
      real(r8), intent(out) :: lmort_collateral ! collateral damage mortality fraction
      real(r8), intent(out) :: lmort_infra      ! infrastructure mortality fraction
      real(r8), intent(out) :: l_degrad         ! fraction of trees that are not killed
                                                ! but suffer from forest degradation (i.e. they
                                                ! are moved to newly-anthro-disturbed secondary
                                                ! forest patch)
      integer, intent(out) :: harvest_tag(:)    ! tag to record the harvest status 
                                                ! for the calculation of harvest debt in C-based
                                                ! harvest mode

      ! Local variables
      integer :: cur_harvest_tag ! the harvest tag of the cohort today
      real(r8) :: harvest_rate ! the final harvest rate to apply to this cohort today
      real(r8) :: state_vector(n_landuse_cats)
      logical  :: site_secondaryland_first_exceeding_min
      real(r8) :: secondary_young_fraction  ! what fraction of secondary land is young secondary land

      ! todo: probably lower the dbhmin default value to 30 cm
      ! todo: change the default logging_event_code to 1 september (-244)
      ! todo: change the default logging_direct_frac to 1.0 for cmip inputs
      ! todo: check outputs against the LUH2 carbon data
      ! todo: eventually set up distinct harvest practices, each with a set of input paramaeters
      ! todo: implement harvested carbon inputs

      ! The transition_landuse_from_off_to_on is for handling the special case of the first timestep after leaving potential
      ! vegetation mode. In this case, all prior historical land-use, including harvest, needs to be applied on that first day.
      ! So logging rates on that day are what is required to deforest exactly the amount of primary lands that will give the
      ! amount of secondary lands dictated by the land use state vector for that year, rather than whatever the continuous
      ! logging rate for that year is supposed to be according to the land use transition matrix.
      if (.not. currentSite%transition_landuse_from_off_to_on) then

         ! Check if the secondaryland exceeds the minimum if in landuse mode
         site_secondaryland_first_exceeding_min = .false.
         if (hlm_use_luh .eq. itrue) then
            call GetLUHStatedata(bc_in, state_vector)
            site_secondaryland_first_exceeding_min =  (state_vector(secondaryland) .gt. currentSite%min_allowed_landuse_fraction) &
                 .and. (.not. currentSite%landuse_vector_gt_min(secondaryland))
         end if

         ! if the total intended area of secondary lands are less than what we can consider without having too-small patches,
         ! or if that was the case until just now, then there is special logic
         if (site_secondaryland_first_exceeding_min) then
            if ( patch_land_use_label .eq. primaryland) then
               harvest_rate = state_vector(secondaryland) / state_vector(primaryland)
               write(fates_log(), *) 'applying state_vector(secondaryland) to plants.', pft_i
            else
               harvest_rate = 0._r8
            endif

            ! For area-based harvest, harvest_tag shall always be fates_bypass_harvest_debt (not applicable).
            harvest_tag = fates_bypass_harvest_debt
            cur_harvest_tag = fates_bypass_harvest_debt

         elseif (logging_time) then

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

               ! if using the classic CLM/ELM surface file, the variable names are:
               ! HARVEST_VH1 = harvest from primary forest
               ! HARVEST_VH2 = harvest from primary non-forest
               ! HARVEST_SH1 = harvest from secondary mature forest
               ! HARVEST_SH2 = harvest from secondary young forest
               ! HARVEST_SH3 = harvest from secondary non-forest (assume this is young for biomass)

               ! if using the direct LUH2 drivers, the variable names are instead (if using area-based logging):
               ! 'primf_harv', 'primn_harv', 'secmf_harv', 'secyf_harv', 'secnf_harv'

               secondary_young_fraction = currentSite%get_secondary_young_fraction()

               ! Get the area-based harvest rates based on info passed to FATES from the boundary condition
               call get_harvest_rate_area (patch_land_use_label, hlm_harvest_catnames, &
                    hlm_harvest_rates, current_fates_landuse_state_vector, secondary_young_fraction, secondary_age, harvest_rate)

               ! For area-based harvest, harvest_tag shall always be 2 (not applicable).
               harvest_tag = fates_bypass_harvest_debt
               cur_harvest_tag = fates_bypass_harvest_debt

               if (fates_global_verbose()) then
                  write(fates_log(), *) 'Successfully Read Harvest Rate from HLM.', hlm_harvest_rates(:), harvest_rate
               end if

            else if (hlm_use_lu_harvest == itrue .and. hlm_harvest_units == hlm_harvest_carbon) then
               ! 2=use carbon from hlm
               ! shall call another subroutine, which transfers biomass/carbon into fraction

               call get_harvest_rate_carbon (patch_land_use_label, hlm_harvest_catnames, &
                    hlm_harvest_rates, secondary_age, harvestable_forest_c, &
                    harvest_rate, harvest_tag, cur_harvest_tag)

               if (fates_global_verbose()) then
                  write(fates_log(), *) 'Successfully Read Harvest Rate from HLM.', hlm_harvest_rates(:), harvest_rate, harvestable_forest_c
               end if

            endif

         else
            harvest_rate = 0._r8
            ! For area-based harvest, harvest_tag shall always be 2 (not applicable).
            harvest_tag = fates_bypass_harvest_debt
            cur_harvest_tag = fates_bypass_harvest_debt
         endif

         ! transfer of area to secondary land is based on overall area affected, not just logged crown area
         ! l_degrad accounts for the affected area between logged crowns
         if(prt_params%woody(pft_i) == itrue)then ! only set logging rates for trees
            if (cur_harvest_tag == fates_no_harvest_debt .or. cur_harvest_tag == fates_bypass_harvest_debt) then
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
         ! the logic below is mainly to prevent conversion of bare ground land; everything else should be primary at this point.
         lmort_direct     = 0.0_r8
         lmort_collateral = 0.0_r8
         lmort_infra      = 0.0_r8
         l_degrad         = 0.0_r8
         if ( patch_land_use_label .eq. primaryland ) then
            call GetInitLanduseHarvestRate(bc_in, currentSite%min_allowed_landuse_fraction, &
                 harvest_rate, currentSite%landuse_vector_gt_min)
            if(prt_params%woody(pft_i) == itrue)then
               lmort_direct     = harvest_rate
            else if (canopy_layer .eq. 1) then
               l_degrad         = harvest_rate
            endif
         else if ( patch_land_use_label .ne. nocomp_bareground_land ) then
            write(fates_log(),*) 'trying to transition away from something that isnt either primary or bare ground,'
            write(fates_log(),*) 'on what should be a first timestep away from potential vegetation. This should not happen.'
            write(fates_log(),*) 'exiting'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         endif
      endif

   end subroutine LoggingMortality_frac


   ! ============================================================================

   subroutine get_harvest_rate_area (patch_land_use_label, hlm_harvest_catnames, hlm_harvest_rates, &
                 current_fates_landuse_state_vector, secondary_young_fraction, secondary_age, harvest_rate)


     ! -------------------------------------------------------------------------------------------
     !
     !  DESCRIPTION:
     !  get the area-based harvest rates based on info passed to FATES from the boundary conditions in.
     !  assumes logging_time == true

      ! Arguments
      real(r8), intent(in) :: hlm_harvest_rates(:) ! annual harvest rate per hlm category
      character(len=64), intent(in) :: hlm_harvest_catnames(:) ! names of hlm harvest categories
      integer, intent(in) :: patch_land_use_label    ! patch level land_use_label
      real(r8), intent(in) :: secondary_age     ! patch level age_since_anthro_disturbance
      real(r8), intent(in) :: current_fates_landuse_state_vector(n_landuse_cats)  ! [m2/m2]
      real(r8), intent(in) :: secondary_young_fraction  ! what fraction of secondary land is young secondary land
      real(r8), intent(out) :: harvest_rate

      ! Local Variables
      integer :: h_index   ! for looping over harvest categories
      integer :: icode   ! Integer equivalent of the event code (parameter file only allows reals)
      real(r8) :: frac_site_primary
      real(r8) :: frac_site_secondary
      real(r8) :: frac_not_bareground

     ! Loop around harvest categories to determine the annual hlm harvest rate for the current cohort based on patch history info
     harvest_rate = 0._r8
     do h_index = 1,hlm_num_lu_harvest_cats
        if (patch_land_use_label .eq. primaryland) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_VH1"  .or. &
                hlm_harvest_catnames(h_index) .eq. "HARVEST_VH2"  .or. &
                hlm_harvest_catnames(h_index) .eq. "primf_harv"  .or. &
                hlm_harvest_catnames(h_index) .eq. "primn_harv") then
              harvest_rate = harvest_rate + hlm_harvest_rates(h_index)
           endif
        else if (patch_land_use_label .eq. secondaryland .and. &
             secondary_age >= secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH1"  .or. &
                hlm_harvest_catnames(h_index) .eq. "secmf_harv") then
              harvest_rate = harvest_rate + hlm_harvest_rates(h_index)
           endif
        else if (patch_land_use_label .eq. secondaryland .and. &
             secondary_age < secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH2" .or. &
                hlm_harvest_catnames(h_index) .eq. "HARVEST_SH3"  .or. &
                hlm_harvest_catnames(h_index) .eq. "secyf_harv"  .or. &
                hlm_harvest_catnames(h_index) .eq. "secnf_harv") then
              harvest_rate = harvest_rate + hlm_harvest_rates(h_index)
           endif
        endif
     end do

     !  Normalize by site-level primary or secondary forest fraction
     !  since harvest_rate is specified as a fraction of the gridcell
     !  also need to put a cap so as not to harvest more primary or secondary area than there is in a gridcell
     !  For secondary, also need to normalize by the young/old fraction.
     !  Lastly, we need to remove the bare ground fraction since the harvest rates are per unit area of the not-bare-ground fraction.
     frac_site_primary = current_fates_landuse_state_vector(primaryland)
     frac_site_secondary = current_fates_landuse_state_vector(secondaryland)
     frac_not_bareground = sum(current_fates_landuse_state_vector(:))
     if (patch_land_use_label .eq. primaryland) then
        if (frac_site_primary .gt. fates_tiny .and. frac_not_bareground .gt. fates_tiny) then
           harvest_rate = min((harvest_rate / (frac_site_primary / frac_not_bareground)),1._r8)
        else
           harvest_rate = 0._r8
        endif
     else if (patch_land_use_label .eq. secondaryland) then
        ! the .gt. -0.5 in the next line is because frac_site_secondary returns -1 if no secondary area.
        if (frac_site_secondary .gt. fates_tiny .and. frac_site_secondary .gt. -0.5_r8 .and. frac_not_bareground .gt. fates_tiny) then
           if (secondary_age .lt. secondary_age_threshold) then
              harvest_rate = min((harvest_rate / ((frac_site_secondary / frac_not_bareground) * secondary_young_fraction)), 1._r8)
           else
              harvest_rate = min((harvest_rate / ((frac_site_secondary / frac_not_bareground) * (1._r8 - secondary_young_fraction))), 1._r8)
           endif
        else
           harvest_rate = 0._r8
        endif
     else
        harvest_rate = 0._r8
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
        if(hlm_current_day.eq.1) then
           harvest_rate = harvest_rate / months_per_year
        end if
     end if

   end subroutine get_harvest_rate_area


   ! ============================================================================

   subroutine get_harvestable_carbon (csite, site_area, hlm_harvest_catnames, harvestable_forest_c )

     !USES:
     use SFParamsMod,  only : SF_val_cwd_frac

     ! -------------------------------------------------------------------------------------------
     !
     !  DESCRIPTION:
     !  get the total carbon availale for harvest for three different harvest categories:
     !  primary forest, secondary mature forest and secondary young forest
     !  under two different scenarios:
     !  harvestable carbon: aggregate all cohorts matching the dbhmin harvest criteria
     !
     !  this subroutine shall be called outside the patch loop
     !  output will be used to estimate the area-based harvest rate (get_harvest_rate_carbon)
     !  for each cohort.

     ! Arguments
     type(ed_site_type), intent(in), target :: csite
     real(r8), intent(in) :: site_area    ! temporary variable
     character(len=64), intent(in) :: hlm_harvest_catnames(:) ! names of hlm harvest categories

     real(r8), intent(out) :: harvestable_forest_c(hlm_num_lu_harvest_cats)

     ! Local Variables
     type(fates_patch_type), pointer  :: currentPatch
     type(fates_cohort_type), pointer :: currentCohort
     real(r8) :: harvestable_patch_c     ! patch level total carbon available for harvest, kgC site-1
     real(r8) :: harvestable_cohort_c    ! cohort level total carbon available for harvest, kgC site-1
     real(r8) :: sapw_m      ! Biomass of sap wood
     real(r8) :: struct_m    ! Biomass of structural organs
     integer :: pft         ! Index of plant functional type
     integer :: h_index     ! for looping over harvest categories

     ! Initialization
     harvestable_forest_c = 0._r8

     ! loop over patches
     currentPatch => csite%oldest_patch
     do while (associated(currentPatch))
        harvestable_patch_c = 0._r8
        currentCohort => currentPatch%tallest

        do while (associated(currentCohort))
           pft = currentCohort%pft

           ! only account for cohorts matching the following conditions
           if(int(prt_params%woody(pft)) == 1)then ! only set logging rates for trees
              sapw_m   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
              struct_m = currentCohort%prt%GetState(struct_organ, carbon12_element)
              ! logging_direct_frac shall be 1 for LUH2 driven simulation and global simulation
              ! in site level study logging_direct_frac shall be surveyed
              ! unit:  [kgC ] = [kgC/plant] * [plant/ha] * [ha/ 10k m2] * [ m2 area ]
              harvestable_cohort_c = logging_direct_frac * ( sapw_m + struct_m ) * &
                     prt_params%allom_agb_frac(currentCohort%pft) * &
                     SF_val_CWD_frac(ncwd) * logging_export_frac * &
                     currentCohort%n * area_inv * site_area

              ! No harvest for trees without canopy 
              if (currentCohort%canopy_layer>=1) then
                 ! logging amount are based on dbh min and max criteria
                 if (currentCohort%dbh >= logging_dbhmin .and. .not. &
                       ((logging_dbhmax < fates_check_param_set) .and. (currentCohort%dbh >= logging_dbhmax )) ) then
                    ! Harvestable C: aggregate cohorts fit the criteria
                    harvestable_patch_c = harvestable_patch_c + harvestable_cohort_c
                 end if
              end if
           end if
           currentCohort => currentCohort%shorter
        end do

        ! judge which category the current patch belong to
        ! since we have not separated forest vs. non-forest
        ! all carbon belongs to the forest categories
        do h_index = 1,hlm_num_lu_harvest_cats
           if (currentPatch%land_use_label .eq. primaryland) then
              ! Primary
              if(hlm_harvest_catnames(h_index) .eq. "HARVEST_VH1") then
                 harvestable_forest_c(h_index) = harvestable_forest_c(h_index) + harvestable_patch_c
              end if
           else if (currentPatch%land_use_label .eq. secondaryland .and. &
                currentPatch%age_since_anthro_disturbance >= secondary_age_threshold) then
              ! Secondary mature
              if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH1") then
                 harvestable_forest_c(h_index) = harvestable_forest_c(h_index) + harvestable_patch_c
              end if
           else if (currentPatch%land_use_label .eq. secondaryland .and. &
                currentPatch%age_since_anthro_disturbance < secondary_age_threshold) then
              ! Secondary young
              if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH2") then
                 harvestable_forest_c(h_index) = harvestable_forest_c(h_index) + harvestable_patch_c
              end if
           end if
        end do
        currentPatch => currentPatch%younger
     end do

   end subroutine get_harvestable_carbon

   ! ============================================================================

   subroutine get_harvest_rate_carbon (patch_land_use_label, hlm_harvest_catnames, &
                 hlm_harvest_rates, secondary_age, harvestable_forest_c, &
                 harvest_rate, harvest_tag, cur_harvest_tag)

     ! -------------------------------------------------------------------------------------------
     !
     !  DESCRIPTION:
     !  get the carbon-based harvest rates based on info passed to FATES from the boundary conditions in.
     !  assumes logging_time == true

      ! Arguments
      real(r8), intent(in) :: hlm_harvest_rates(:) ! annual harvest rate per hlm category
      character(len=64), intent(in) :: hlm_harvest_catnames(:) ! names of hlm harvest categories
      integer, intent(in) :: patch_land_use_label    ! patch level land_use_label
      real(r8), intent(in) :: secondary_age     ! patch level age_since_anthro_disturbance
      real(r8), intent(in) :: harvestable_forest_c(:)  ! site level forest c matching criteria available for harvest, kgC site-1
      real(r8), intent(out) :: harvest_rate      ! area fraction
      integer,  intent(inout) :: harvest_tag(:)  ! This harvest tag can be raused to patch level but since all
                                                 ! logging functions happen within cohort loop we can only put the 
                                                 ! calculation here. Can think about optimizing the logging calculation
                                                 ! in the future.
      integer,  intent(out), optional :: cur_harvest_tag  ! harvest tag of the current cohort

      ! Local Variables
      integer :: h_index   ! for looping over harvest categories
      integer :: icode   ! Integer equivalent of the event code (parameter file only allows reals)
      real(r8) :: harvest_rate_c    ! Temporary variable, kgC site-1
      real(r8) :: harvest_rate_supply  ! Temporary variable, kgC site-1

     ! This subroutine follows the same logic of get_harvest_rate_area
     ! Loop over harvest categories to determine the hlm harvest rate demand and actual harvest rate for the 
     ! current cohort based on patch history info

     ! Initialize local variables
     harvest_rate = 0._r8
     harvest_rate_c = 0._r8
     harvest_rate_supply = 0._r8
     harvest_tag(:) = fates_bypass_harvest_debt

     ! Since we have five harvest categories from forcing data but in FATES non-forest harvest
     ! is merged with forest harvest, we only have three logging type in FATES (primary, secondary
     ! mature and secondary young).
     ! Get the harvest rate from HLM
     do h_index = 1,hlm_num_lu_harvest_cats
        if (patch_land_use_label .eq. primaryland) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_VH1"  .or. &
                hlm_harvest_catnames(h_index) .eq. "HARVEST_VH2") then
              harvest_rate_c = harvest_rate_c + hlm_harvest_rates(h_index)
           endif
        else if (patch_land_use_label .eq. secondaryland .and. &
             secondary_age >= secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH1") then
              harvest_rate_c = harvest_rate_c + hlm_harvest_rates(h_index)
           endif
        else if (patch_land_use_label .eq. secondaryland .and. &
             secondary_age < secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH2" .or. &
                hlm_harvest_catnames(h_index) .eq. "HARVEST_SH3") then
              harvest_rate_c = harvest_rate_c + hlm_harvest_rates(h_index)
           endif
        endif
     end do

     ! Determine harvest status (succesful or not)
     ! Here only three categories are used
     do h_index = 1,hlm_num_lu_harvest_cats
        if (patch_land_use_label .eq. primaryland) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_VH1" ) then
              if(harvestable_forest_c(h_index) >= harvest_rate_c) then
                 harvest_rate_supply = harvest_rate_supply + harvestable_forest_c(h_index)
                 harvest_tag(h_index) = fates_no_harvest_debt
              else
                 harvest_tag(h_index) = fates_with_harvest_debt
              end if
           end if
        else if (patch_land_use_label .eq. secondaryland .and. &
              secondary_age >= secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH1" ) then
              if(harvestable_forest_c(h_index) >= harvest_rate_c) then
                 harvest_rate_supply = harvest_rate_supply + harvestable_forest_c(h_index)
                 harvest_tag(h_index) = fates_no_harvest_debt
              else
                 harvest_tag(h_index) = fates_with_harvest_debt
              end if
           end if
        else if (patch_land_use_label .eq. secondaryland .and. &
              secondary_age < secondary_age_threshold) then
           if(hlm_harvest_catnames(h_index) .eq. "HARVEST_SH2" ) then
               if(harvestable_forest_c(h_index) >= harvest_rate_c) then
                  harvest_rate_supply = harvest_rate_supply + harvestable_forest_c(h_index)
                  harvest_tag(h_index) = fates_no_harvest_debt
               else
                  harvest_tag(h_index) = fates_with_harvest_debt
               end if
           end if
        end if
     end do

     ! If any harvest category available, assign to cur_harvest_tag and trigger logging event
     if(present(cur_harvest_tag))then
       cur_harvest_tag = minval(harvest_tag)
     end if

     ! Transfer carbon-based harvest rate to area-based harvest rate
     if (harvest_rate_supply > rsnbl_math_prec .and. harvest_rate_supply > harvest_rate_c) then
        harvest_rate = harvest_rate_c / harvest_rate_supply
     else
        ! If we force harvest rate to 1 when we don't have enough C, we will produce
        ! primary patch with no area, which cannot be terminated under nocomp mode.
        ! So we still keep the harvest rate to 0 for now.
        harvest_rate = 0._r8
     end if

     ! Prevent the generation of tiny secondary patches
     if(harvest_rate < 1e-8) harvest_rate = 0._r8

     ! For carbon-based harvest rate, normalizing by site-level primary or secondary forest fraction
     ! is not needed

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

   end subroutine get_harvest_rate_carbon

   ! ============================================================================

   subroutine logging_litter_fluxes(currentSite, currentPatch, newPatch, patch_site_areadis, bc_in)

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
      use SFParamsMod,       only : SF_val_cwd_frac
      use EDtypesMod,        only : area
      use EDtypesMod,        only : ed_site_type
      use FatesCohortMod,    only : fates_cohort_type
      use FatesConstantsMod, only : rsnbl_math_prec
      use FatesAllometryMod, only : carea_allom


      ! !ARGUMENTS:
      type(ed_site_type)  , intent(inout), target  :: currentSite 
      type(fates_patch_type) , intent(inout), target  :: currentPatch
      type(fates_patch_type) , intent(inout), target  :: newPatch
      real(r8)            , intent(in)             :: patch_site_areadis
      type(bc_in_type)    , intent(in)             :: bc_in


      !LOCAL VARIABLES:
      type(fates_cohort_type), pointer   :: currentCohort
      type(site_massbal_type), pointer   :: site_mass
      type(elem_diag_type), pointer      :: elflux_diags
      type(litter_type),pointer          :: new_litt
      type(litter_type),pointer          :: cur_litt

      real(r8) :: direct_dead            ! Mortality count through direct logging
      real(r8) :: indirect_dead          ! Mortality count through: impacts, infrastructure and collateral damage
      real(r8) :: trunk_product_site     ! flux of carbon in trunk products exported off site      [ kgC/site ] 
                                         ! (note we are accumulating over the patch, but scale is site level)
      real(r8) :: delta_litter_stock     ! flux of carbon in total litter flux                     [ kgC/site ]
      real(r8) :: delta_biomass_stock    ! total flux of carbon through mortality (litter+product) [ kgC/site ]
      real(r8) :: delta_individual       ! change in plant number through mortality [ plants/site ]
      real(r8) :: leaf_litter            ! Leafy biomass transferred through mortality [kgC/site]
      real(r8) :: root_litter            ! Rooty + storage biomass transferred through mort [kgC/site]
      real(r8) :: ag_wood                ! above ground wood mass [kg]
      real(r8) :: bg_wood                ! below ground wood mass [kg]
      real(r8) :: remainder_area         ! current patch's remaining area after donation [m2]
      real(r8) :: leaf_m                 ! leaf element mass [kg]
      real(r8) :: fnrt_m                 ! fineroot element mass [kg]
      real(r8) :: sapw_m                 ! sapwood element mass [kg]
      real(r8) :: store_m                ! storage element mass [kg]
      real(r8) :: struct_m               ! structure element mass [kg]
      real(r8) :: repro_m                ! reproductive mass [kg]
      real(r8) :: retain_frac            ! fraction of litter retained in the donor patch
      real(r8) :: retain_m2              ! area normalization for litter mass destined to old patch [m-2]
      real(r8) :: donate_m2              ! area normalization for litter mass destined to new patch [m-2]
      real(r8) :: dcmpy_frac             ! fraction going into each decomposability pool
      integer  :: dcmpy                  ! index for decomposability pools
      integer  :: element_id             ! parteh global element index
      integer  :: pft                    ! pft index
      integer  :: c                      ! cwd index
      integer  :: nlevsoil               ! number of soil layers
      integer  :: ilyr                   ! soil layer loop index
      integer  :: el                     ! elemend loop index
      real(r8) :: SF_val_CWD_frac_adj(4) !Updated wood partitioning to CWD based on dbh

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

      if(remainder_area > rsnbl_math_prec) then
         retain_m2 = retain_frac/remainder_area
         donate_m2 = (1.0_r8-retain_frac)/newPatch%area
      else
         retain_m2 = 0._r8
         donate_m2 = 1._r8/newPatch%area
      end if
  
  
      do el = 1,num_elements

         element_id = element_list(el)
         site_mass => currentSite%mass_balance(el)
         elflux_diags=> currentSite%flux_diags%elem(el)
         cur_litt  => currentPatch%litter(el)   ! Litter pool of "current" patch
         new_litt  => newPatch%litter(el)        ! Litter pool of "new" patch
         
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

               if(prt_params%woody(pft) == itrue) then
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

            call set_root_fraction(currentSite%rootfrac_scr, pft, &
                 currentSite%zi_soil, &
                 bc_in%max_rooting_depth_index_col)
         
            ag_wood = (direct_dead+indirect_dead) * (struct_m + sapw_m ) * &
                  prt_params%allom_agb_frac(currentCohort%pft)
            bg_wood = (direct_dead+indirect_dead) * (struct_m + sapw_m ) * &
                  (1._r8 - prt_params%allom_agb_frac(currentCohort%pft))
	
	    !adjust how wood is partitioned between the cwd classes based on cohort dbh
            call adjust_SF_CWD_frac(currentCohort%dbh,ncwd,SF_val_CWD_frac,SF_val_CWD_frac_adj) 
            
	    do c = 1,ncwd-1
               
               new_litt%ag_cwd(c)     = new_litt%ag_cwd(c) + &
                     ag_wood * SF_val_CWD_frac_adj(c) * donate_m2
               cur_litt%ag_cwd(c)     = cur_litt%ag_cwd(c) + &
                     ag_wood * SF_val_CWD_frac_adj(c) * retain_m2

               do ilyr = 1,nlevsoil
                  
                  new_litt%bg_cwd(c,ilyr) = new_litt%bg_cwd(c,ilyr) + &
                        bg_wood * currentSite%rootfrac_scr(ilyr) * &
                        SF_val_CWD_frac_adj(c) * donate_m2
                  
                  cur_litt%bg_cwd(c,ilyr) = cur_litt%bg_cwd(c,ilyr) + &
                        bg_wood * currentSite%rootfrac_scr(ilyr) * &
                        SF_val_CWD_frac_adj(c) * retain_m2
               end do

               
               ! Diagnostics on fluxes into the AG and BG CWD pools
               elflux_diags%cwd_ag_input(c) = elflux_diags%cwd_ag_input(c) + & 
                    SF_val_CWD_frac_adj(c) * ag_wood
               
               elflux_diags%cwd_bg_input(c) = elflux_diags%cwd_bg_input(c) + & 
                    SF_val_CWD_frac_adj(c) * bg_wood
            
               ! Diagnostic specific to resource management code
               if( element_id .eq. carbon12_element) then
                   delta_litter_stock  = delta_litter_stock  + &
                         (ag_wood + bg_wood) * SF_val_CWD_frac_adj(c)
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
                  SF_val_CWD_frac_adj(ncwd) * donate_m2

            cur_litt%ag_cwd(ncwd) = cur_litt%ag_cwd(ncwd) + ag_wood * &
                  SF_val_CWD_frac_adj(ncwd) * retain_m2
            
            do ilyr = 1,nlevsoil
               
               new_litt%bg_cwd(ncwd,ilyr) = new_litt%bg_cwd(ncwd,ilyr) + &
                     bg_wood * currentSite%rootfrac_scr(ilyr) * &
                     SF_val_CWD_frac_adj(ncwd) * donate_m2
               
               cur_litt%bg_cwd(ncwd,ilyr) = cur_litt%bg_cwd(ncwd,ilyr) + &
                     bg_wood * currentSite%rootfrac_scr(ilyr) * &
                     SF_val_CWD_frac_adj(ncwd) * retain_m2

            end do

            elflux_diags%cwd_ag_input(ncwd) = elflux_diags%cwd_ag_input(ncwd) + & 
                 SF_val_CWD_frac_adj(ncwd) * ag_wood
            
            elflux_diags%cwd_bg_input(ncwd) = elflux_diags%cwd_bg_input(ncwd) + & 
                 SF_val_CWD_frac_adj(ncwd) * bg_wood

            if( element_id .eq. carbon12_element) then
                delta_litter_stock  = delta_litter_stock + &
                      (ag_wood+bg_wood) * SF_val_CWD_frac_adj(ncwd)
            end if

            ! ---------------------------------------------------------------------------------------
            ! Handle below-ground trunk flux for directly logged trees (c = ncwd)
            ! ----------------------------------------------------------------------------------------
            
            bg_wood = direct_dead * (struct_m + sapw_m ) * SF_val_CWD_frac_adj(ncwd) * &
                  (1._r8 - prt_params%allom_agb_frac(currentCohort%pft))

            do ilyr = 1,nlevsoil
                new_litt%bg_cwd(ncwd,ilyr) = new_litt%bg_cwd(ncwd,ilyr) + &
                      bg_wood * currentSite%rootfrac_scr(ilyr) * &
                      donate_m2
                
                cur_litt%bg_cwd(ncwd,ilyr) = cur_litt%bg_cwd(ncwd,ilyr) + &
                      bg_wood * currentSite%rootfrac_scr(ilyr) * &
                      retain_m2
            end do
            
            elflux_diags%cwd_bg_input(ncwd) = elflux_diags%cwd_bg_input(ncwd) + &
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
                  SF_val_CWD_frac_adj(ncwd)

            trunk_product_site = trunk_product_site + &
                  ag_wood * logging_export_frac

            ! This is for checking the total mass balance [kg/site/day]
            site_mass%wood_product_harvest(pft) = site_mass%wood_product_harvest(pft) + &
                  ag_wood * logging_export_frac

            new_litt%ag_cwd(ncwd) = new_litt%ag_cwd(ncwd) + ag_wood * &
                  (1._r8-logging_export_frac)*donate_m2
            
            cur_litt%ag_cwd(ncwd) = cur_litt%ag_cwd(ncwd) + ag_wood * &
                  (1._r8-logging_export_frac)*retain_m2

            ! ---------------------------------------------------------------------------
            ! Handle fluxes of leaf, root and storage carbon into litter pools. 
            !  (none of these are exported)
            ! ---------------------------------------------------------------------------
            
            leaf_litter = (direct_dead+indirect_dead)*(leaf_m + repro_m)
            root_litter = (direct_dead+indirect_dead)*(fnrt_m + store_m)


            do dcmpy=1,ndcmpy

               dcmpy_frac = GetDecompyFrac(pft,leaf_organ,dcmpy)

               new_litt%leaf_fines(dcmpy) = new_litt%leaf_fines(dcmpy) + &
                    leaf_litter * donate_m2 * dcmpy_frac
               
               cur_litt%leaf_fines(dcmpy) = cur_litt%leaf_fines(dcmpy) + &
                    leaf_litter * retain_m2 * dcmpy_frac

               dcmpy_frac = GetDecompyFrac(pft,fnrt_organ,dcmpy)
               do ilyr = 1,nlevsoil
                  new_litt%root_fines(dcmpy,ilyr) = new_litt%root_fines(dcmpy,ilyr) + &
                       root_litter * currentSite%rootfrac_scr(ilyr) * dcmpy_frac * &
                       donate_m2
                  
                  cur_litt%root_fines(dcmpy,ilyr) = cur_litt%root_fines(dcmpy,ilyr) + &
                       root_litter * currentSite%rootfrac_scr(ilyr) * dcmpy_frac * &
                       retain_m2
               end do
            end do
               
            ! track as diagnostic fluxes
            elflux_diags%surf_fine_litter_input(pft) = elflux_diags%surf_fine_litter_input(pft) + & 
                 leaf_litter
            
            elflux_diags%root_litter_input(pft) = elflux_diags%root_litter_input(pft) + & 
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

         ! Amount of trunk mass exported off site [kg/m2]
         elflux_diags%exported_harvest = elflux_diags%exported_harvest + &
              trunk_product_site * area_inv

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
               currentCohort%pft,currentCohort%crowndamage,currentCohort%c_area)
         currentCohort => currentCohort%taller
      enddo
      
      return
   end subroutine logging_litter_fluxes


  ! =====================================================================================

   subroutine UpdateHarvestC(currentSite,bc_out)

      ! ----------------------------------------------------------------------------------
      ! Added by Shijie Shu.
      ! This subroutine is called when logging is completed and need to update 
      ! Harvested C flux in HLM.
      ! ----------------------------------------------------------------------------------
      use EDtypesMod             , only : ed_site_type
      use PRTGenericMod          , only : element_pos
      use PRTGenericMod          , only : carbon12_element
      use FatesInterfaceTypesMod , only : bc_out_type
  
      ! Arguments
      type(ed_site_type), intent(inout), target :: currentSite     ! site structure
      type(bc_out_type), intent(inout)          :: bc_out
  
      integer :: icode
      integer :: i_pft
      real(r8) :: unit_trans_factor
  

      ! Flush the older value before update
      bc_out%hrv_deadstemc_to_prod10c = 0._r8
      bc_out%hrv_deadstemc_to_prod100c = 0._r8
  
      ! Calculate the unit transfer factor (from kgC m-2 day-1 to gC m-2 s-1)
      unit_trans_factor = g_per_kg * days_per_sec

      ! harvest-associated wood product pools
      do i_pft = 1,numpft
         bc_out%hrv_deadstemc_to_prod10c = bc_out%hrv_deadstemc_to_prod10c + &
              currentSite%mass_balance(element_pos(carbon12_element))%wood_product_harvest(i_pft) * &
              area_inv * EDPftvarcon_inst%harvest_pprod10(i_pft) * unit_trans_factor
         bc_out%hrv_deadstemc_to_prod100c = bc_out%hrv_deadstemc_to_prod100c + &
              currentSite%mass_balance(element_pos(carbon12_element))%wood_product_harvest(i_pft) * &
              area_inv * (1._r8 - EDPftvarcon_inst%harvest_pprod10(i_pft)) * unit_trans_factor
      end do

      ! land-use-change-associated wood product pools
      do i_pft = 1,numpft
         bc_out%hrv_deadstemc_to_prod10c = bc_out%hrv_deadstemc_to_prod10c + &
              currentSite%mass_balance(element_pos(carbon12_element))%wood_product_landusechange(i_pft) * &
              area_inv * EDPftvarcon_inst%landusechange_pprod10(i_pft) * unit_trans_factor
         bc_out%hrv_deadstemc_to_prod100c = bc_out%hrv_deadstemc_to_prod100c + &
              currentSite%mass_balance(element_pos(carbon12_element))%wood_product_landusechange(i_pft) * &
              area_inv * (1._r8 - EDPftvarcon_inst%landusechange_pprod10(i_pft)) * unit_trans_factor
      end do

      return

   end subroutine UpdateHarvestC

   subroutine get_harvest_debt(site_in, bc_in, harvest_tag)

      !
      ! !DESCRIPTION:
      !
      ! Calculate if we have harvest debt for primary and secondary land
      ! Harvest debt is the accumulated total carbon 
      ! deficiency once the carbon amount available for harvest 
      ! is smaller than the harvest rate of forcing data.
      ! Harvest debt is calculated on site level
      ! TODO: we can define harvest debt as a fraction of the 
      ! harvest rate in the future
      ! Note: Non-forest harvest is accounted for under forest
      ! harvest, thus the harvest tag for non-forest is not applicable (= 2)
      !
      ! !ARGUMENTS:
      type(ed_site_type) , intent(inout), target :: site_in
      type(bc_in_type),    intent(in)         :: bc_in
      integer  :: harvest_tag(hlm_num_lu_harvest_cats)

      ! !LOCAL VARIABLES:
      integer  :: h_index
      real(r8) :: harvest_debt_pri
      real(r8) :: harvest_debt_sec_mature
      real(r8) :: harvest_debt_sec_young

      if(logging_time) then

         ! Initialize the local variables
         harvest_debt_pri = 0._r8
         harvest_debt_sec_mature = 0._r8
         harvest_debt_sec_young = 0._r8
 
         ! First we need to get harvest rate for all three categories
         do h_index = 1, hlm_num_lu_harvest_cats
            ! Primary forest harvest rate
            if(bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_VH1" .or. &
                bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_VH2" ) then
                  harvest_debt_pri = harvest_debt_pri + bc_in%hlm_harvest_rates(h_index)
            else if(bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_SH1") then
                harvest_debt_sec_mature = harvest_debt_sec_mature + bc_in%hlm_harvest_rates(h_index)
            else if(bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_SH2" .or. &
                     bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_SH3") then
                harvest_debt_sec_young = harvest_debt_sec_young + bc_in%hlm_harvest_rates(h_index)
            end if
         end do
         ! Next we get the harvest debt through the harvest tag 
         do h_index = 1, hlm_num_lu_harvest_cats
            if (harvest_tag(h_index) .eq. fates_with_harvest_debt) then
               if(bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_VH1") then
                  site_in%resources_management%harvest_debt = site_in%resources_management%harvest_debt + &
                      harvest_debt_pri
               else if(bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_SH1") then
                  site_in%resources_management%harvest_debt = site_in%resources_management%harvest_debt + &
                      harvest_debt_sec_mature
                  site_in%resources_management%harvest_debt_sec = site_in%resources_management%harvest_debt_sec + &
                      harvest_debt_sec_mature
               else if(bc_in%hlm_harvest_catnames(h_index) .eq. "HARVEST_SH2") then
                  site_in%resources_management%harvest_debt = site_in%resources_management%harvest_debt + &
                      harvest_debt_sec_young
                  site_in%resources_management%harvest_debt_sec = site_in%resources_management%harvest_debt_sec + &
                      harvest_debt_sec_young
               end if
            end if
         end do
      end if

   end subroutine get_harvest_debt

end module EDLoggingMortalityMod
