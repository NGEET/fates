module FatesForestManagementMod

   ! ============================================================================
   ! Created by Shijie Shu on Sep 8. 2025
   ! This module contains all subroutines and functions related to planned 
   ! wood harvest and improved forest management.
   !-----------------------------------------------------------------------------------
   ! Improved forest management practies can be quantified as modifiers to change the
   ! logging disturbance rates at cohort-patch-site levels, thus quantify strategies as
   ! different priorties regarding size/age for harvest activity. 
   ! Site: site_in%resources_management%harvest_wp_scale. Adjust the overall harvest rate
   ! to match the demand.
   ! Patch: currentPatch%harvest_rate_scale. Modifier considering age priority.
   ! Cohort [LOCAL]: harvest_rate_scale_cohort. Modifier considering size priority.
   ! ============================================================================

   ! Uses
   use EDParamsMod            , only : logging_age_preference, logging_preference_options
   use EDParamsMod            , only : logging_ifm_harvest_scale
   use EDTypesMod             , only : AREA, AREA_INV
   use EDTypesMod             , only : ed_site_type
   use FatesPatchMod          , only : fates_patch_type
   use FatesCohortMod         , only : fates_cohort_type
   use FatesGlobals           , only : fates_log
   use FatesConstantsMod      , only : r8 => fates_r8
   use FatesConstantsMod      , only : logging_no_age_preference, logging_oldest_first
   use FatesConstantsMod      , only : logging_uniform_size, logging_double_rotation, logging_quadruple_rotation
   use FatesConstantsMod      , only : logging_logistic_size, logging_inv_logistic_size, logging_gaussian_size
   use FatesConstantsMod      , only : fates_unset_r8, fates_tiny
   use FatesConstantsMod      , only : primaryland, secondaryland
   use FatesConstantsMod      , only : min_harvest_rate, min_nocomp_pftfrac_perlanduse
   use FatesConstantsMod      , only : hlm_harvest_carbon
   use FatesLitterMod         , only : ncwd
   use FatesInterfaceTypesMod , only : bc_in_type
   use FatesInterfaceTypesMod , only : hlm_num_lu_harvest_cats, numpft
   use PRTParametersMod       , only : prt_params
   use PRTGenericMod          , only : carbon12_element
   use PRTGenericMod          , only : leaf_organ
   use PRTGenericMod          , only : fnrt_organ
   use PRTGenericMod          , only : sapw_organ
   use PRTGenericMod          , only : store_organ
   use PRTGenericMod          , only : repro_organ
   use PRTGenericMod          , only : struct_organ

   ! Subroutines
   public :: get_site_harvest_rate_scale
   public :: get_patch_harvest_rate_scale
   public :: get_harvestable_carbon
   private :: get_harvestable_patch_carbon
   private :: get_harvested_cohort_carbon

 contains

   ! ============================================================================

   subroutine get_site_harvest_rate_scale(site_in, bc_in, harvestable_forest_c, frac_site_primary, frac_site_secondary_mature)

      ! !USES:
      use EDLoggingMortalityMod , only : logging_time, LoggingMortality_frac

      ! !ARGUMENTS:
      type(ed_site_type) , intent(inout) :: site_in
      type(bc_in_type) , intent(in) :: bc_in
      real(r8), intent(in) :: harvestable_forest_c(hlm_num_lu_harvest_cats)
      real(r8), intent(in) :: frac_site_primary
      real(r8), intent(in) :: frac_site_secondary_mature
      !
      ! !LOCAL VARIABLES:
      type (fates_patch_type) , pointer :: currentPatch
      type (fates_cohort_type), pointer :: currentCohort

      real(r8) :: lmort_direct
      real(r8) :: lmort_collateral
      real(r8) :: lmort_infra
      real(r8) :: l_degrad
      real(r8) :: harvested_cohort_c   ! Variable to store cohort level harvested C
      integer :: harvest_tag(hlm_num_lu_harvest_cats)

      real(r8) :: target_wp
      real(r8) :: actual_wp
      integer :: it

      ! !PARAMETERS
      integer, parameter :: max_iteration = 5 ! Maximum iterations to obtain target harvest wood product. Currently I don't
                                              ! plan to move it into parameter files

      ! ==================================================================================================
      !  DESCRIPTION:
      !
      ! For size-dependet harvest, cetain size class might have fewer or no harvest if the available 
      ! inventory is less than the demand.
      ! We design the following algorithm to force the other nearby size class to fullfill these demands:
      ! Loop through cohorts to preview the summed harvested product C for both no harvest size priority,
      ! i.e., target wood product C
      ! (target_wp), and with size priority, i.e., actual harvested C  (actual_wp)
      ! Then calculate the ratio of target_wp to actual_wp as adjustment ratio
      ! The adjustment ratio will be applied as additional scaling factor to modify cohort-level harvest mortality
      ! One time adjustment might still not able to match the demand, we iterate this process no more than 5 times.
      ! ==================================================================================================

      !Initialize
      site_in%resources_management%harvest_wp_scale = 1._r8

      !Only run this part when considering size dependent priority
      if (logging_time .and. logging_preference_options >= logging_logistic_size) then

         iteration_loop: do it = 1, max_iteration

            target_wp = 0._r8
            actual_wp = 0._r8

            currentPatch => site_in%oldest_patch
            do while (associated(currentPatch))

               currentCohort => currentPatch%shortest
               do while(associated(currentCohort))

                  ! Target harvest product
                  call LoggingMortality_frac(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_layer, &
                        lmort_direct,lmort_collateral,lmort_infra,l_degrad,&
                        bc_in%hlm_harvest_rates, &
                        bc_in%hlm_harvest_catnames, &
                        bc_in%hlm_harvest_units, &
                        currentPatch%land_use_label, &
                        currentPatch%age_since_anthro_disturbance, &
                        frac_site_primary, &
                        frac_site_secondary_mature, &
                        harvestable_forest_c, &
                        currentPatch%harvest_rate_scale, &
                        harvest_tag, logging_ifm_harvest_scale, logging_uniform_size)

                  call get_harvested_cohort_carbon(currentCohort, lmort_direct, harvested_cohort_c)

                  target_wp = target_wp + harvested_cohort_c

                  ! Actual harvest product
                  call LoggingMortality_frac(currentCohort%pft, currentCohort%dbh, currentCohort%canopy_layer, &
                        lmort_direct,lmort_collateral,lmort_infra,l_degrad,&
                        bc_in%hlm_harvest_rates, &
                        bc_in%hlm_harvest_catnames, &
                        bc_in%hlm_harvest_units, &
                        currentPatch%land_use_label, &
                        currentPatch%age_since_anthro_disturbance, &
                        frac_site_primary, &
                        frac_site_secondary_mature, &
                        harvestable_forest_c, &
                        currentPatch%harvest_rate_scale, &
                        harvest_tag, site_in%resources_management%harvest_wp_scale, logging_preference_options)

                  call get_harvested_cohort_carbon(currentCohort, lmort_direct, harvested_cohort_c)

                  actual_wp = actual_wp + harvested_cohort_c

                  currentCohort => currentCohort%taller

               end do

               currentPatch => currentPatch%younger

            end do

            ! adjustment ratio
            if( target_wp / (actual_wp + 1e-7_r8) < 1.01 .and. target_wp / (actual_wp + 1e-7_r8) > 0.99 ) exit iteration_loop

            site_in%resources_management%harvest_wp_scale = site_in%resources_management%harvest_wp_scale * target_wp/(actual_wp+1e-7_r8)

            write(fates_log(),*) 'See adjustment factor for harvest scale:', site_in%resources_management%harvest_wp_scale

         end do iteration_loop

      end if  ! logging_time

   end subroutine get_site_harvest_rate_scale

   ! ============================================================================

   subroutine get_patch_harvest_rate_scale(site_in, bc_in, harvestable_forest_c, frac_site_primary, frac_site_secondary_mature)

      ! !USES:

      ! !ARGUMENTS:
      type(ed_site_type) , intent(inout) :: site_in
      type(bc_in_type) , intent(in) :: bc_in
      real(r8), intent(in) :: harvestable_forest_c(hlm_num_lu_harvest_cats)
      real(r8), intent(in) :: frac_site_primary
      real(r8), intent(in) :: frac_site_secondary_mature
      
      ! !LOCAL VARIABLES:
      type (fates_patch_type) , pointer :: currentPatch
      type (fates_cohort_type), pointer :: currentCohort

      integer :: npatches        ! Number of patches
      integer :: ipatch, curpft ! index, used in min-heap merge
      integer, allocatable :: jpatch(:)   ! index array, used in min-heap merge
      integer, allocatable :: pftlen(:)   ! array length, used in min-heap merge
      real(r8), allocatable :: minheap(:)            ! minimum heap array, used in min-heap merge
      integer, allocatable :: order_of_patches(:)    ! Record a copy of patch order 
      integer, allocatable :: order_of_patches_per_pft(:,:)  ! Record a copy of patch order for each pft
      real(r8), allocatable :: age_patches(:)   ! Record a copy of patch age
      real(r8), allocatable :: age_patches_per_pft(:,:)   ! Record a copy of patch age for each pft

      real(r8) :: frac_site      ! Temporary variable
      real(r8) :: harvest_rate     ! Temporary variable
      real(r8) :: harvestable_patch_c
      real(r8) :: remain_harvest_rate(hlm_num_lu_harvest_cats)    ! Temporary variable
      real(r8) :: temp_totc    ! for debug test, can be removed in released version
      integer :: h_index     ! index for looping over harvest categories
      integer :: harvest_tag(hlm_num_lu_harvest_cats)

      ! !PARAMETERS

      ! Initialize the rest of local variables
      npatches = 0
      frac_site = 1._r8
      remain_harvest_rate = fates_unset_r8  ! set to negative value to indicate uninitialized

      ! =====================================================================
      !  DESCRIPTION:
      !
      !  This modifier is for age_prioritized forest management
      !  The following code will calculate patch level harvest_rate_scale 
      !  for each secondary patch
      !  There're currently 2 options:
      !  1) uniformly harvest each patch (harvest_rate_scale = 1._r8)
      !  2) harvest the oldest patch first

      ! Initialize harvest_rate_scale
      currentPatch => site_in%oldest_patch

      do while (associated(currentPatch))
         ! Initialize harvest rate scale to 1.0 first
         currentPatch%harvest_rate_scale = 1._r8

         currentPatch => currentPatch%younger
      end do

      ! =====================================================================
      ! "Sort" patch based on patch age since anthropegenic disturbance
      ! =====================================================================
      ! Patch sequence in the data list (patchno) follows pseudo-age defined 
      ! in function GetPseudoPatchAge to simplify LUC related calculation.
      ! Here we don't modify the actual patch sequence but use another age 
      ! tag currentPatch%order_age_since_anthro to record patch age for 
      ! age-prioritized wood harvest strategy.
      !
      ! In FATES, PFT-patch relation can be different under diffrent modes. 
      ! This can be separated into 2 cases:
      ! Case 1: For most of FATES modes, we directly use pseudo-age since 
      ! there's no paired 1 patch - 1 PFT relation
      ! Case 2: For the more general case in global simulation that harvest 
      ! all PFTs equally under no competition mode, we need to sort again 
      ! using, actual patch age. For each PFT the patch sequence is already 
      ! sorted, here we need to loop over PFTs and merge these sorted 
      ! sequences into one sorted array through min-heap merge

      if (logging_time .and. logging_age_preference == logging_oldest_first) then

         ! First loop to count total number of patches and initialize order
         ! using patchno
         ! For most of the cases the process is over after this step
         currentPatch => site_in%oldest_patch
         npatches = 0

         do while (associated(currentPatch))
            npatches = npatches + 1
            currentPatch%order_age_since_anthro = currentPatch%patchno

            currentPatch => currentPatch%younger
         end do

         ! For case 2, i.e., nocomp mode only
         if (hlm_use_nocomp.eq.itrue) then

            allocate(order_of_patches(1:npatches))
            allocate(order_of_patches_per_pft(1:numpft,1:npatches))
            allocate(age_patches(1:npatches))
            allocate(age_patches_per_pft(1:numpft,1:npatches))
            allocate(minheap(1:numpft))
            allocate(jpatch(1:numpft))
            allocate(pftlen(1:numpft))

            ! Second loop to assign patch sequence, age and number of pacthes in each PFT to array
            ! thus can avoid for loop over data list
            currentPatch => site_in%oldest_patch
            ipatch = 0
            pftlen(:) = 0
            do while (associated(currentPatch))
               ipatch = ipatch + 1
               order_of_patches(ipatch) = currentPatch%order_age_since_anthro
               age_patches(ipatch) = currentPatch%age

               do curpft=1,numpft
                  if(currentPatch%nocomp_pft_label .eq. curpft) then
                     pftlen(curpft) = pftlen(curpft) + 1
                     order_of_patches_per_pft(curpft, pftlen(curpft)) = currentPatch%order_age_since_anthro
                     age_patches_per_pft(curpft, pftlen(curpft)) = currentPatch%age
                  end if
               end do

               currentPatch => currentPatch%younger
            end do

            ! Third loop to initailize min-heap, first element of each array
            minheap(:) = 0._r8
            do curpft=1,numpft
               if(pftlen(curpft) > 0) then
                  minheap(curpft) = age_patches_per_pft(curpft, 1)
               end if
            end do

            ! Fourth loop to get min from min-heap then insert 'age_patches_per_pft' and record order in 'order_of_patches'
            ! no need to loop over data list here
            jpatch(:) = 1
            do ipatch = 1, npacthes
               ! Find the minimum and assign the order
               curpft = minloc(minheap, dim=1)
               order_of_patches(ipatch) = order_of_patches_per_pft(curpft, jpatch(curpft))
               ! Move to the next element and insert into minheap array
               jpatch(curpft) = jpatch(curpft) + 1
               if(jpatch(curpft) <= pftlen(curpft)) then
                  minheap(curpft) = age_patches_per_pft(ipft, jpatch(curpft))
               else
                  minheap(curpft) = 0._r8
               end if
            end do

            ! Fifthth loop to assign order back to patches
            do ipatch = 1, npatches
               currentPatch => site_in%oldest_patch

               ! Look for the patch by usig patchno
               inner_loop: do while (associated(currentPatch))

                  if(order_of_patches(ipatch) .eq. currentPatch%patchno) then
                     currentPatch%order_age_since_anthro = ipatch
                     exit inner_loop
                  end if

                  currentPatch => currentPatch%younger
               end do inner_loop
            end do

            ! Clean temperary variables
            deallocate(order_of_patches)
            deallocate(order_of_patches_per_pft)
            deallocate(age_patches)
            deallocate(age_patches_per_pft)
            deallocate(minheap)
            deallocate(jpatch)
            deallocate(pftlen)

         end if  ! nocomp

         ! Loop through all patches and calculate the harvest rate scale based
         ! on oldest first strategy
         target_loop: do ipatch = 1, npatches

            currentPatch => site_in%oldest_patch
            search_loop: do while (associated(currentPatch))

               found_patch: if(currentPatch%order_age_since_anthro .eq. ipatch) then

                  ! Here age priority is based on age since anthropogenic disturbance
                  ! thus we skip the primary land since they all have this age = 0
                  if_secondary: if(currentPatch%land_use_label .eq. secondaryland) then
                      if(currentPatch%age_since_anthro_disturbance >= secondary_age_threshold) then
                         h_index = 3
                         frac_site = frac_site_secondary_mature
                      else
                         h_index = 4
                         frac_site = 1._r8 - frac_site_primary - frac_site_secondary_mature
                      end if
                      ! Obtain uniform harvest rate (area fraction) of the corresponding patch
                      if(bc_in%hlm_harvest_units == hlm_harvest_carbon) then
                          call get_harvest_rate_carbon (currentPatch%land_use_label, bc_in%hlm_harvest_catnames, &
                               bc_in%hlm_harvest_rates, currentPatch%age_since_anthro_disturbance, harvestable_forest_c, &
                               harvest_rate, harvest_tag)
                      else
                          call get_harvest_rate_area (currentPatch%land_use_label, bc_in%hlm_harvest_catnames, &
                               bc_in%hlm_harvest_rates, frac_site_primary, frac_site_secondary_mature, &
                               currentPatch%age_since_anthro_disturbance, harvest_rate)
                      end if
                      ! For the first time, initialize remain_harvest_rate
                      if(remain_harvest_rate(h_index) < 0) then
                         if(bc_in%hlm_harvest_units == hlm_harvest_carbon) then
                            ! In kgC ha-1, no need to scale
                            remain_harvest_rate(h_index) = bc_in%hlm_harvest_rates(h_index)
                         else
                            ! In area fraction (0 - 1) after scaling to the area of per patch land use type
                            remain_harvest_rate(h_index) = harvest_rate
                         end if
                      end if
                      ! Only calculate harvest rate scale for non-zero harvest rate
                      if (harvest_rate > min_harvest_rate) then
                         if(bc_in%hlm_harvest_units == hlm_harvest_carbon) then
                            ! Compare the patch level harvestable carbon and see if larger than remain_harvest_rate
                            ! Note the scaling factor applied on area fraction as harvest unit
                            call get_harvestable_patch_carbon(currentPatch, harvestable_patch_c)
                            if(harvestable_patch_c >= remain_harvest_rate(h_index)) then
                                currentPatch%harvest_rate_scale = remain_harvest_rate(h_index) / &
                                    (harvestable_patch_c + fates_tiny) / harvest_rate
                                remain_harvest_rate(h_index) = 0._r8
                            else
                                ! harvest almost 100%, leave 0.01% to prevent removing the patch completely, which may cause
                                ! model crashing under nocomp mode
                                currentPatch%harvest_rate_scale = (1._r8 - min_nocomp_pftfrac_perlanduse) / harvest_rate
                                remain_harvest_rate(h_index) = remain_harvest_rate(h_index) - (1._r8 - min_nocomp_pftfrac_perlanduse) * harvestable_patch_c
                            end if
                         else
                            ! Compare the patch area and see if larger than remain_harvest_rate
                            if((currentPatch%area/(AREA*frac_site)) >= remain_harvest_rate(h_index)) then
                                currentPatch%harvest_rate_scale = remain_harvest_rate(h_index) / &
                                    (currentPatch%area/(AREA*frac_site) + fates_tiny) / harvest_rate
                                remain_harvest_rate(h_index) = 0._r8
                            else
                                ! leave 1% to prevent removing the patch completely
                                currentPatch%harvest_rate_scale = (1._r8 - min_nocomp_pftfrac_perlanduse) / harvest_rate
                                remain_harvest_rate(h_index) = remain_harvest_rate(h_index) - (1._r8 - min_nocomp_pftfrac_perlanduse) * (currentPatch%area/(AREA*frac_site))
                            end if
                         end if
                      end if

                      ! For diagnosing the harvestable carbon
                      ! These codes can be removed in released version
                      temp_totc = 0._r8
                      currentCohort => currentPatch%tallest

                      do while (associated(currentCohort))

                         if (currentCohort%canopy_layer>=1) then
                            temp_totc = temp_totc + (currentCohort%prt%GetState(sapw_organ, carbon12_element) + &
                                currentCohort%prt%GetState(struct_organ, carbon12_element)) * currentCohort%n * AREA_INV
                         end if
                         currentCohort => currentCohort%shorter

                      end do

                      write(fates_log(),*) 'See patch number:', currentPatch%patchno
                      write(fates_log(),*) 'See patch age:', currentPatch%age
                      write(fates_log(),*) 'See patch total carbon (kgc m-2):', temp_totc
                      write(fates_log(),*) 'See patch age sec:', currentPatch%age_since_anthro_disturbance
                      write(fates_log(),*) 'See patch order sec:', currentPatch%order_age_since_anthro
                      write(fates_log(),*) 'See harvest index:', h_index
                      write(fates_log(),*) 'See harvest rate scale:', currentPatch%harvest_rate_scale
                      write(fates_log(),*) 'See original harvest rate:', bc_in%hlm_harvest_rates(h_index)
                      write(fates_log(),*) 'See harvest rate:', harvest_rate
                      write(fates_log(),*) 'See harvestable_forest_c:', harvestable_forest_c
                      write(fates_log(),*) 'See site fraction:', frac_site
                      write(fates_log(),*) 'See sec mature fraction:', frac_site_secondary_mature
                      write(fates_log(),*) 'See sec young fraction:', 1 - frac_site_primary - frac_site_secondary_mature
                      write(fates_log(),*) 'See area:', currentPatch%area/AREA
                      write(fates_log(),*) '==================== Next patch ================'

                  end if if_secondary

                  ! Exit current inner loop to find the next old ipatch
                  exit search_loop

               end if found_patch

               currentPatch => currentPatch%younger
            end do search_loop

         end do target_loop

      end if  ! logging_time

   end subroutine get_patch_harvest_rate_scale

   ! ============================================================================

   subroutine get_harvestable_carbon (csite, site_area, hlm_harvest_catnames, harvestable_forest_c )

     !USES:
     use SFParamsMod,  only : SF_val_cwd_frac

     ! -------------------------------------------------------------------------------------------
     !
     !  DESCRIPTION:
     !  get the total site level carbon availale for harvest for three different harvest categories:
     !  primary forest, secondary mature forest and secondary young forest
     !  under two different scenarios:
     !  harvestable carbon: aggregate all cohorts matching the dbhmin harvest criteria
     !
     !  this subroutine shall be called outside the patch loop
     !  output will be used to estimate the area-based harvest rate (get_harvest_rate_carbon)
     !  for each cohort.

     ! Arguments
     type(ed_site_type), intent(in), target :: csite
     real(r8), intent(in) :: site_area    ! The total site area
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
           if(int(prt_params%woody(pft)) == 1) then ! only set logging rates for trees
              sapw_m   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
              struct_m = currentCohort%prt%GetState(struct_organ, carbon12_element)
              ! logging_direct_frac shall be 1 for LUH2 driven simulation and global simulation
              ! in site level study logging_direct_frac shall be surveyed
              ! unit:  [kgC ] = [kgC/plant] * [plant/ha] * [ha/ 10k m2] * [ m2 area ]
              harvestable_cohort_c = logging_direct_frac * ( sapw_m + struct_m ) * &
                     prt_params%allom_agb_frac(currentCohort%pft) * &
                     SF_val_CWD_frac(ncwd) * logging_export_frac * &
                     currentCohort%n! * AREA_INV * site_area

              ! No harvest for trees without canopy (exclude 0)
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

   subroutine get_harvestable_patch_carbon ( currentPatch, harvestable_patch_c )

     !USES:
     use SFParamsMod,  only : SF_val_cwd_frac

     ! -------------------------------------------------------------------------------------------
     !
     !  DESCRIPTION:
     !  get the total carbon (ha-1) availale for harvest in current Patch
     !  harvestable carbon: aggregate all cohorts matching the dbhmin and dbhmax harvest criteria
     !
     !  this subroutine shall be called inside the patch loop
     !  output will be used for diagnosis or for scaling patch level harvest rate

     ! Arguments
     type(fates_patch_type) , intent(in), target  :: currentPatch

     real(r8), intent(out) :: harvestable_patch_c

     ! Local Variables
     type(fates_cohort_type), pointer :: currentCohort
     real(r8) :: harvestable_cohort_c    ! cohort level total carbon available for harvest, kgC site-1
     real(r8) :: sapw_m      ! Biomass of sap wood
     real(r8) :: struct_m    ! Biomass of structural organs
     integer :: pft         ! Index of plant functional type

     ! loop over cohorts
     harvestable_patch_c = 0._r8
     currentCohort => currentPatch%tallest

     do while (associated(currentCohort))
        pft = currentCohort%pft

        ! only account for cohorts matching the following conditions
        if(int(prt_params%woody(pft)) == 1) then ! only set logging rates for trees
           sapw_m   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
           struct_m = currentCohort%prt%GetState(struct_organ, carbon12_element)
           ! logging_direct_frac shall be 1 for LUH2 driven simulation and global simulation
           ! in site level study logging_direct_frac shall be surveyed
           ! unit:  [kgC ] = [kgC/plant] * [plant/ha] * [ha/ 10k m2] * [ m2 area ]
           harvestable_cohort_c = logging_direct_frac * ( sapw_m + struct_m ) * &
                  prt_params%allom_agb_frac(currentCohort%pft) * &
                  SF_val_CWD_frac(ncwd) * logging_export_frac * &
                  currentCohort%n! * AREA_INV * site_area

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

   end subroutine get_harvestable_patch_carbon


   ! ============================================================================

   subroutine get_harvested_cohort_carbon ( currentCohort, lmort_direct, harvested_cohort_c )

     !USES:
     use SFParamsMod,  only : SF_val_cwd_frac


     ! -------------------------------------------------------------------------------------------
     !
     !  DESCRIPTION:
     !  get the harvested carbon in current Cohort through lmort_direct
     !
     !  this subroutine shall be called inside the cohort loop
     !  output will be used for adjusting IFM wood harvest with size priority

     ! Arguments
     type(fates_cohort_type) , intent(in), target  :: currentCohort
     real(r8), intent(in)  :: lmort_direct

     real(r8), intent(out) :: harvested_cohort_c

     ! Local Variables
     real(r8) :: sapw_m      ! Biomass of sap wood
     real(r8) :: struct_m    ! Biomass of structural organs

     ! only account for cohorts matching the following conditions
     sapw_m   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
     struct_m = currentCohort%prt%GetState(struct_organ, carbon12_element)

     harvested_cohort_c = lmort_direct * ( sapw_m + struct_m ) * &
            prt_params%allom_agb_frac(currentCohort%pft) * &
            SF_val_CWD_frac(ncwd) * logging_export_frac * &
            currentCohort%n

   end subroutine get_harvested_cohort_carbon

   ! ============================================================================

end module FatesForestManagementMod
