module FatesCumulativeMemoryMod
   !---~---
   !   This module contains procedures that update multiple cumulative/memory variables
   ! that drive leaf phenology, and potentially mortality, disturbances, and management.
   ! These sub-routines are related to FatesRunningMeanMod, however, they need a separate
   ! module because many of these procedures require access to EDTypesMod, which in turn
   ! depends on FatesRunningMeanMod.
   !---~---



   use EDBtranMod            , only : check_layer_water
   use EDParamsMod           , only : ED_val_phen_chiltemp
   use EDTypesMod            , only : area_inv
   use EDTypesMod            , only : ed_site_type
   use EDTypesMod            , only : num_vegtemp_mem
   use EDTypesMod            , only : numWaterMem
   use EDTypesMod            , only : phen_cstat_iscold
   use FatesAllometryMod     , only : set_root_fraction
   use FatesConstantsMod     , only : ndays_per_year
   use FatesConstantsMod     , only : nearzero
   use FatesConstantsMod     , only : r8   => fates_r8
   use FatesConstantsMod     , only : tfrz => t_water_freeze_k_1atm
   use FatesInterfaceTypesMod, only : bc_in_type
   use FatesInterfaceTypesMod, only : hlm_day_of_year
   use FatesInterfaceTypesMod, only : numpft
   use FatesPatchMod         , only : fates_patch_type

   implicit none
   private

   !---~---
   !   The main sub-routine that updates all cumulative and memory-related variables. 
   ! This is the only sub-routine that needs to be public, the specific procedures that
   ! update time, temperature and moisture memories are called by this sub-routine.
   !---~---
   public :: UpdateCumulativeMemoryVars

   !---~---
   !   Imposed soil matric potential lower bound for frozen or excessively dry soils,
   ! used when computing water stress.
   !---~---
   real(r8), public, parameter :: smp_lwr_bound = -1000000._r8 

contains

   subroutine UpdateCumulativeMemoryVars(currentSite,bc_in)
      !---~---
      !   This sub-routine updates multiple variables that retain cumulative or "memory"
      ! related information. These are mostly used for phenology, but some of the variables
      ! can be useful for disturbances, mortality and management.
      !---~---

      ! Arguments
      type(ed_site_type), intent(inout), target :: currentSite
      type(bc_in_type),   intent(in)            :: bc_in


      ! Update elapsed time
      call UpdatePhenologyDate(currentSite)

      ! Update temperature-related variables. Currently this is just the temperature
      !    memory, but it will eventually contain growing degree days and number of
      !    chilling days.
      call UpdateCumulativeThermal(currentSite,bc_in)

      ! Update moisture-related memory variables.
      call UpdateMemoryMoisture(currentSite,bc_in)

      return
   end subroutine UpdateCumulativeMemoryVars






   subroutine UpdatePhenologyDate(currentSite)
      !---~---
      !   This sub-routine simply updates the number of elapsed days in this simulation.
      ! The first day of the simulation is 1, and it continues monotonically, 
      ! indefinitely.
      !---~---

      ! Arguments
      type(ed_site_type), intent(inout), target :: currentSite

      !---~---
      !    Advance elapsed time. The only reason this is a site variable instead of a 
      ! global variable is that we need to save this information to the restart file,
      ! and we do not have global scalars in the restart file.
      !---~---
      currentSite%phen_model_date = currentSite%phen_model_date + 1

      return
   end subroutine UpdatePhenologyDate






   subroutine UpdateCumulativeThermal(currentSite,bc_in)
      !---~---
      !   Currently this sub-routine updates the canopy temperature memory variable only.
      ! In the future, this will also host updates to growing degree days and number of 
      ! chilling days.
      !---~---

      ! Arguments
      type(ed_site_type), intent(inout), target :: currentSite
      type(bc_in_type),   intent(in)            :: bc_in

      ! Local variables
      type(fates_patch_type), pointer :: cpatch    ! Current patch
      real(r8)   :: temp_in_C         ! Daily averaged temperature in degC
      real(r8)   :: temp_wgt          ! Canopy area weighting factor for daily average
                                      !    vegetation temperature calculation

      ! Find the average temperature for the previous 24 hours
      temp_in_C = 0._r8
      temp_wgt = 0._r8
      cpatch => CurrentSite%oldest_patch
      do while( associated(cpatch) )
         temp_in_C = temp_in_C + cpatch%tveg24%GetMean()*cpatch%total_canopy_area
         temp_wgt = temp_wgt + cpatch%total_canopy_area
         cpatch => cpatch%younger
      end do
      if ( temp_wgt > nearzero ) then
         temp_in_C = temp_in_C/temp_wgt - tfrz
      else
         ! If there is no canopy area, we use the vegetation temperature
         ! of the first patch, which is the forcing air temperature
         ! as defined in CLM/ELM. The forcing air temperature
         ! should be the same across all patches. (Although
         ! it is unlikely there are more than 1 in this scenario)
         temp_in_C = CurrentSite%oldest_patch%tveg24%GetMean() - tfrz
      end if


      ! Record temperature for the last 10 days to the memory variable.
      currentSite%vegtemp_memory(2:num_vegtemp_mem) = currentSite%vegtemp_memory(1:num_vegtemp_mem-1)
      currentSite%vegtemp_memory(1) = temp_in_C


      return
   end subroutine UpdateCumulativeThermal






   subroutine UpdateMemoryMoisture(currentSite,bc_in)
      !---~---
      !   This sub-routine updates the moisture-related memory variables. These variables
      ! are used for leaf phenology, and may become useful for other modules too.
      !---~---

      ! Arguments
      type(ed_site_type), intent(inout), target :: currentSite
      type(bc_in_type),   intent(in)            :: bc_in

      ! Local variables
      integer  :: ipft              ! Plant Functional Type index
      integer  :: i_wmem            ! Loop counter for water memory days
      integer  :: j                 ! Soil layer index
      integer  :: nlevroot          ! Number of rooting levels to consider
      real(r8) :: rootfrac_notop    ! Total rooting fraction excluding the top soil layer



      ! The soil memory variables are defined for each PFT
      pft_memory_loop: do ipft=1,numpft

         !    Update soil moisture information memory (we always track the last 10 days).
         ! We shift the memory from days to the previous day, and make room for current day
         do i_wmem = numWaterMem,2,-1
            currentSite%liqvol_memory(i_wmem,ipft) = currentSite%liqvol_memory(i_wmem-1,ipft)
            currentSite%smp_memory   (i_wmem,ipft) = currentSite%smp_memory   (i_wmem-1,ipft)
         end do

         ! Find the rooting depth distribution for PFT
         call set_root_fraction( currentSite%rootfrac_scr, ipft, currentSite%zi_soil, &
                                 bc_in%max_rooting_depth_index_col )
         nlevroot = max(2,min(ubound(currentSite%zi_soil,1),bc_in%max_rooting_depth_index_col))


         !    The top most layer is typically very thin (~ 2cm) and dries rather quickly. Despite
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

         !    Set the memory to be the weighted average of the liquid volumetric soil water,
         ! using the root fraction of each layer (except the topmost one) as the weighting
         ! factor.
         currentSite%liqvol_memory(1,ipft) = &
            sum( bc_in%h2o_liqvol_sl(2:nlevroot) * currentSite%rootfrac_scr(2:nlevroot) ) / &
            rootfrac_notop

         !    For the soil matric potential, we must check whether or not to include the layer
         ! based on the layer temperature.
         currentSite%smp_memory   (1,ipft) = 0._r8
         root_loop: do j = 2,nlevroot
            if ( check_layer_water(bc_in%h2o_liqvol_sl(j),bc_in%tempk_sl(j)) ) then
               currentSite%smp_memory   (1,ipft) = &
                  currentSite%smp_memory   (1,ipft) + &
                  bc_in%smp_sl(j) * currentSite%rootfrac_scr(j) / rootfrac_notop
            else
               ! Nominal extreme suction for frozen or unreasonably dry soil
               currentSite%smp_memory   (1,ipft) = &
                  currentSite%smp_memory   (1,ipft) + &
                  smp_lwr_bound * currentSite%rootfrac_scr(j)  / rootfrac_notop
            end if
         end do root_loop
      end do pft_memory_loop


      return
   end subroutine UpdateMemoryMoisture
end module FatesCumulativeMemoryMod
