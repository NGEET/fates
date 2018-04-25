module EDCanopyStructureMod

  ! =====================================================================================
  ! Code to determine whether the canopy is closed, and which plants are either in the 
  ! understorey or overstorey. This is obviosuly far too complicated for it's own good 
  ! =====================================================================================

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesGlobals          , only : fates_log
  use EDPftvarcon           , only : EDPftvarcon_inst
  use FatesAllometryMod     , only : carea_allom
  use EDCohortDynamicsMod   , only : copy_cohort, terminate_cohorts, fuse_cohorts
  use FatesAllometryMod     , only : tree_lai
  use FatesAllometryMod     , only : tree_sai
  use EDtypesMod            , only : ed_site_type, ed_patch_type, ed_cohort_type, ncwd
  use EDTypesMod            , only : nclmax
  use EDTypesMod            , only : nlevleaf
  use EDtypesMod            , only : AREA
  use FatesGlobals          , only : endrun => fates_endrun
  use FatesInterfaceMod     , only : hlm_days_per_year
  use FatesInterfaceMod     , only : numpft


  ! CIME Globals
  use shr_log_mod           , only : errMsg => shr_log_errMsg

  implicit none
  private

  public :: canopy_structure
  public :: canopy_spread
  public :: calc_areaindex
  public :: canopy_summarization
  public :: update_hlm_dynamics

  logical, parameter :: DEBUG=.false.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  real(r8), parameter :: area_target_precision = 1.0E-6_r8  ! Area conservation must be within this tolerance
  real(r8), parameter :: area_check_precision  = 1.0E-4_r8  ! Area conservation checks must be within this tolerance
  integer, parameter  :: max_layer_iterations  = 100        ! Don't let these loop hang indefinitely

  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

  ! ============================================================================
   subroutine canopy_structure( currentSite , bc_in )
      !
      ! !DESCRIPTION:
      ! create cohort instance
      !
      ! This routine allocates the 'canopy_layer' attribute to each cohort
      ! All top leaves in the same canopy layer get the same light resources.
      ! The first canopy layer is the 'canopy' or 'overstorey'. The second is the 'understorey'.
      ! More than two layers is not permitted at the moment
      ! Seeds germinating into the 3rd or higher layers are automatically removed. 
      !
      ! ------Perfect Plasticity-----
      ! The idea of these canopy layers derives originally from Purves et al. 2009
      ! Their concept is that, given enoughplasticity in canopy position, size, shape and depth
      ! all of the gound area will be filled perfectly by leaves, and additional leaves will have
      ! to exist in the understorey. 
      ! Purves et al. use the concept of 'Z*' to assume that the height required to attain a place in the
      ! canopy is spatially uniform. In this implementation, described in Fisher et al. (2010, New Phyt) we
      ! extent that concept to assume that position in the canopy has some random element, and that BOTH height
      ! and chance combine to determine whether trees get into the canopy. 
      ! Thus, when the canopy is closed and there is excess area, some of it must be demoted
      ! If we demote -all- the trees less than a given height, there is a massive advantage in being the cohort that is 
      ! the biggest when the canopy is closed. 
      ! In this implementation, the amount demoted, ('weight') is a function of the height weighted by the competitive exclusion
      ! parameter (ED_val_comp_excln). 

      ! Complexity in this routine results from a few things. 
      ! Firstly, the complication of the demotion amount sometimes being larger than the cohort area (for a very small, short cohort)
      ! Second, occasionaly, disturbance (specifically fire) can cause the canopy layer to become less than closed, 
      ! without changing the area of the patch. If this happens, then some of the plants in the lower layer need to be 'promoted' so 
      ! all of the routine has to happen in both the downwards and upwards directions. 
      !
      ! The order of events here is therefore:
      ! (The entire subroutine has a single outer 'patch' loop. 
      ! Section 1: figure out the total area, and whether there are >1 canopy layers at all. 
      !
      ! Sorts out cohorts into canopy and understorey layers...                              
      !
      ! !USES:

      use EDParamsMod, only : ED_val_comp_excln
      use EDtypesMod , only : ncwd
      use EDTypesMod , only : min_patch_area
      use EDTypesMod , only : val_check_ed_vars
      use FatesInterfaceMod, only : bc_in_type
      !
      ! !ARGUMENTS    
      type(ed_site_type) , intent(inout), target   :: currentSite
      type(bc_in_type), intent(in)                 :: bc_in

      !
      ! !LOCAL VARIABLES:
      type(ed_patch_type) , pointer :: currentPatch
      type(ed_cohort_type), pointer :: currentCohort
      integer  :: i_lyr                  ! current layer index
      integer  :: z                      ! Current number of canopy layers. (1= canopy, 2 = understorey) 
      integer  :: ipft
      real(r8) :: arealayer(nclmax+2)    ! Amount of plant area currently in each canopy layer
      integer  :: patch_area_counter     ! count iterations used to solve canopy areas
      logical  :: area_not_balanced      ! logical controlling if the patch layer areas
                                         ! have successfully been redistributed
      integer  :: return_code            ! math checks on variables will return>0 if problems exist
      integer, parameter  :: max_patch_iterations = 100
      

      !----------------------------------------------------------------------

      currentPatch => currentSite%oldest_patch    
      ! 
      ! zero site-level demotion / promotion tracking info
      currentSite%demotion_rate(:) = 0._r8
      currentSite%promotion_rate(:) = 0._r8
      currentSite%demotion_carbonflux = 0._r8
      currentSite%promotion_carbonflux = 0._r8

      !
      ! Section 1: Check  total canopy area.    
      !
      do while (associated(currentPatch)) ! Patch loop    

         ! ------------------------------------------------------------------------------
         ! Perform numerical checks on some cohort and patch structures
         ! ------------------------------------------------------------------------------

         call val_check_ed_vars(currentPatch,'co_n:co_dbh:pa_area',return_code)
         ! No need to make error message, already generated in math_check_ed_vars
         if(return_code>0) call endrun(msg=errMsg(sourcefile, __LINE__))

         ! canopy layer has a special bounds check
         currentCohort => currentPatch%tallest
         do while (associated(currentCohort))
            if( currentCohort%canopy_layer < 1 .or. currentCohort%canopy_layer > nclmax+1 ) then 
               write(fates_log(),*) 'lat:',currentSite%lat
               write(fates_log(),*) 'lon:',currentSite%lon
               write(fates_log(),*) 'BOGUS CANOPY LAYER: ',currentCohort%canopy_layer
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
            currentCohort => currentCohort%shorter
         enddo


         ! Does any layer have excess area in it? Keep going until it does not...
         patch_area_counter = 0
         area_not_balanced = .true.
         
         do while(area_not_balanced)
            
            ! ---------------------------------------------------------------------------
            ! Demotion Phase: Identify upper layers that are too full, and demote them to
            ! the layers below.
            ! ---------------------------------------------------------------------------
            
            ! Calculate how many layers we have in this canopy
            ! This also checks the understory to see if its crown 
            ! area is large enough to warrant a temporary sub-understory layer
            z = NumPotentialCanopyLayers(currentPatch,currentSite%spread,include_substory=.true.)
            
            do i_lyr = 1,z ! Loop around the currently occupied canopy layers. 
               call DemoteFromLayer(currentSite, currentPatch, i_lyr)
            end do
            
            ! Remove cohorts that are incredibly sparse
            call terminate_cohorts(currentSite, currentPatch, 1)
            
            call fuse_cohorts(currentPatch, bc_in)
            
            ! Remove cohorts for various other reasons
            call terminate_cohorts(currentSite, currentPatch, 2)

            
            ! ---------------------------------------------------------------------------------------
            ! Promotion Phase: Identify if any upper-layers are underful and layers below them
            ! have cohorts that can be split and promoted to the layer above.
            ! ---------------------------------------------------------------------------------------
            
            ! Re-calculate Number of layers without the false substory
            z = NumPotentialCanopyLayers(currentPatch,currentSite%spread,include_substory=.false.)

            ! We only promote if we have at least two layers
            if (z>1) then
               
               do i_lyr=1,z-1 
                  call PromoteIntoLayer(currentSite, currentPatch, i_lyr)
               end do
               
               ! Remove cohorts that are incredibly sparse
               call terminate_cohorts(currentSite, currentPatch, 1)
               
               call fuse_cohorts(currentPatch, bc_in)
               
               ! Remove cohorts for various other reasons
               call terminate_cohorts(currentSite, currentPatch, 2)
               
            end if
            
            ! ---------------------------------------------------------------------------------------
            ! Check on Layer Area (if the layer differences are not small
            ! Continue trying to demote/promote. Its possible on the first pass through,
            ! that cohort fusion has nudged the areas a little bit.
            ! ---------------------------------------------------------------------------------------
            
            z = NumPotentialCanopyLayers(currentPatch,currentSite%spread,include_substory=.false.)
            area_not_balanced = .false.
            do i_lyr = 1,z
               call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer(i_lyr))
               if( (arealayer(i_lyr)-currentPatch%area)  >  area_check_precision )then
                  area_not_balanced = .true.
               endif
            enddo
            
            ! ---------------------------------------------------------------------------------------
            ! Gracefully exit if too many iterations have gone by
            ! ---------------------------------------------------------------------------------------
            
            patch_area_counter = patch_area_counter + 1
            if(patch_area_counter > max_patch_iterations) then
               write(fates_log(),*) 'PATCH AREA CHECK NOT CLOSING'
               write(fates_log(),*) 'patch area:',currentpatch%area
               write(fates_log(),*) 'lat:',currentSite%lat
               write(fates_log(),*) 'lon:',currentSite%lon
	       write(fates_log(),*) 'spread:',currentSite%spread
               currentCohort => currentPatch%tallest
               do while (associated(currentCohort))  
                  write(fates_log(),*) 'coh ilayer:',currentCohort%canopy_layer
                  write(fates_log(),*) 'coh dbh:',currentCohort%dbh
                  write(fates_log(),*) 'coh pft:',currentCohort%pft
                  write(fates_log(),*) 'coh n:',currentCohort%n
                  write(fates_log(),*) 'coh carea:',currentCohort%c_area
		  ipft=currentCohort%pft
		  write(fates_log(),*) 'maxh:',EDPftvarcon_inst%allom_dbh_maxheight(ipft)
                  write(fates_log(),*) 'lmode: ',EDPftvarcon_inst%allom_lmode(ipft)
		  write(fates_log(),*) 'd2bl2: ',EDPftvarcon_inst%allom_d2bl2(ipft)
		  write(fates_log(),*) 'd2bl_ediff: ',EDPftvarcon_inst%allom_blca_expnt_diff(ipft)
		  write(fates_log(),*) 'd2ca_min: ',EDPftvarcon_inst%allom_d2ca_coefficient_min(ipft)
		  write(fates_log(),*) 'd2ca_max: ',EDPftvarcon_inst%allom_d2ca_coefficient_max(ipft)
                  currentCohort => currentCohort%shorter
               enddo
               call endrun(msg=errMsg(sourcefile, __LINE__))
            end if
            
         enddo ! do while(area_not_balanced)
            
            
         ! Set current canopy layer occupancy indicator. 
         currentPatch%NCL_p = min(nclmax,z)    

         ! -------------------------------------------------------------------------------------------
         ! if we are using "strict PPA", then calculate a z_star value as 
         ! the height of the smallest tree in the canopy 
         ! loop from top to bottom and locate the shortest cohort in level 1 whose shorter 
         ! neighbor is in level 2 set zstar as the ehight of that shortest level 1 cohort
         ! -------------------------------------------------------------------------------------------
         
         if ( ED_val_comp_excln .lt. 0.0_r8) then
            currentPatch%zstar = 0._r8
            currentCohort => currentPatch%tallest
            do while (associated(currentCohort))
               if(currentCohort%canopy_layer .eq. 2)then
                  if (associated(currentCohort%taller)) then
                     if (currentCohort%taller%canopy_layer .eq. 1 ) then
                        currentPatch%zstar = currentCohort%taller%hite
                     endif
                  endif
               endif
               currentCohort => currentCohort%shorter
            enddo
         endif
         
         currentPatch => currentPatch%younger
      enddo !patch

      return
   end subroutine canopy_structure

   
   ! ==============================================================================================
   

   subroutine DemoteFromLayer(currentSite,currentPatch,i_lyr)

      use EDParamsMod, only : ED_val_comp_excln
      use EDtypesMod , only : ncwd
      use SFParamsMod, only : SF_val_CWD_frac

      ! !ARGUMENTS
      type(ed_site_type), intent(inout), target  :: currentSite
      type(ed_patch_type), intent(inout), target :: currentPatch
      integer, intent(in)                        :: i_lyr   ! Current canopy layer of interest

      ! !LOCAL VARIABLES:
      type(ed_cohort_type), pointer :: currentCohort,copyc
      integer  :: i_cwd                  ! Index for CWD pool
      real(r8) :: cc_loss
      real(r8) :: lossarea
      real(r8) :: newarea
      real(r8) :: weight                 ! The amount of the total lost area from this cohort
      real(r8) :: sumdiff
      real(r8) :: sum_weights
      real(r8) :: arealayer              ! the area of the current canopy layer
      real(r8) :: rankordered_area_sofar ! the amount of total canopy area occupied by 
                                         ! cohorts upto this point
      integer  :: layer_area_counter

      ! First, determine how much total canopy area we have in this layer

      call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer)


      layer_area_counter = 0
      do while( (arealayer-currentPatch%area) > area_target_precision ) 

         ! Is this layer currently over-occupied? 
         ! In that case, we need to work out which cohorts to demote. 

         sumdiff  = 0.0_r8    
         rankordered_area_sofar = 0.0_r8
         currentCohort => currentPatch%tallest 
         do while (associated(currentCohort))

            call carea_allom(currentCohort%dbh,currentCohort%n, &
                  currentSite%spread,currentCohort%pft,currentCohort%c_area)

            if(arealayer > currentPatch%area.and.currentCohort%canopy_layer == i_lyr)then
               if (ED_val_comp_excln .ge. 0.0_r8 ) then
                  ! normal (stochastic) case. weight cohort demotion by inverse size to a constant power
                  currentCohort%excl_weight = 1.0_r8/(currentCohort%dbh**ED_val_comp_excln)  
               else
                  ! deterministic ranking case. only demote cohorts in smallest size classes
                  if ( (rankordered_area_sofar + currentCohort%c_area) .gt. currentPatch%area ) then
                     currentCohort%excl_weight = min(currentCohort%c_area, &
                           rankordered_area_sofar + currentCohort%c_area - currentPatch%area)
                  else
                     currentCohort%excl_weight = 0.0_r8
                  endif
                  rankordered_area_sofar = rankordered_area_sofar + currentCohort%c_area
               endif
               sumdiff = sumdiff + currentCohort%excl_weight
            endif
            currentCohort => currentCohort%shorter  
         enddo !currentCohort


         lossarea = arealayer - currentPatch%area  !how much do we have to lose?
         sum_weights = 0.0_r8
         currentCohort => currentPatch%tallest    !start from the tallest cohort

         ! Correct the demoted cohorts for  
         if (ED_val_comp_excln .ge. 0.0_r8 ) then
            do while (associated(currentCohort))
               if(currentCohort%canopy_layer  ==  i_lyr) then
                  weight = currentCohort%excl_weight/sumdiff
                  currentCohort%excl_weight = min(currentCohort%c_area/lossarea, weight)
                  sum_weights = sum_weights + currentCohort%excl_weight
               endif
               currentCohort => currentCohort%shorter      
            enddo
         endif

         currentCohort => currentPatch%tallest
         do while (associated(currentCohort))      
            if(currentCohort%canopy_layer == i_lyr)then !All the trees in this layer need to lose some area...
               if (ED_val_comp_excln .ge. 0.0_r8 ) then
                  weight = currentCohort%excl_weight/sum_weights
                  cc_loss = lossarea*weight !what this cohort has to lose. 
               else
                  ! in deterministic ranking mode, cohort loss is not renormalized
                  cc_loss = currentCohort%excl_weight
               endif
               if (cc_loss > 0._r8) then
                  !-----------Split and copy boundary cohort-----------------!
                  if(cc_loss < currentCohort%c_area)then


                     ! Make a copy of the current cohort.  The copy and the original
                     ! conserve total number density of the original.  The copy
                     ! remains in the upper-story.  The original is the one
                     ! demoted to the understory

                     ! n.b this needs to happen BEFORE the cohort goes into the new layer, 
                     ! otherwise currentPatch%spread(i_lyr+1) will be higher and the area will change...!!! 

                     allocate(copyc)
                     call copy_cohort(currentCohort, copyc) !

                     newarea = currentCohort%c_area - cc_loss
                     copyc%n = currentCohort%n*newarea/currentCohort%c_area   !
                     currentCohort%n = currentCohort%n - (currentCohort%n*newarea/currentCohort%c_area) !     

                     copyc%canopy_layer = i_lyr !the taller cohort is the copy

                     ! Demote the current cohort to the understory.
                     currentCohort%canopy_layer = i_lyr + 1 

                     ! keep track of number and biomass of demoted cohort
                     currentSite%demotion_rate(currentCohort%size_class) = &
                           currentSite%demotion_rate(currentCohort%size_class) + currentCohort%n
                     currentSite%demotion_carbonflux = currentSite%demotion_carbonflux + &
                           currentCohort%b_total() * currentCohort%n

                     !kill the ones which go into canopy layers that are not allowed... (default nclmax=2) 
                     if(i_lyr+1 > nclmax)then 
                        
                        !put the litter from the terminated cohorts into the fragmenting pools
                        do i_cwd=1,ncwd
                           
                           currentPatch%CWD_AG(i_cwd)  = currentPatch%CWD_AG(i_cwd) + &
                                 (currentCohort%bdead+currentCohort%bsw) * &
                                 EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * &
                                 SF_val_CWD_frac(i_cwd)*currentCohort%n/currentPatch%area  
                           
                           currentPatch%CWD_BG(i_cwd)  = currentPatch%CWD_BG(i_cwd) + &
                                 (currentCohort%bdead+currentCohort%bsw) * &
                                 (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) * &
                                 SF_val_CWD_frac(i_cwd)*currentCohort%n/currentPatch%area !litter flux per m2.
                           
                        enddo
                        
                        currentPatch%leaf_litter(currentCohort%pft)  = &
                              currentPatch%leaf_litter(currentCohort%pft) + (currentCohort%bl)* &
                              currentCohort%n/currentPatch%area ! leaf litter flux per m2.
                        
                        currentPatch%root_litter(currentCohort%pft)  = &
                              currentPatch%root_litter(currentCohort%pft) + &
                              (currentCohort%br+currentCohort%bstore)*currentCohort%n/currentPatch%area
                        
                        ! keep track of the above fluxes at the site level as a CWD/litter input flux (in kg / site-m2 / yr)
                        do i_cwd=1,ncwd
                           currentSite%CWD_AG_diagnostic_input_carbonflux(i_cwd) = &
                                 currentSite%CWD_AG_diagnostic_input_carbonflux(i_cwd) &
                                 + currentCohort%n*(currentCohort%bdead+currentCohort%bsw) * &
                                 SF_val_CWD_frac(i_cwd) * EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) &
                                 * hlm_days_per_year / AREA
                           currentSite%CWD_BG_diagnostic_input_carbonflux(i_cwd) = &
                                 currentSite%CWD_BG_diagnostic_input_carbonflux(i_cwd) &
                                 + currentCohort%n*(currentCohort%bdead+currentCohort%bsw) * &
                                 SF_val_CWD_frac(i_cwd) * (1.0_r8 -  &
                                 EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) * hlm_days_per_year / AREA
                        enddo

                        currentSite%leaf_litter_diagnostic_input_carbonflux(currentCohort%pft) = &
                              currentSite%leaf_litter_diagnostic_input_carbonflux(currentCohort%pft) +  &
                              currentCohort%n * (currentCohort%bl) * hlm_days_per_year  / AREA
                        currentSite%root_litter_diagnostic_input_carbonflux(currentCohort%pft) = &
                              currentSite%root_litter_diagnostic_input_carbonflux(currentCohort%pft) + &
                              currentCohort%n * (currentCohort%br+currentCohort%bstore) * hlm_days_per_year  / AREA

                        currentCohort%n = 0.0_r8
                        currentCohort%c_area = 0._r8
                     else  
                        call carea_allom(currentCohort%dbh,currentCohort%n, &
                              currentSite%spread,currentCohort%pft,currentCohort%c_area)
                     endif
                     
                     call carea_allom(copyc%dbh,copyc%n,currentSite%spread,copyc%pft,copyc%c_area)


                     !----------- Insert copy into linked list ------------------------!                         
                     copyc%shorter => currentCohort
                     if(associated(currentCohort%taller))then
                        copyc%taller => currentCohort%taller
                        currentCohort%taller%shorter => copyc
                     else
                        currentPatch%tallest => copyc
                        copyc%taller => null()
                     endif
                     currentCohort%taller => copyc

                  else ! matches: if(cc_loss < currentCohort%c_area)then

                     currentCohort%canopy_layer = i_lyr + 1 !the whole cohort becomes demoted

                     ! keep track of number and biomass of demoted cohort
                     currentSite%demotion_rate(currentCohort%size_class) = &
                           currentSite%demotion_rate(currentCohort%size_class) + currentCohort%n
                     currentSite%demotion_carbonflux = currentSite%demotion_carbonflux + &
                           currentCohort%b_total() * currentCohort%n

                     !kill the ones which go into canopy layers that are not allowed... (default nclmax=2) 
                     if(i_lyr+1 > nclmax)then  

                        !put the litter from the terminated cohorts into the fragmenting pools
                        do i_cwd=1,ncwd

                           currentPatch%CWD_AG(i_cwd)  = currentPatch%CWD_AG(i_cwd) + &
                                 (currentCohort%bdead+currentCohort%bsw) * &
                                 EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) * &
                                 SF_val_CWD_frac(i_cwd)*currentCohort%n/currentPatch%area           
                           currentPatch%CWD_BG(i_cwd)  = currentPatch%CWD_BG(i_cwd) + &
                                 (currentCohort%bdead+currentCohort%bsw) * &
                                 (1.0_r8-EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) * &
                                 SF_val_CWD_frac(i_cwd)*currentCohort%n/currentPatch%area !litter flux per m2.

                        enddo

                        currentPatch%leaf_litter(currentCohort%pft)  = &
                              currentPatch%leaf_litter(currentCohort%pft) + currentCohort%bl* &
                              currentCohort%n/currentPatch%area ! leaf litter flux per m2.

                        currentPatch%root_litter(currentCohort%pft)  = &
                              currentPatch%root_litter(currentCohort%pft) + &
                              (currentCohort%br+currentCohort%bstore)*currentCohort%n/currentPatch%area

                        ! keep track of the above fluxes at the site level as a CWD/litter input flux (in kg / site-m2 / yr)
                        do i_cwd=1,ncwd
                           currentSite%CWD_AG_diagnostic_input_carbonflux(i_cwd) = &
                                 currentSite%CWD_AG_diagnostic_input_carbonflux(i_cwd) &
                                 + currentCohort%n*(currentCohort%bdead+currentCohort%bsw) * &
                                 SF_val_CWD_frac(i_cwd) * EDPftvarcon_inst%allom_agb_frac(currentCohort%pft) &
                                 * hlm_days_per_year / AREA
                           currentSite%CWD_BG_diagnostic_input_carbonflux(i_cwd) = &
                                 currentSite%CWD_BG_diagnostic_input_carbonflux(i_cwd) &
                                 + currentCohort%n*(currentCohort%bdead+currentCohort%bsw) * &
                                 SF_val_CWD_frac(i_cwd) * (1.0_r8 -  EDPftvarcon_inst%allom_agb_frac(currentCohort%pft)) &
                                 * hlm_days_per_year / AREA
                        enddo

                        currentSite%leaf_litter_diagnostic_input_carbonflux(currentCohort%pft) = &
                              currentSite%leaf_litter_diagnostic_input_carbonflux(currentCohort%pft) +  &
                              currentCohort%n * (currentCohort%bl) * hlm_days_per_year  / AREA
                        currentSite%root_litter_diagnostic_input_carbonflux(currentCohort%pft) = &
                              currentSite%root_litter_diagnostic_input_carbonflux(currentCohort%pft) + &
                              currentCohort%n * (currentCohort%br+currentCohort%bstore) * hlm_days_per_year  / AREA

                        currentCohort%n      = 0.0_r8
                        currentCohort%c_area = 0._r8

                     else  
                        call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread,currentCohort%pft,currentCohort%c_area)
                     endif

                  endif ! matches: if (cc_loss < currentCohort%c_area)then
               endif    ! matches: if (cc_loss > 0._r8) then

               !----------- End of cohort splitting ------------------------------!             
            endif !canopy layer = i

            currentCohort => currentCohort%shorter

         enddo !currentCohort 


         ! Update the area calculations of the current layer
         ! And the layer below that may or may not had recieved
         ! Demotions

         call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer)

         layer_area_counter=layer_area_counter+1
         if(layer_area_counter > max_layer_iterations) then
            write(fates_log(),*) 'Layer demotion area not closing,i_lyr: ',i_lyr
            write(fates_log(),*) 'patch area:',currentpatch%area
            write(fates_log(),*) 'lat:',currentSite%lat
            write(fates_log(),*) 'lon:',currentSite%lon
            write(fates_log(),*) 'arealayer:',arealayer
            currentCohort => currentPatch%tallest
            do while (associated(currentCohort))  
               write(fates_log(),*) '---------------'
               write(fates_log(),*) 'coh ilayer:',currentCohort%canopy_layer
               write(fates_log(),*) 'coh dbh:',currentCohort%dbh
               write(fates_log(),*) 'coh pft:',currentCohort%pft
               write(fates_log(),*) 'coh n:',currentCohort%n
               write(fates_log(),*) 'coh carea:',currentCohort%c_area
               currentCohort => currentCohort%shorter
            enddo
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if


      enddo ! matches do while( (arealayer-currentPatch%area) > area_trim_tolerance ) 



   end subroutine DemoteFromLayer

   ! ==============================================================================================

   subroutine PromoteIntoLayer(currentSite,currentPatch,i_lyr)
      
      ! -------------------------------------------------------------------------------------------
      ! Check whether the intended 'full' layers are actually filling all the space.
      ! If not, promote some fraction of cohorts upwards.
      ! THIS SECTION MIGHT BE TRIGGERED BY A FIRE OR MORTALITY EVENT, FOLLOWED BY A PATCH FUSION, 
      ! SO THE TOP LAYER IS NO LONGER FULL.
      ! -------------------------------------------------------------------------------------------
      
      use EDParamsMod, only : ED_val_comp_excln

      ! !ARGUMENTS
      type(ed_site_type), intent(inout), target  :: currentSite
      type(ed_patch_type), intent(inout), target :: currentPatch
      integer, intent(in)                        :: i_lyr   ! Current canopy layer of interest

      ! !LOCAL VARIABLES:
      type(ed_cohort_type), pointer :: currentCohort
      type(ed_cohort_type), pointer :: copyc

      real(r8) :: sumdiff
      real(r8) :: promarea
      real(r8) :: newarea
      real(r8) :: sum_weights
      real(r8) :: weight
      real(r8) :: cc_gain
      real(r8) :: arealayer_current      ! area (m2) of the current canopy layer
      real(r8) :: arealayer_below        ! area (m2) of the layer below the current layer
      real(r8) :: rankordered_area_sofar ! the amount of total canopy area occupied by cohorts upto this point
      integer  :: layer_area_counter
      logical  :: layer_below_exists     ! If enough of the layer below exists

      
      

      call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer_current)
      call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr+1,arealayer_below)


      layer_area_counter = 0
      layer_below_exists = .true.
      do while( (arealayer_current-currentPatch%area) < -area_target_precision .and. layer_below_exists )
         

         ! Promote all cohorts from layer below if that whole layer has area smaller
         ! than the tolerance on the gains needed into current layer
         ! -------------------------------------------------------------------------------------

         if(arealayer_below <= area_target_precision)then
            currentCohort => currentPatch%tallest 
            do while (associated(currentCohort))            
               if(currentCohort%canopy_layer == i_lyr+1)then !look at the cohorts in the canopy layer below... 
                  currentCohort%canopy_layer = i_lyr   
                  call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread, &
                        currentCohort%pft,currentCohort%c_area)
                  ! keep track of number and biomass of promoted cohort
                  currentSite%promotion_rate(currentCohort%size_class) = &
                        currentSite%promotion_rate(currentCohort%size_class) + currentCohort%n
                  currentSite%promotion_carbonflux = currentSite%promotion_carbonflux + &
                        currentCohort%b_total() * currentCohort%n

               endif
               currentCohort => currentCohort%shorter   
            enddo
            arealayer_below = 0.0_r8
            call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer_current)
         endif                 
         
         sumdiff = 0.0_r8    
         rankordered_area_sofar = 0.0_r8
         ! figure out with what weighting we need to promote cohorts.
         ! This is the opposite of the demotion weighting... 
         currentCohort => currentPatch%tallest 
         do while (associated(currentCohort))
            call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread, &
                  currentCohort%pft,currentCohort%c_area)
            if(currentCohort%canopy_layer == i_lyr+1)then !look at the cohorts in the canopy layer below... 
               if (ED_val_comp_excln .ge. 0.0_r8 ) then
                  ! normal (stochastic) case, as above.
                  currentCohort%prom_weight = currentCohort%dbh**ED_val_comp_excln   !as opposed to 1/(dbh^C_e) 
               else
                  ! deterministic case, as above, but inverse, so only take tallest cohorts from i+1 canopy layer
                  if ( rankordered_area_sofar .lt. currentPatch%area - arealayer_current  ) then
                     currentCohort%prom_weight = max(min(currentCohort%c_area, &
                           currentPatch%area - arealayer_current - rankordered_area_sofar ), 0._r8)
                  else
                     currentCohort%prom_weight = 0.0_r8
                  endif
                  rankordered_area_sofar = rankordered_area_sofar + currentCohort%c_area
               endif
               sumdiff = sumdiff + currentCohort%prom_weight
            endif
            currentCohort => currentCohort%shorter  
         enddo !currentCohort


         promarea    =  currentPatch%area - arealayer_current ! how much do we need to gain?
         sum_weights = 0.0_r8

         if (ED_val_comp_excln .ge. 0.0_r8 ) then
            currentCohort => currentPatch%tallest    !start from the tallest cohort
            do while (associated(currentCohort))
               if(currentCohort%canopy_layer  ==  i_lyr+1) then !still looking at the layer beneath. 
                  weight = currentCohort%prom_weight/sumdiff
                  if(promarea > 0._r8)then    
                     currentCohort%prom_weight = min(currentCohort%c_area/promarea, weight)
                  else
                     currentCohort%prom_weight = 0._r8
                  endif
                  sum_weights = sum_weights + currentCohort%prom_weight
               endif
               currentCohort => currentCohort%shorter      
            enddo
         endif
         
         currentCohort => currentPatch%tallest
         do while (associated(currentCohort))      
            if(currentCohort%canopy_layer == i_lyr+1)then !All the trees in this layer need to promote some area upwards... 
               
               if (ED_val_comp_excln .ge. 0.0_r8) then
                  ! normal mode, renormalize areas
                  weight = currentCohort%prom_weight/sum_weights
                  cc_gain = promarea*weight !what this cohort has to promote. 
               else
                  ! in deterministic ranking mode, cohort loss is not renormalized   
                  cc_gain = currentCohort%prom_weight
               endif
               if ( cc_gain > 0._r8 ) then
                  !-----------Split and copy boundary cohort-----------------!
                  if(cc_gain < currentCohort%c_area)then
                     allocate(copyc)
                     
                     call copy_cohort(currentCohort, copyc) !makes an identical copy...
                     ! n.b this needs to happen BEFORE the cohort goes into the new layer, otherwise currentPatch
                     ! %spread(+1) will be higher and the area will change...!!!
                     
                     newarea = currentCohort%c_area - cc_gain !new area of existing cohort
                     copyc%n = currentCohort%n*cc_gain/currentCohort%c_area   !number of individuals in promoted cohort. 
                     ! number of individuals in cohort remaining in understorey    
                     currentCohort%n = currentCohort%n - (currentCohort%n*cc_gain/currentCohort%c_area) 
                     
                     currentCohort%canopy_layer = i_lyr + 1 ! keep current cohort in the understory.        
                     copyc%canopy_layer = i_lyr             ! promote copy to the higher canopy layer. 
                     
                     ! keep track of number and biomass of promoted cohort
                     currentSite%promotion_rate(copyc%size_class) = &
                           currentSite%promotion_rate(copyc%size_class) + copyc%n
                     currentSite%promotion_carbonflux = currentSite%promotion_carbonflux + &
                           copyc%b_total() * copyc%n
                         
                     ! seperate cohorts. 
                     ! needs to be a very small number to avoid causing non-linearity issues with c_area. 
                     ! is this really required? 
                     currentCohort%dbh = currentCohort%dbh - 0.000000000001_r8 
                     copyc%dbh = copyc%dbh + 0.000000000001_r8
                     
                     call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread, &
                           currentCohort%pft,currentCohort%c_area)
                     call carea_allom(copyc%dbh,copyc%n,currentSite%spread,copyc%pft,copyc%c_area)
                         
                     !----------- Insert copy into linked list ------------------------!                         
                     copyc%shorter => currentCohort
                     if(associated(currentCohort%taller))then
                        copyc%taller => currentCohort%taller
                        currentCohort%taller%shorter => copyc
                     else
                        currentPatch%tallest => copyc
                        copyc%taller => null()
                     endif
                     currentCohort%taller => copyc                  
                  else             ! if(cc_gain < currentCohort%c_area)then
                     currentCohort%canopy_layer = i_lyr  !the whole cohort becomes promoted
                     
                     ! update area AFTER we sum up the losses. the cohort may shrink at this point,
                     ! if the upper canopy spread is smaller. this shold be dealt with by the 'excess area' loop.  
                     call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread, &
                           currentCohort%pft,currentCohort%c_area)

                     ! keep track of number and biomass of promoted cohort
                     currentSite%promotion_rate(currentCohort%size_class) = &
                           currentSite%promotion_rate(currentCohort%size_class) + currentCohort%n
                     currentSite%promotion_carbonflux = currentSite%promotion_carbonflux + &
                           currentCohort%b_total() * currentCohort%n
                     
                  endif          ! if(cc_gain < currentCohort%c_area)then

               endif             ! if ( cc_gain > 0._r8 ) then
               
            endif   ! if(currentCohort%canopy_layer == i_lyr+1) then
            currentCohort => currentCohort%shorter
         enddo !currentCohort 

         call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer_current)
         call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr+1,arealayer_below)

         ! Only continue trying to promote if
         ! there is enough canopy area in the layer below
         ! to promote cohorts.  Make sure it is not just more than zero,
         ! but larger than the precision in the area comparison criteria
         if( arealayer_below > area_target_precision) then
            layer_below_exists = .true.
         else
            layer_below_exists = .false.
         end if

         layer_area_counter = layer_area_counter + 1
         if(layer_area_counter > max_layer_iterations) then
            write(fates_log(),*) 'Layer promotion area not closing,i_lyr: ',i_lyr
            write(fates_log(),*) 'patch area:',currentpatch%area
            write(fates_log(),*) 'lat:',currentSite%lat
            write(fates_log(),*) 'lon:',currentSite%lon
            write(fates_log(),*) 'arealayer_current:',arealayer_current
            write(fates_log(),*) 'arealayer_below:',arealayer_below
            currentCohort => currentPatch%tallest
            do while (associated(currentCohort))  
               write(fates_log(),*) '---------------'
               write(fates_log(),*) 'coh ilayer:',currentCohort%canopy_layer
               write(fates_log(),*) 'coh dbh:',currentCohort%dbh
               write(fates_log(),*) 'coh pft:',currentCohort%pft
               write(fates_log(),*) 'coh n:',currentCohort%n
               write(fates_log(),*) 'coh carea:',currentCohort%c_area
               currentCohort => currentCohort%shorter
            enddo
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
         
      enddo ! do while( (arealayer_current-currentPatch%area) 
      !            < -0.000001_r8 .and. layer_below_exists ) then
      
      return
   end subroutine PromoteIntoLayer

  ! ============================================================================

  subroutine canopy_spread( currentSite )
    !
    ! !DESCRIPTION:
    !  Calculates the spatial spread of tree canopies based on canopy closure.                             
    !
    ! !USES:
    use EDTypesMod        , only : AREA
    use EDParamsMod, only : ED_val_canopy_closure_thresh    
    !
    ! !ARGUMENTS    
    type (ed_site_type), intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type), pointer :: currentCohort
    type (ed_patch_type) , pointer :: currentPatch
    real(r8) :: sitelevel_canopyarea  ! Amount of canopy in top layer at the site level
    real(r8) :: inc                   ! Arbitrary daily incremental change in canopy area 
    integer  :: z
    !----------------------------------------------------------------------

    inc = 0.05_r8

    currentPatch => currentSite%oldest_patch

    sitelevel_canopyarea = 0.0_r8   
    do while (associated(currentPatch))

       !calculate canopy area in each patch...
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))
          call carea_allom(currentCohort%dbh,currentCohort%n,currentSite%spread,currentCohort%pft,currentCohort%c_area)
          if(EDPftvarcon_inst%woody(currentCohort%pft) .eq. 1 .and. currentCohort%canopy_layer .eq. 1 ) then
             sitelevel_canopyarea = sitelevel_canopyarea + currentCohort%c_area
          endif
          currentCohort => currentCohort%shorter
       enddo

       currentPatch => currentPatch%younger

    enddo !currentPatch

    !If the canopy area is approaching closure, squash the tree canopies and make them taller and thinner
    if( sitelevel_canopyarea/AREA .gt. ED_val_canopy_closure_thresh ) then
       currentSite%spread = currentSite%spread - inc
    else 
       currentSite%spread = currentSite%spread + inc 
    endif

    ! put within bounds to make sure it stays between 0 and 1
    currentSite%spread = max(min(currentSite%spread, 1._r8), 0._r8)

  end subroutine canopy_spread


  ! =====================================================================================

  subroutine canopy_summarization( nsites, sites, bc_in )

     ! ----------------------------------------------------------------------------------
     ! Much of this routine was once ed_clm_link minus all the IO and history stuff
     ! ---------------------------------------------------------------------------------

    use FatesInterfaceMod    , only : bc_in_type
    use EDPatchDynamicsMod   , only : set_patchno
    use EDPatchDynamicsMod   , only : set_root_fraction
    use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
    use EDtypesMod           , only : area
    use EDPftvarcon            , only : EDPftvarcon_inst

    ! !ARGUMENTS    
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort
    integer  :: s
    integer  :: ft                                      ! plant functional type
    integer  :: ifp
    integer  :: patchn                                  ! identification number for each patch. 
    real(r8) :: canopy_leaf_area                        ! total amount of leaf area in the vegetated area. m2.  

    !----------------------------------------------------------------------

    if ( DEBUG ) then
       write(fates_log(),*) 'in canopy_summarization'
    endif

    do s = 1,nsites
       
       ! --------------------------------------------------------------------------------
       ! Set the patch indices (this is usefull mostly for communicating with a host or 
       ! driving model.  Loops through all patches and sets cpatch%patchno to the integer 
       ! order of oldest to youngest where the oldest is 1.
       ! --------------------------------------------------------------------------------
       call set_patchno( sites(s) )

       currentPatch => sites(s)%oldest_patch

       do while(associated(currentPatch))
          
          call set_root_fraction(currentPatch,bc_in(s)%zi_sisl)

          !zero cohort-summed variables. 
          currentPatch%total_canopy_area = 0.0_r8
          currentPatch%total_tree_area = 0.0_r8
          canopy_leaf_area = 0.0_r8
          
          !update cohort quantitie s                                  
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             
             ft = currentCohort%pft

             
             ! Update the cohort's index within the size bin classes
             ! Update the cohort's index within the SCPF classification system
             call sizetype_class_index(currentCohort%dbh,currentCohort%pft, &
                                       currentCohort%size_class,currentCohort%size_by_pft_class)

             call carea_allom(currentCohort%dbh,currentCohort%n,sites(s)%spread,&
                  currentCohort%pft,currentCohort%c_area)
             currentCohort%treelai = tree_lai(currentCohort%bl, currentCohort%status_coh, &
                  currentCohort%pft, currentCohort%c_area, currentCohort%n )

             canopy_leaf_area = canopy_leaf_area + currentCohort%treelai *currentCohort%c_area
                  
             if(currentCohort%canopy_layer==1)then
                currentPatch%total_canopy_area = currentPatch%total_canopy_area + currentCohort%c_area
                if(EDPftvarcon_inst%woody(ft)==1)then
                   currentPatch%total_tree_area = currentPatch%total_tree_area + currentCohort%c_area
                endif
             endif
             
             ! Check for erroneous zero values. 
             if(currentCohort%dbh <= 0._r8 .or. currentCohort%n == 0._r8)then
                write(fates_log(),*) 'FATES: dbh or n is zero in canopy_summarization', &
                      currentCohort%dbh,currentCohort%n
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif
             if(currentCohort%pft == 0.or.currentCohort%canopy_trim <= 0._r8)then
                write(fates_log(),*) 'FATES: PFT or trim is zero in canopy_summarization', &
                      currentCohort%pft,currentCohort%canopy_trim
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif
             if( (currentCohort%bsw + currentCohort%bl + currentCohort%br) <= 0._r8)then
                write(fates_log(),*) 'FATES: alive biomass is zero in canopy_summarization', &
                      currentCohort%bsw + currentCohort%bl + currentCohort%br
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif

             currentCohort => currentCohort%taller
             
          enddo ! ends 'do while(associated(currentCohort))
          
          if ( currentPatch%total_canopy_area-currentPatch%area > 0.000001_r8 ) then
             if ( currentPatch%total_canopy_area-currentPatch%area > 0.001_r8 ) then
                write(fates_log(),*) 'FATES: canopy area bigger than area', &
                     currentPatch%total_canopy_area ,currentPatch%area
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             currentPatch%total_canopy_area = currentPatch%area
          endif

          currentPatch => currentPatch%younger
       end do !patch loop
            
       call leaf_area_profile(sites(s),bc_in(s)%snow_depth_si,bc_in(s)%frac_sno_eff_si) 
       
    end do ! site loop
    
    return
  end subroutine canopy_summarization
 
 ! =====================================================================================

 subroutine leaf_area_profile( currentSite , snow_depth_si, frac_sno_eff_si)
    
    ! -----------------------------------------------------------------------------------
    ! This subroutine calculates how leaf and stem areas are distributed 
    ! in vertical and horizontal space.
    !
    ! The following cohort level diagnostics are updated here:
    ! 
    ! currentCohort%treelai    ! LAI per unit crown area  (m2/m2)
    ! currentCohort%treesai    ! SAI per unit crown area  (m2/m2)
    ! currentCohort%lai        ! LAI per unit canopy area (m2/m2)
    ! currentCohort%sai        ! SAI per unit canopy area (m2/m2)
    ! currentCohort%NV         ! The number of discrete vegetation
    !                          ! layers needed to describe this crown
    !
    ! The following patch level diagnostics are updated here:
    ! 
    ! currentPatch%canopy_layer_tai(cl)    ! TAI of each canopy layer
    ! currentPatch%ncan(cl,ft)             ! number of vegetation layers needed
    !                                      ! in this patch's pft/canopy-layer 
    ! currentPatch%nrad(cl,ft)             ! same as ncan, but does not include
    !                                      ! layers occluded by snow
    !                                      ! CURRENTLY SAME AS NCAN
    ! currentPatch%canopy_mask(cl,ft)      ! are there canopy elements in this pft-layer?
    !                                      ! (This is redundant with nrad though...)
    ! currentPatch%tlai_profile(cl,ft,iv)  ! m2 of leaves per m2 of the PFT's footprint
    ! currentPatch%elai_profile(cl,ft,iv)  ! non-snow covered m2 of leaves per m2 of PFT footprint
    ! currentPatch%tsai_profile(cl,ft,iv)  ! m2 of stems per m2 of PFT footprint
    ! currentPatch%esai_profile(cl,ft,iv)  ! non-snow covered m2 of stems per m2 of PFT footprint
    ! currentPatch%canopy_area_profile(cl,ft,iv)  ! Fractional area of leaf layer 
    !                                             ! relative to vegetated area
    ! currentPatch%layer_height_profile(cl,ft,iv) ! Elevation of layer in m
    !
    ! -----------------------------------------------------------------------------------

    ! !USES:

    use EDtypesMod           , only : area, dinc_ed, hitemax, n_hite_bins
  
    !
    ! !ARGUMENTS    
    type(ed_site_type)     , intent(inout) :: currentSite
    real(r8)               , intent(in)    :: snow_depth_si
    real(r8)               , intent(in)    :: frac_sno_eff_si

    !
    ! !LOCAL VARIABLES:
    type (ed_patch_type)  , pointer :: currentPatch
    type (ed_cohort_type) , pointer :: currentCohort
    real(r8) :: remainder                !Thickness of layer at bottom of canopy. 
    real(r8) :: fleaf                    ! fraction of cohort incepting area that is leaves.  
    integer  :: ft                       ! Plant functional type index. 
    integer  :: iv                       ! Vertical leaf layer index   
    integer  :: cl                       ! Canopy layer index
    real(r8) :: fraction_exposed         ! how much of this layer is not covered by snow?
    real(r8) :: layer_top_hite           ! notional top height of this canopy layer (m)
    real(r8) :: layer_bottom_hite        ! notional bottom height of this canopy layer (m)
    integer  :: smooth_leaf_distribution ! is the leaf distribution this option (1) or not (0)
    real(r8) :: frac_canopy(N_HITE_BINS) ! amount of canopy in each height class
    real(r8) :: patch_lai                ! LAI summed over the patch in m2/m2 of canopy area
    real(r8) :: minh(N_HITE_BINS)        ! minimum height in height class (m)
    real(r8) :: maxh(N_HITE_BINS)        ! maximum height in height class (m)
    real(r8) :: dh                       ! vertical detph of height class (m)
    real(r8) :: min_chite                ! bottom of cohort canopy  (m)
    real(r8) :: max_chite                ! top of cohort canopy      (m)
    real(r8) :: lai                      ! summed lai for checking m2 m-2
    real(r8) :: snow_depth_avg           ! avg snow over whole site
    
    !----------------------------------------------------------------------



    smooth_leaf_distribution = 0

    ! Here we are trying to generate a profile of leaf area, indexed by 'z' and by pft
    ! We assume that each point in the canopy recieved the light attenuated by the average
    ! leaf area index above it, irrespective of PFT identity... 
    ! Each leaf is defined by how deep in the canopy it is, in terms of LAI units.  (FIX(RF,032414), GB)
    
    currentPatch => currentSite%oldest_patch   
    do while(associated(currentPatch))

       ! --------------------------------------------------------------------------------
       ! Calculate tree and canopy areas. 
       ! calculate tree lai and sai.
       ! --------------------------------------------------------------------------------

       currentPatch%canopy_layer_tai(:)         = 0._r8
       currentPatch%ncan(:,:)                   = 0 
       currentPatch%nrad(:,:)                   = 0 
       patch_lai                                = 0._r8
       currentPatch%tlai_profile(:,:,:)         = 0._r8
       currentPatch%tsai_profile(:,:,:)         = 0._r8  
       currentPatch%elai_profile(:,:,:)         = 0._r8
       currentPatch%esai_profile(:,:,:)         = 0._r8 
       currentPatch%layer_height_profile(:,:,:) = 0._r8
       currentPatch%canopy_area_profile(:,:,:)  = 0._r8       
       currentPatch%canopy_mask(:,:)            = 0

       ! ------------------------------------------------------------------------------
       ! It is remotely possible that in deserts we will not have any canopy
       ! area, ie not plants at all...
       ! ------------------------------------------------------------------------------
       
       if (currentPatch%total_canopy_area > tiny(currentPatch%total_canopy_area)) then

       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 

          ft = currentCohort%pft
          cl = currentCohort%canopy_layer

          currentCohort%treelai = tree_lai(currentCohort%bl, currentCohort%status_coh, currentCohort%pft, &
               currentCohort%c_area, currentCohort%n )
          currentCohort%treesai = tree_sai(currentCohort%dbh, currentCohort%pft, currentCohort%canopy_trim, &
               currentCohort%c_area, currentCohort%n)

          currentCohort%lai =  currentCohort%treelai *currentCohort%c_area/currentPatch%total_canopy_area 
          currentCohort%sai =  currentCohort%treesai *currentCohort%c_area/currentPatch%total_canopy_area  

          ! Number of actual vegetation layers in this cohort's crown
          currentCohort%NV =  ceiling((currentCohort%treelai+currentCohort%treesai)/dinc_ed)  

          currentPatch%ncan(cl,ft) = max(currentPatch%ncan(cl,ft),currentCohort%NV)

          patch_lai = patch_lai + currentCohort%lai

!          currentPatch%canopy_layer_tai(cl) = currentPatch%canopy_layer_tai(cl) + &
!                currentCohort%lai + currentCohort%sai

          do cl = 1,nclmax-1
             if(currentCohort%canopy_layer == cl)then
                currentPatch%canopy_layer_tai(cl) = currentPatch%canopy_layer_tai(cl) + &
                     currentCohort%lai + currentCohort%sai
             endif
          enddo

          currentCohort => currentCohort%taller 
          
       enddo !currentCohort

       if(smooth_leaf_distribution == 1)then

          ! -----------------------------------------------------------------------------
          ! we are going to ignore the concept of canopy layers, and put all of the leaf 
          ! area into height banded bins.  using the same domains as we had before, except 
          ! that CL always = 1
          ! -----------------------------------------------------------------------------
          
          ! this is a crude way of dividing up the bins. Should it be a function of actual maximum height? 
          dh = 1.0_r8*(HITEMAX/N_HITE_BINS) 
          do iv = 1,N_HITE_BINS  
             if (iv == 1) then
                minh(iv) = 0.0_r8
                maxh(iv) = dh
             else 
                minh(iv) = (iv-1)*dh
                maxh(iv) = (iv)*dh
             endif
          enddo
          
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))  
             ft = currentCohort%pft
             min_chite = currentCohort%hite - currentCohort%hite * EDPftvarcon_inst%crown(ft)
             max_chite = currentCohort%hite  
             do iv = 1,N_HITE_BINS  
                frac_canopy(iv) = 0.0_r8
                ! this layer is in the middle of the canopy
                if(max_chite > maxh(iv).and.min_chite < minh(iv))then 
                   frac_canopy(iv)= min(1.0_r8,dh / (currentCohort%hite*EDPftvarcon_inst%crown(ft)))
                   ! this is the layer with the bottom of the canopy in it. 
                elseif(min_chite < maxh(iv).and.min_chite > minh(iv).and.max_chite > maxh(iv))then 
                   frac_canopy(iv) = (maxh(iv) -min_chite ) / (currentCohort%hite*EDPftvarcon_inst%crown(ft))
                   ! this is the layer with the top of the canopy in it. 
                elseif(max_chite > minh(iv).and.max_chite < maxh(iv).and.min_chite < minh(iv))then 
                   frac_canopy(iv) = (max_chite - minh(iv)) / (currentCohort%hite*EDPftvarcon_inst%crown(ft))
                elseif(max_chite < maxh(iv).and.min_chite > minh(iv))then !the whole cohort is within this layer. 
                   frac_canopy(iv) = 1.0_r8
                endif
                
                ! no m2 of leaf per m2 of ground in each height class
                currentPatch%tlai_profile(1,ft,iv) = currentPatch%tlai_profile(1,ft,iv) + frac_canopy(iv) * &
                      currentCohort%lai
                currentPatch%tsai_profile(1,ft,iv) = currentPatch%tsai_profile(1,ft,iv) + frac_canopy(iv) * &
                      currentCohort%sai
                
                !snow burial
                !write(fates_log(), *) 'calc snow'
                snow_depth_avg = snow_depth_si * frac_sno_eff_si
                if(snow_depth_avg  > maxh(iv))then
                   fraction_exposed = 0._r8
                endif
                if(snow_depth_avg < minh(iv))then
                   fraction_exposed = 1._r8
                endif
                if(snow_depth_avg>= minh(iv).and.snow_depth_avg <= maxh(iv))then !only partly hidden... 
                   fraction_exposed =  max(0._r8,(min(1.0_r8,(snow_depth_avg-minh(iv))/dh)))
                endif
                fraction_exposed = 1.0_r8
                ! no m2 of leaf per m2 of ground in each height class
                ! FIX(SPM,032414) these should be uncommented this and double check
                
                if ( DEBUG ) write(fates_log(), *) 'leaf_area_profile()', currentPatch%elai_profile(1,ft,iv)
                
                currentPatch%elai_profile(1,ft,iv) = currentPatch%tlai_profile(1,ft,iv) * fraction_exposed
                currentPatch%esai_profile(1,ft,iv) = currentPatch%tsai_profile(1,ft,iv) * fraction_exposed
                
                if ( DEBUG ) write(fates_log(), *) 'leaf_area_profile()', currentPatch%elai_profile(1,ft,iv)
                
             enddo ! (iv) hite bins
             
             currentCohort => currentCohort%taller
             
          enddo !currentCohort 
          
          ! -----------------------------------------------------------------------------
          ! Perform a leaf area conservation check on the LAI profile
          lai = 0.0_r8
          do ft = 1,numpft
             lai = lai+ sum(currentPatch%tlai_profile(1,ft,:))
          enddo
          
          if(lai > patch_lai)then
             write(fates_log(), *) 'FATES: problem with lai assignments'
             call endrun(msg=errMsg(sourcefile, __LINE__))
          endif
          
          
       else ! smooth leaf distribution  

          ! -----------------------------------------------------------------------------
          ! Standard canopy layering model.
          ! Go through all cohorts and add their leaf area 
          ! and canopy area to the accumulators. 
          ! -----------------------------------------------------------------------------

             
             currentCohort => currentPatch%shortest
             do while(associated(currentCohort))   
                
                ft = currentCohort%pft 
                cl = currentCohort%canopy_layer
                
                ! ----------------------------------------------------------------
                ! How much of each tree is stem area index? Assuming that there is 
                ! This may indeed be zero if there is a sensecent grass
                ! ----------------------------------------------------------------
                
                if( (currentCohort%treelai+currentCohort%treesai) > 0._r8)then    
                   fleaf = currentCohort%lai / (currentCohort%lai + currentCohort%sai) 
                else
                   fleaf = 0._r8
                endif
                
                ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                ! SNOW BURIAL IS CURRENTLY TURNED OFF
                ! WHEN IT IS TURNED ON, IT WILL HAVE TO BE COMPARED
                ! WITH SNOW HEIGHTS CALCULATED BELOW.
                ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                
                currentPatch%nrad(cl,ft) = currentPatch%ncan(cl,ft) 

                if (currentPatch%nrad(cl,ft) > nlevleaf ) then
                   write(fates_log(), *) 'Number of radiative leaf layers is larger'
                   write(fates_log(), *) ' than the maximum allowed.'
                   write(fates_log(), *) ' cl: ',cl
                   write(fates_log(), *) ' ft: ',ft
                   write(fates_log(), *) ' nlevleaf: ',nlevleaf
                   write(fates_log(), *) ' currentPatch%nrad(cl,ft): ', currentPatch%nrad(cl,ft)
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if


                ! --------------------------------------------------------------------------
                ! Whole layers.  Make a weighted average of the leaf area in each layer 
                ! before dividing it by the total area. Fill up layer for whole layers.  
                ! --------------------------------------------------------------------------
                
                do iv = 1,currentCohort%NV
                   
                   ! This loop builds the arrays that define the effective (not snow covered)
                   ! and total (includes snow covered) area indices for leaves and stems
                   ! We calculate the absolute elevation of each layer to help determine if the layer
                   ! is obscured by snow.
                   
                   layer_top_hite = currentCohort%hite - &
                         ( dble(iv-1.0)/currentCohort%NV * currentCohort%hite *  &
                         EDPftvarcon_inst%crown(currentCohort%pft) )
                   
                   layer_bottom_hite = currentCohort%hite - &
                         ( dble(iv)/currentCohort%NV * currentCohort%hite * &
                         EDPftvarcon_inst%crown(currentCohort%pft) )
                   
                   fraction_exposed = 1.0_r8
                   snow_depth_avg = snow_depth_si * frac_sno_eff_si
                   if(snow_depth_avg  > layer_top_hite)then
                      fraction_exposed = 0._r8
                   endif
                   if(snow_depth_avg < layer_bottom_hite)then
                      fraction_exposed = 1._r8
                   endif
                   if( snow_depth_avg>= layer_bottom_hite .and. &
                         snow_depth_avg <= layer_top_hite) then !only partly hidden...
                      fraction_exposed =  max(0._r8,(min(1.0_r8,(snow_depth_avg-layer_bottom_hite)/ &
                         (layer_top_hite-layer_bottom_hite ))))
                   endif
                   
                   ! =========== OVER-WRITE =================
                   fraction_exposed= 1.0_r8
                   ! =========== OVER-WRITE =================
                   
                   if(iv==currentCohort%NV) then
                      remainder = (currentCohort%treelai + currentCohort%treesai) - &
                            (dinc_ed*dble(currentCohort%NV-1.0_r8))
                      if(remainder > dinc_ed )then
                         write(fates_log(), *)'ED: issue with remainder', &
                               currentCohort%treelai,currentCohort%treesai,dinc_ed, & 
                               currentCohort%NV,remainder
                         call endrun(msg=errMsg(sourcefile, __LINE__))
                      endif
                   else
                      remainder = dinc_ed
                   end if
                   
                   currentPatch%tlai_profile(cl,ft,iv) = currentPatch%tlai_profile(cl,ft,iv) + &
                         remainder * fleaf * currentCohort%c_area/currentPatch%total_canopy_area
                   
                   currentPatch%elai_profile(cl,ft,iv) = currentPatch%elai_profile(cl,ft,iv) + &
                         remainder * fleaf * currentCohort%c_area/currentPatch%total_canopy_area * &
                         fraction_exposed
                   
                   currentPatch%tsai_profile(cl,ft,iv) = currentPatch%tsai_profile(cl,ft,iv) + &
                         remainder * (1._r8 - fleaf) * currentCohort%c_area/currentPatch%total_canopy_area
                   
                   currentPatch%esai_profile(cl,ft,iv) = currentPatch%esai_profile(cl,ft,iv) + &
                         remainder * (1._r8 - fleaf) * currentCohort%c_area/currentPatch%total_canopy_area * &
                         fraction_exposed
                   
                   currentPatch%canopy_area_profile(cl,ft,iv) = currentPatch%canopy_area_profile(cl,ft,iv) + &
                         currentCohort%c_area/currentPatch%total_canopy_area
                   
                   currentPatch%layer_height_profile(cl,ft,iv) = currentPatch%layer_height_profile(cl,ft,iv) + &
                         (remainder * fleaf * currentCohort%c_area/currentPatch%total_canopy_area * &
                         (layer_top_hite+layer_bottom_hite)/2.0_r8) !average height of layer. 
                   
                end do
                
                currentCohort => currentCohort%taller
                
             enddo !cohort
         
             ! --------------------------------------------------------------------------
            
             ! If there is an upper-story, the top canopy layer
             ! should have a value of exactly 1.0 in its top leaf layer
             ! --------------------------------------------------------------------------
             
             if ( (currentPatch%NCL_p > 1) .and. &
                  (sum(currentPatch%canopy_area_profile(1,:,1)) < 0.9999 )) then
                write(fates_log(), *) 'FATES: canopy_area_profile was less than 1 at the canopy top'
                write(fates_log(), *) 'cl: ',1
                write(fates_log(), *) 'iv: ',1
                write(fates_log(), *) 'sum(cpatch%canopy_area_profile(1,:,1)): ', &
                     sum(currentPatch%canopy_area_profile(1,:,1))
                currentCohort => currentPatch%shortest
                do while(associated(currentCohort))
                   if(currentCohort%canopy_layer==1)then
                      write(fates_log(), *) 'FATES: cohorts',currentCohort%dbh,currentCohort%c_area, &
                           currentPatch%total_canopy_area,currentPatch%area
                      write(fates_log(), *) 'ED: fracarea', currentCohort%pft, &
                           currentCohort%c_area/currentPatch%total_canopy_area
                   endif
                   currentCohort => currentCohort%taller  
                enddo !currentCohort
                call endrun(msg=errMsg(sourcefile, __LINE__))
                
             end if
         

             ! --------------------------------------------------------------------------
             ! In the following loop we are now normalizing the effective and
             ! total area profiles to convert from units of leaf/stem area per vegetated
             ! canopy area, into leaf/stem area per area of their own radiative column
             ! which is typically the footprint of all cohorts contained in the canopy
             ! layer x pft bins.
             ! Also perform some checks on area normalization.
             ! Check the area of each leaf layer, across pfts.
             ! It should never be larger than 1 or less than 0.
             ! --------------------------------------------------------------------------

             do cl = 1,currentPatch%NCL_p
                do iv = 1,currentPatch%ncan(cl,ft)
                   
                   if( DEBUG .and. sum(currentPatch%canopy_area_profile(cl,:,iv)) > 1.0001_r8 ) then
                      
                      write(fates_log(), *) 'FATES: A canopy_area_profile exceeded 1.0'
                      write(fates_log(), *) 'cl: ',cl
                      write(fates_log(), *) 'iv: ',iv
                      write(fates_log(), *) 'sum(cpatch%canopy_area_profile(cl,:,iv)): ', &
                            sum(currentPatch%canopy_area_profile(cl,:,iv))
                       currentCohort => currentPatch%shortest
                       do while(associated(currentCohort))
                          if(currentCohort%canopy_layer==cl)then
                             write(fates_log(), *) 'FATES: cohorts in layer cl = ',cl, &
                                  currentCohort%dbh,currentCohort%c_area, &
                                  currentPatch%total_canopy_area,currentPatch%area
                             write(fates_log(), *) 'ED: fracarea', currentCohort%pft, &
                                  currentCohort%c_area/currentPatch%total_canopy_area
                          endif
                          currentCohort => currentCohort%taller  
                       enddo !currentCohort
                       call endrun(msg=errMsg(sourcefile, __LINE__))
                    end if
                end do
                   
                do ft = 1,numpft
                   do iv = 1,currentPatch%ncan(cl,ft)

                      if( currentPatch%canopy_area_profile(cl,ft,iv) > &
                          tiny(currentPatch%canopy_area_profile(cl,ft,iv)) )then
                         
                         currentPatch%tlai_profile(cl,ft,iv) = currentPatch%tlai_profile(cl,ft,iv) / &
                               currentPatch%canopy_area_profile(cl,ft,iv)
                         
                         currentPatch%tsai_profile(cl,ft,iv) = currentPatch%tsai_profile(cl,ft,iv) / &
                               currentPatch%canopy_area_profile(cl,ft,iv)
                         
                         currentPatch%elai_profile(cl,ft,iv) = currentPatch%elai_profile(cl,ft,iv) / &
                               currentPatch%canopy_area_profile(cl,ft,iv)
                         
                         currentPatch%esai_profile(cl,ft,iv) = currentPatch%esai_profile(cl,ft,iv) / &
                               currentPatch%canopy_area_profile(cl,ft,iv)
                      end if
                      
                      if(currentPatch%tlai_profile(cl,ft,iv)>tiny(currentPatch%tlai_profile(cl,ft,iv)))then
                         currentPatch%layer_height_profile(cl,ft,iv) = currentPatch%layer_height_profile(cl,ft,iv) &
                               /currentPatch%tlai_profile(cl,ft,iv)
                      end if
                      
                   enddo
                   
                enddo
             enddo
             
             ! --------------------------------------------------------------------------
             ! Set the mask that identifies which PFT x can-layer combinations have
             ! scattering elements in them.
             ! --------------------------------------------------------------------------

             do cl = 1,currentPatch%NCL_p
                do ft = 1,numpft
                   do  iv = 1, currentPatch%nrad(cl,ft)
                      if(currentPatch%canopy_area_profile(cl,ft,iv) > 0._r8)then
                         currentPatch%canopy_mask(cl,ft) = 1     
                      endif
                   end do !iv
                enddo !ft
             enddo ! loop over cl
             
          endif !leaf distribution
          
       end if
       
       currentPatch => currentPatch%younger 
       
    enddo !patch       
    
    return
 end subroutine leaf_area_profile

 ! ======================================================================================

  subroutine update_hlm_dynamics(nsites,sites,fcolumn,bc_out)

     ! ----------------------------------------------------------------------------------
     ! The purpose of this routine is to package output boundary conditions related
     ! to vegetation coverage to the host land model.
     ! ----------------------------------------------------------------------------------

     use EDTypesMod        , only : ed_patch_type, ed_cohort_type, &
                                    ed_site_type, AREA
     use FatesInterfaceMod , only : bc_out_type
     use EDPftvarcon       , only : EDPftvarcon_inst


     !
     ! !ARGUMENTS    
     integer,            intent(in)            :: nsites
     type(ed_site_type), intent(inout), target :: sites(nsites)
     integer,            intent(in)            :: fcolumn(nsites)
     type(bc_out_type),  intent(inout)         :: bc_out(nsites)

     ! Locals
     type (ed_cohort_type) , pointer :: currentCohort
     integer :: s, ifp, c, p
     type (ed_patch_type)  , pointer :: currentPatch
     real(r8) :: bare_frac_area
     real(r8) :: total_patch_area
     real(r8) :: weight  ! Weighting for cohort variables in patch


     do s = 1,nsites

        ifp = 0
        total_patch_area = 0._r8 
        currentPatch => sites(s)%oldest_patch
        c = fcolumn(s)
        do while(associated(currentPatch))
           ifp = ifp+1

           if ( currentPatch%total_canopy_area-currentPatch%area > 0.000001_r8 ) then
              write(fates_log(),*) 'ED: canopy area bigger than area',currentPatch%total_canopy_area ,currentPatch%area
              currentPatch%total_canopy_area = currentPatch%area
           endif


           if (associated(currentPatch%tallest)) then
              bc_out(s)%htop_pa(ifp) = currentPatch%tallest%hite
           else
              ! FIX(RF,040113) - should this be a parameter for the minimum possible vegetation height?
              bc_out(s)%htop_pa(ifp) = 0.1_r8
           endif
           
           bc_out(s)%hbot_pa(ifp) = max(0._r8, min(0.2_r8, bc_out(s)%htop_pa(ifp)- 1.0_r8))

           ! Use leaf area weighting for all cohorts in the patch to define the characteristic
           ! leaf width used by the HLM
           ! ----------------------------------------------------------------------------
!           bc_out(s)%dleaf_pa(ifp) = 0.0_r8
!           if(currentPatch%lai>1.0e-9_r8) then
!              currentCohort => currentPatch%shortest
!              do while(associated(currentCohort))
!                 weight = min(1.0_r8,currentCohort%lai/currentPatch%lai)
!                 bc_out(s)%dleaf_pa(ifp) = bc_out(s)%dleaf_pa(ifp) + &
!                       EDPftvarcon_inst%dleaf(currentCohort%pft)*weight
!                 currentCohort => currentCohort%taller  
!              enddo
!           end if

           ! Roughness length and displacement height are not PFT properties, they are
           ! properties of the canopy assemblage.  Defining this needs an appropriate model.
           ! Right now z0 and d are pft level parameters.  For the time being we will just
           ! use the 1st index until a suitable model is defined. (RGK 04-2017)
           ! -----------------------------------------------------------------------------
           bc_out(s)%z0m_pa(ifp)    = EDPftvarcon_inst%z0mr(1) * bc_out(s)%htop_pa(ifp)
           bc_out(s)%displa_pa(ifp) = EDPftvarcon_inst%displar(1) * bc_out(s)%htop_pa(ifp)
           bc_out(s)%dleaf_pa(ifp)  = EDPftvarcon_inst%dleaf(1)


           ! We are assuming here that grass is all located underneath tree canopies. 
           ! The alternative is to assume it is all spatial distinct from tree canopies.
           ! In which case, the bare area would have to be reduced by the grass area...
           ! currentPatch%total_canopy_area/currentPatch%area is fraction of this patch cover by plants 
           ! currentPatch%area/AREA is the fraction of the soil covered by this patch. 
           
           bc_out(s)%canopy_fraction_pa(ifp) = min(1.0_r8,currentPatch%total_canopy_area/currentPatch%area) * &
                 (currentPatch%area/AREA)

           bare_frac_area = (1.0_r8-min(1.0_r8,currentPatch%total_canopy_area/currentPatch%area))* &
                 (currentPatch%area/AREA)
           
           total_patch_area = total_patch_area + bc_out(s)%canopy_fraction_pa(ifp) + bare_frac_area
    
           ! Calculate area indices for output boundary to HLM
           ! It is assumed that cpatch%canopy_area_profile and cpat%xai_profiles
           ! have been updated (ie ed_leaf_area_profile has been called since dynamics has been called)

           bc_out(s)%elai_pa(ifp) = calc_areaindex(currentPatch,'elai')
           bc_out(s)%tlai_pa(ifp) = calc_areaindex(currentPatch,'tlai')
           bc_out(s)%esai_pa(ifp) = calc_areaindex(currentPatch,'esai')
           bc_out(s)%tsai_pa(ifp) = calc_areaindex(currentPatch,'tsai')

           ! Fraction of vegetation free of snow. This is used to flag those
           ! patches which shall under-go photosynthesis
           ! INTERF-TODO: we may want to stop using frac_veg_nosno_alb and let
           ! FATES internal variables decide if photosynthesis is possible
           ! we are essentially calculating it inside FATES to tell the 
           ! host to tell itself when to do things (circuitous). Just have
           ! to determine where else it is used

           if ((bc_out(s)%elai_pa(ifp) + bc_out(s)%esai_pa(ifp)) > 0._r8) then
              bc_out(s)%frac_veg_nosno_alb_pa(ifp) = 1.0_r8
           else
              bc_out(s)%frac_veg_nosno_alb_pa(ifp) = 0.0_r8
           end if
           
           currentPatch => currentPatch%younger
        end do
        
        if(abs(total_patch_area-1.0_r8)>1e-9)then
           write(fates_log(),*) 'total area is wrong in update_hlm_dynamics',total_patch_area
        endif
        

     end do


  end subroutine update_hlm_dynamics

  ! =====================================================================================

  function calc_areaindex(cpatch,ai_type) result(ai)

     ! ----------------------------------------------------------------------------------
     ! This subroutine calculates the exposed leaf area index of a patch
     ! this is the square meters of leaf per square meter of ground area
     ! It does so by integrating over the depth and functional type profile of leaf area
     ! which are per area of crown.  This value has to be scaled by crown area to convert
     ! to ground area.
     ! ----------------------------------------------------------------------------------
     
     ! Arguments
     type(ed_patch_type),intent(in), target :: cpatch
     character(len=*),intent(in)            :: ai_type

     integer :: cl,ft
     real(r8) :: ai
     ! TODO: THIS MIN LAI IS AN ARTIFACT FROM TESTING LONG-AGO AND SHOULD BE REMOVED
     ! THIS HAS BEEN KEPT THUS FAR TO MAINTAIN B4B IN TESTING OTHER COMMITS
     real(r8),parameter :: ai_min = 0.1_r8
     real(r8),pointer   :: ai_profile

     ai = 0._r8
     if     (trim(ai_type) == 'elai') then
        do cl = 1,cpatch%NCL_p
           do ft = 1,numpft
              ai = ai + sum(cpatch%canopy_area_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                    cpatch%elai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
           enddo
        enddo
     elseif (trim(ai_type) == 'tlai') then
        do cl = 1,cpatch%NCL_p
           do ft = 1,numpft
              ai = ai + sum(cpatch%canopy_area_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                    cpatch%tlai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
           enddo
        enddo
     elseif (trim(ai_type) == 'esai') then
         do cl = 1,cpatch%NCL_p
           do ft = 1,numpft
              ai = ai + sum(cpatch%canopy_area_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                    cpatch%esai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
           enddo
        enddo
     elseif (trim(ai_type) == 'tsai') then
        do cl = 1,cpatch%NCL_p
           do ft = 1,numpft
              ai = ai + sum(cpatch%canopy_area_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                    cpatch%tsai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
           enddo
        enddo
     else

        write(fates_log(),*) 'Unsupported area index sent to calc_areaindex'
        call endrun(msg=errMsg(sourcefile, __LINE__))
     end if
     
     ai = max(ai_min,ai)
          
     return

  end function calc_areaindex

  ! ===============================================================================================
  
  subroutine CanopyLayerArea(currentPatch,site_spread,layer_index,layer_area)
     
     ! --------------------------------------------------------------------------------------------
     ! This function calculates the total crown area footprint for a desired layer of the canopy
     ! within a patch.
     ! The return units are the same as patch%area, which is m2
     ! ---------------------------------------------------------------------------------------------
                   
     ! Arguments
     type(ed_patch_type),intent(inout), target   :: currentPatch
     real(r8),intent(in)                         :: site_spread
     integer,intent(in)                          :: layer_index
     real(r8),intent(inout)                      :: layer_area

     type(ed_cohort_type), pointer :: currentCohort
     
     
     layer_area = 0.0_r8
     currentCohort => currentPatch%tallest
     do while (associated(currentCohort))
        call carea_allom(currentCohort%dbh,currentCohort%n,site_spread, &
              currentCohort%pft,currentCohort%c_area)
        if (currentCohort%canopy_layer .eq. layer_index) then
           layer_area = layer_area + currentCohort%c_area
        end if
        currentCohort => currentCohort%shorter
     enddo
     return
  end subroutine CanopyLayerArea

  ! ===============================================================================================
  
  function NumPotentialCanopyLayers(currentPatch,site_spread,include_substory) result(z)

     ! --------------------------------------------------------------------------------------------
     ! Calculate the number of canopy layers in this patch.
     ! This simple call only determines total layering by querying the cohorts
     ! which layer they are in, it doesn't do any size evaluation.
     ! It may also, optionally, account for the temporary "substory", which is the imaginary
     ! layer below the understory which will be needed to temporarily accomodate demotions from
     ! the understory in the event the understory has reached maximum allowable area.
     ! --------------------------------------------------------------------------------------------

     type(ed_patch_type),target   :: currentPatch
     real(r8),intent(in)          :: site_spread
     logical                      :: include_substory

     type(ed_cohort_type),pointer :: currentCohort
     
     integer :: z
     real(r8) :: c_area
     real(r8) :: arealayer

     z = 1
     currentCohort => currentPatch%tallest
     do while (associated(currentCohort))  
        z = max(z,currentCohort%canopy_layer)
        currentCohort => currentCohort%shorter
     enddo

     if(include_substory)then
        arealayer = 0.0
        currentCohort => currentPatch%tallest
        do while (associated(currentCohort))  
           if(currentCohort%canopy_layer == z) then
              call carea_allom(currentCohort%dbh,currentCohort%n,site_spread,currentCohort%pft,c_area)
              arealayer = arealayer + c_area
           end if
           currentCohort => currentCohort%shorter
        enddo
        
        ! Does the bottom layer have more than a full canopy? 
        ! If so we need to make another layer.
        if(arealayer > currentPatch%area)then
           z = z + 1
        endif
     end if
     
  end function NumPotentialCanopyLayers

end module EDCanopyStructureMod
