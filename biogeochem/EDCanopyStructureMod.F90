module EDCanopyStructureMod

  ! ============================================================================
  ! Code to determine whether the canopy is closed, and which plants are either in the understorey or overstorey
  ! This is obviosuly far too complicated for it's own good and needs re-writing.  
  ! ============================================================================

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesGlobals          , only : fates_log
  use EDPftvarcon             , only : EDPftvarcon_inst
  use EDGrowthFunctionsMod  , only : c_area
  use EDCohortDynamicsMod   , only : copy_cohort, terminate_cohorts, fuse_cohorts
  use EDtypesMod            , only : ed_site_type, ed_patch_type, ed_cohort_type, ncwd
  use EDTypesMod            , only : nclmax
  use EDTypesMod            , only : nlevcan
  use EDTypesMod            , only : numpft_ed
  use FatesGlobals          , only : endrun => fates_endrun

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

  ! 10/30/09: Created by Rosie Fisher
  ! ============================================================================

contains

  ! ============================================================================
  subroutine canopy_structure( currentSite )
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

    use EDParamsMod, only : ED_val_comp_excln, ED_val_ag_biomass
    use SFParamsMod, only : SF_val_cwd_frac
    use EDtypesMod , only : ncwd, min_patch_area
    !
    ! !ARGUMENTS    
    type(ed_site_type) , intent(inout), target   :: currentSite
    !
    ! !LOCAL VARIABLES:
    type(ed_patch_type) , pointer :: currentPatch
    type(ed_cohort_type), pointer :: currentCohort,copyc
    integer  :: i,j
    integer  :: c     ! index for CWD
    integer  :: z     ! Current number of canopy layers. (1= canopy, 2 = understorey) 
    real(r8) :: checkarea
    real(r8) :: cc_loss
    real(r8) :: lossarea
    real(r8) :: newarea
    real(r8) :: arealayer(nlevcan) ! Amount of plant area currently in each canopy layer
    real(r8) :: sumdiff(nlevcan)   ! The total of the exclusion weights for all cohorts in layer z 
    real(r8) :: weight                ! The amount of the total lost area that comes from this cohort
    real(r8) :: sum_weights(nlevcan)
    real(r8) :: new_total_area_check
    real(r8) :: missing_area, promarea,cc_gain,sumgain
    integer  :: promswitch,lower_cohort_switch
    real(r8) :: sumloss,excess_area
    integer  :: count_mi
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
    new_total_area_check = 0._r8
    do while (associated(currentPatch)) ! Patch loop    

       if (currentPatch%area .gt. min_patch_area) then  ! avoid numerical weirdness that shouldn't be happening anyway

       excess_area = 1.0_r8   
        
       ! Does any layer have excess area in it? Keep going until it does not...
       
       do while(excess_area > 0.000001_r8)

          ! Calculate the area currently in each canopy layer. 
          z = 1 
          arealayer = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))  
             currentCohort%c_area = c_area(currentCohort) ! Reassess cohort area. 
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
             z = max(z,currentCohort%canopy_layer) ! What is the current number of canopy layers? 
             currentCohort => currentCohort%shorter
          enddo
          
          ! Does the bottom layer have more than a full canopy? If so we need to make another layer.
          
          if(arealayer(z) > currentPatch%area)then  ! Do we have too much area in either layer?  
              !write(fates_log(),*) 'CANOPY CLOSURE', z
             z = z + 1
          endif

          currentPatch%NCL_p = min(nclmax,z)   ! Set current canopy layer occupancy indicator.  

          do i = 1,z ! Loop around the currently occupied canopy layers. 
             
             do while((arealayer(i)-currentPatch%area) > 0.000001_r8) 
             ! Is this layer currently over-occupied? 
             ! In that case, we need to work out which cohorts to demote. 
             
                sumloss = 0.0_r8
                new_total_area_check = 0.0_r8
                sumdiff(i) = 0.0_r8    
                currentCohort => currentPatch%tallest 
                do while (associated(currentCohort))
                   currentCohort%c_area = c_area(currentCohort)
                   if(arealayer(i) > currentPatch%area.and.currentCohort%canopy_layer == i)then
                      currentCohort%excl_weight = 1.0_r8/(currentCohort%dbh**ED_val_comp_excln)  
                      sumdiff(i) = sumdiff(i) + currentCohort%excl_weight
                   endif
                   currentCohort => currentCohort%shorter  
                enddo !currentCohort

                lossarea = arealayer(i) - currentPatch%area  !how much do we have to lose?
                sum_weights(i) = 0.0_r8
                currentCohort => currentPatch%tallest    !start from the tallest cohort
                 
                ! Correct the demoted cohorts for  
                do while (associated(currentCohort))
                   if(currentCohort%canopy_layer  ==  i) then
                      weight = currentCohort%excl_weight/sumdiff(i)     
                      currentCohort%excl_weight = min(currentCohort%c_area/lossarea, weight)
                      sum_weights(i) = sum_weights(i) + currentCohort%excl_weight
                   endif
                   currentCohort => currentCohort%shorter      
                enddo

                currentCohort => currentPatch%tallest
                do while (associated(currentCohort))      
                   if(currentCohort%canopy_layer == i)then !All the trees in this layer need to lose some area...
                      weight = currentCohort%excl_weight/sum_weights(i)
                      cc_loss = lossarea*weight !what this cohort has to lose. 
                      !-----------Split and copy boundary cohort-----------------!
                      if(cc_loss < currentCohort%c_area)then
                         allocate(copyc)

                         call copy_cohort(currentCohort, copyc) !makes an identical copy...
                         ! n.b this needs to happen BEFORE the cohort goes into the new layer, 
                         ! otherwise currentPatch%spread(i+1) will be higher and the area will change...!!! 
                         sumloss = sumloss + cc_loss 

                         newarea = currentCohort%c_area - cc_loss
                         copyc%n = currentCohort%n*newarea/currentCohort%c_area   !
                         currentCohort%n = currentCohort%n - (currentCohort%n*newarea/currentCohort%c_area) !     

                         copyc%canopy_layer = i !the taller cohort is the copy
                         currentCohort%canopy_layer = i + 1 !demote the current cohort to the understory.           
                         ! seperate cohorts. 
                         ! - 0.000000000001_r8 !needs to be a very small number to avoid 
                         ! causing non-linearity issues with c_area.  is this really required? 
                         currentCohort%dbh = currentCohort%dbh 
                         copyc%dbh = copyc%dbh !+ 0.000000000001_r8

                         ! keep track of number and biomass of demoted cohort
                         currentSite%demotion_rate(currentCohort%size_class) = &
                              currentSite%demotion_rate(currentCohort%size_class) + currentCohort%n
                         currentSite%demotion_carbonflux = currentSite%demotion_carbonflux + &
                              currentCohort%b * currentCohort%n

                         !kill the ones which go into canopy layers that are not allowed... (default nclmax=2) 
                         if(i+1 > nclmax)then 
                           !put the litter from the terminated cohorts into the fragmenting pools
                          ! write(fates_log(),*) '3rd canopy layer'
                            do c=1,ncwd

                               currentPatch%CWD_AG(c)  = currentPatch%CWD_AG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    ED_val_ag_biomass * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area  
         
                               currentPatch%CWD_BG(c)  = currentPatch%CWD_BG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    (1.0_r8-ED_val_ag_biomass) * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area !litter flux per m2.

                            enddo

                                currentPatch%leaf_litter(currentCohort%pft)  = &
                                     currentPatch%leaf_litter(currentCohort%pft) + (currentCohort%bl)* &
                                          currentCohort%n/currentPatch%area ! leaf litter flux per m2.

                                currentPatch%root_litter(currentCohort%pft)  = &
                                     currentPatch%root_litter(currentCohort%pft) + &
                                     (currentCohort%br+currentCohort%bstore)*currentCohort%n/currentPatch%area
   									
                            currentCohort%n = 0.0_r8
                            currentCohort%c_area = 0._r8
                         else  
                            currentCohort%c_area = c_area(currentCohort)       
                         endif
                         copyc%c_area = c_area(copyc)
                         new_total_area_check = new_total_area_check+copyc%c_area

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
                      else
                         currentCohort%canopy_layer = i + 1 !the whole cohort becomes demoted
                         sumloss = sumloss + currentCohort%c_area 

                         ! keep track of number and biomass of demoted cohort
                         currentSite%demotion_rate(currentCohort%size_class) = &
                              currentSite%demotion_rate(currentCohort%size_class) + currentCohort%n
                         currentSite%demotion_carbonflux = currentSite%demotion_carbonflux + &
                              currentCohort%b * currentCohort%n

                         !kill the ones which go into canopy layers that are not allowed... (default nclmax=2) 
                         if(i+1 > nclmax)then  

                           !put the litter from the terminated cohorts into the fragmenting pools
                            do c=1,ncwd

                               currentPatch%CWD_AG(c)  = currentPatch%CWD_AG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    ED_val_ag_biomass * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area           
                               currentPatch%CWD_BG(c)  = currentPatch%CWD_BG(c) + (currentCohort%bdead+currentCohort%bsw) * &
                                    (1.0_r8-ED_val_ag_biomass) * &
                                    SF_val_CWD_frac(c)*currentCohort%n/currentPatch%area !litter flux per m2.

                            enddo

                                currentPatch%leaf_litter(currentCohort%pft)  = &
                                     currentPatch%leaf_litter(currentCohort%pft) + currentCohort%bl* &
                                          currentCohort%n/currentPatch%area ! leaf litter flux per m2.

                                currentPatch%root_litter(currentCohort%pft)  = &
                                     currentPatch%root_litter(currentCohort%pft) + &
                                     (currentCohort%br+currentCohort%bstore)*currentCohort%n/currentPatch%area
                            currentCohort%n = 0.0_r8
                            currentCohort%c_area = 0._r8

                         else  
                            currentCohort%c_area = c_area(currentCohort)       
                         endif

                         !write(fates_log(),*) 'demoting whole cohort', currentCohort%c_area,cc_loss, &
                              !currentCohort%canopy_layer,currentCohort%dbh

                      endif

                      !----------- End of cohort splitting ------------------------------!             
                   endif !canopy layer = i

                   currentCohort => currentCohort%shorter

                enddo !currentCohort 
                
                arealayer(i) = arealayer(i) - sumloss
                !Update arealayer for diff calculations of layer below. 
                arealayer(i + 1) = arealayer(i + 1) + sumloss 

             enddo !arealayer loop
             if(arealayer(i)-currentPatch%area > 0.00001_r8)then
                write(fates_log(),*) 'lossarea problem', lossarea,sumloss,z,currentPatch%patchno
             endif

          enddo !z  

          z = 1
          arealayer = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))  
             currentCohort%c_area = c_area(currentCohort)
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
             z = max(z,currentCohort%canopy_layer)
             currentCohort => currentCohort%shorter
          enddo

          !does the bottom layer have more than a full canopy? If so we need to make another layer.
          if(arealayer(z) > currentPatch%area)then
             z = z + 1
          endif
          excess_area = 0.0_r8
          do j=1,z
             if(arealayer(j) > currentPatch%area)then 
                excess_area = arealayer(j)-currentPatch%area
             endif
          enddo
          currentPatch%ncl_p = min(z,nclmax)

       enddo !is there still excess area in any layer?      

       call fuse_cohorts(currentPatch)
       call terminate_cohorts(currentSite, currentPatch)

       ! ----------- Check cohort area ------------------------------!
       do i = 1,z
          checkarea = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             if(currentCohort%canopy_layer == i)then
                checkarea = checkarea + c_area(currentCohort)
             endif

             currentCohort=>currentCohort%shorter

          enddo

       enddo ! 


       ! ----------- Check whether the intended 'full' layers are actually filling all the space.
       ! If not, promote some fraction of cohorts upwards ------------------------------!     
       ! THIS SECTION MIGHT BE TRIGGERED BY A FIRE OR MORTALITY EVENT, FOLLOWED BY A PATCH FUSION, 
       ! SO THE TOP LAYER IS NO LONGER FULL...  

       promswitch = 0 

       missing_area=1.0_r8    
       count_mi = 0
       !does any layer have excess area in it? keep going until it does not...
       do while(missing_area > 0.000001_r8.and.z > 1) 
          count_mi = count_mi +1
          do i = 1,z-1 ! if z is greater than one, there is a possibility of too many plants in the understorey. 
             lower_cohort_switch = 1
             ! is the area of the layer less than the area of the patch, if it is supposed to be closed (z>1) 
             do while((arealayer(i)-currentPatch%area) < -0.000001_r8.and.lower_cohort_switch == 1) 

                if(arealayer(i+1) <= 0.000001_r8)then
                   currentCohort => currentPatch%tallest 
                   arealayer = 0._r8
                   do while (associated(currentCohort))            
                      if(currentCohort%canopy_layer == i+1)then !look at the cohorts in the canopy layer below... 
                         currentCohort%canopy_layer = i   
                         currentCohort%c_area = c_area(currentCohort)

                         ! keep track of number and biomass of promoted cohort
                         currentSite%promotion_rate(currentCohort%size_class) = &
                              currentSite%promotion_rate(currentCohort%size_class) + currentCohort%n
                         currentSite%promotion_carbonflux = currentSite%promotion_carbonflux + &
                              currentCohort%b * currentCohort%n

                        ! write(fates_log(),*) 'promoting very small cohort', currentCohort%c_area,currentCohort%canopy_layer
                      endif
                      arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer)+currentCohort%c_area 
                      currentCohort => currentCohort%shorter   
                   enddo

                endif !promoting all of the small amount of area in the lower layers. 


                lower_cohort_switch = 0 
                sumgain = 0.0_r8
                sumdiff(i) = 0.0_r8    
                ! figure out with what weighting we need to promote cohorts.
                ! This is the opposite of the demotion weighting... 
                currentCohort => currentPatch%tallest 
                do while (associated(currentCohort))
                   currentCohort%c_area = c_area(currentCohort)
                   if(currentCohort%canopy_layer == i+1)then !look at the cohorts in the canopy layer below... 
                      currentCohort%prom_weight = currentCohort%dbh**ED_val_comp_excln   !as opposed to 1/(dbh^C_e) 
                      sumdiff(i) = sumdiff(i) + currentCohort%prom_weight
                   endif
                   currentCohort => currentCohort%shorter  
                enddo !currentCohort

                promarea =  currentPatch%area -arealayer(i) !how much do we need to gain?
                sum_weights(i) = 0.0_r8
                currentCohort => currentPatch%tallest    !start from the tallest cohort

                do while (associated(currentCohort))
                   if(currentCohort%canopy_layer  ==  i+1) then !still looking at the layer beneath. 
                      weight = currentCohort%prom_weight/sumdiff(i)
                      if(promarea > 0._r8)then    
                         currentCohort%prom_weight = min(currentCohort%c_area/promarea, weight)
                      else
                         currentCohort%prom_weight = 0._r8
                      endif
                      sum_weights(i) = sum_weights(i) + currentCohort%prom_weight
                   endif
                   currentCohort => currentCohort%shorter      
                enddo

                currentCohort => currentPatch%tallest
                do while (associated(currentCohort))      
                   if(currentCohort%canopy_layer == i+1)then !All the trees in this layer need to promote some area upwards... 
                      lower_cohort_switch = 1
                      weight = currentCohort%prom_weight/sum_weights(i)
                      cc_gain = promarea*weight !what this cohort has to promote. 
                      !-----------Split and copy boundary cohort-----------------!
                      if(cc_gain < currentCohort%c_area)then
                         allocate(copyc)

                         call copy_cohort(currentCohort, copyc) !makes an identical copy...
                         ! n.b this needs to happen BEFORE the cohort goes into the new layer, otherwise currentPatch
                         ! %spread(+1) will be higher and the area will change...!!!
                         sumgain = sumgain + cc_gain


                         newarea = currentCohort%c_area - cc_gain !new area of existing cohort
                         copyc%n = currentCohort%n*cc_gain/currentCohort%c_area   !number of individuals in promoted cohort. 
                         ! number of individuals in cohort remaining in understorey    
                         currentCohort%n = currentCohort%n - (currentCohort%n*cc_gain/currentCohort%c_area) 

                         currentCohort%canopy_layer = i+1  !keep current cohort in the understory.        
                         copyc%canopy_layer = i ! promote copy to the higher canopy layer. 

                         ! keep track of number and biomass of promoted cohort
                         currentSite%promotion_rate(copyc%size_class) = &
                              currentSite%promotion_rate(copyc%size_class) + copyc%n
                         currentSite%promotion_carbonflux = currentSite%promotion_carbonflux + &
                              copyc%b * copyc%n

                         ! seperate cohorts. 
                         ! needs to be a very small number to avoid causing non-linearity issues with c_area. 
                         ! is this really required? 
                         currentCohort%dbh = currentCohort%dbh - 0.000000000001_r8 
                         copyc%dbh = copyc%dbh + 0.000000000001_r8

                         currentCohort%c_area = c_area(currentCohort)          
                         copyc%c_area = c_area(copyc)

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
                      else
                         currentCohort%canopy_layer = i  !the whole cohort becomes promoted
                         sumgain = sumgain + currentCohort%c_area  !inserting deliberate mistake to see how far we make it... 
                         ! update area AFTER we sum up the losses. the cohort may shrink at this point,
                         ! if the upper canopy spread is smaller. this shold be dealt with by the 'excess area' loop.  
                         currentCohort%c_area = c_area(currentCohort) 

                         ! keep track of number and biomass of promoted cohort
                         currentSite%promotion_rate(currentCohort%size_class) = &
                              currentSite%promotion_rate(currentCohort%size_class) + currentCohort%n
                         currentSite%promotion_carbonflux = currentSite%promotion_carbonflux + &
                              currentCohort%b * currentCohort%n

                         promswitch = 1

                        ! write(fates_log(),*) 'promoting whole cohort', currentCohort%c_area,cc_gain,currentCohort%canopy_layer, &
                              !currentCohort%pft,currentPatch%patchno

                      endif
                      if(promswitch == 1)then
                        ! write(fates_log(),*) 'cohort loop',currentCohort%pft,currentPatch%patchno
                      endif
                      !----------- End of cohort splitting ------------------------------!             
                   else
                      if(promswitch == 1)then
                         ! write(fates_log(),*) 'cohort list',currentCohort%pft, &
                             ! currentCohort%canopy_layer,currentCohort%c_area
                      endif
                   endif

                   currentCohort => currentCohort%shorter
                enddo !currentCohort 
                arealayer(i) = arealayer(i) + sumgain
                arealayer(i + 1) = arealayer(i + 1) - sumgain !Update arealayer for diff calculations of layer below. 

                if(promswitch == 1)then
                  ! write(fates_log(),*) 'arealayer loop',arealayer(1:3),currentPatch%area,promarea,sumgain, &
                        !currentPatch%patchno,z,i,lower_cohort_switch
                endif
                if(promswitch == 1.and.associated(currentPatch%tallest))then
                   ! write(fates_log(),*) 'cohorts',currentCohort%pft,currentPatch%patchno, &
                        !currentCohort%c_area
                endif
             enddo !arealayer loop

             if(currentPatch%area-arealayer(i) < 0.000001_r8)then
                !write(fates_log(),*) 'gainarea problem',sumgain,arealayer(i),currentPatch%area,z, &
                     !currentPatch%patchno,currentPatch%area - arealayer(i),i,missing_area,count_mi
             endif
             if(promswitch == 1)then
               ! write(fates_log(),*) 'z loop',arealayer(1:3),currentPatch%patchno,z
             endif
          enddo !z  

          z = 1
          arealayer = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))  
             currentCohort%c_area = c_area(currentCohort)
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
             z = max(z,currentCohort%canopy_layer)
             currentCohort => currentCohort%shorter
          enddo

          missing_area = 0.0_r8
          do j=1,z-1
             if(arealayer(j) < currentPatch%area)then  !this is the amount of area that we still have spare in this layer. 
                missing_area = currentPatch%area - arealayer(j)
                if(missing_area <= 0.000001_r8.and.missing_area > 0._r8)then
                   missing_area = 0.0_r8
                  ! write(fates_log(),*) 'correcting MI',j,currentPatch%area - arealayer(j)
                endif
             endif
          enddo
          currentPatch%ncl_p = min(z,nclmax)
          if(promswitch == 1)then
            ! write(fates_log(),*) 'missingarea loop',arealayer(1:3),currentPatch%patchno,missing_area,z
          endif
       enddo !is there still not enough canopy area in any layer?         

       call fuse_cohorts(currentPatch)
       call terminate_cohorts(currentSite, currentPatch)

       if(promswitch == 1)then
          !write(fates_log(),*) 'going into cohort check'
       endif
       ! ----------- Check cohort area ------------------------------!
       do i = 1,z
          checkarea = 0.0_r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             if(currentCohort%canopy_layer == i)then
                checkarea = checkarea + c_area(currentCohort)
             endif

             currentCohort => currentCohort%shorter

          enddo

          if(((checkarea-currentPatch%area)) > 0.0001)then
             write(fates_log(),*) 'problem with canopy area', checkarea, currentPatch%area, checkarea - currentPatch%area, &
                  i, z, missing_area 
             currentCohort => currentPatch%tallest
             do while (associated(currentCohort))
             if(currentCohort%canopy_layer == i)then
                write(fates_log(),*) 'c_areas in top layer', c_area(currentCohort)
             endif
             currentCohort => currentCohort%shorter

          enddo
                  
          endif

          if ( i  >  1) then
             if ( (arealayer(i)  -  arealayer(i-1) )>1e-11 ) then
                write(fates_log(),*) 'smaller top layer than bottom layer ',arealayer(i),arealayer(i-1), &
                     currentPatch%area,currentPatch%spread(i-1:i)
             endif
          endif
       enddo ! 

       if(promswitch == 1)then 
         ! write(fates_log(),*) 'end patch loop',currentSite%clmgcell
       endif

       else !terminate  logic to only do if patch_area_sufficiently large
          write(fates_log(),*) 'canopy_structure: patch area too small.', currentPatch%area
       end if 
       

       currentPatch => currentPatch%younger
    enddo !patch  

    if(promswitch == 1)then
      ! write(fates_log(),*) 'end  canopy structure',currentSite%clmgcell
    endif

  end subroutine canopy_structure

  ! ============================================================================
  subroutine canopy_spread( currentSite )
    !
    ! !DESCRIPTION:
    !  Calculates the spatial spread of tree canopies based on canopy closure.                             
    !
    ! !USES:
    use EDParamsMod , only : ED_val_maxspread, ED_val_minspread 
    !
    ! !ARGUMENTS    
    type (ed_site_type), intent(inout), target :: currentSite
    !
    ! !LOCAL VARIABLES:
    type (ed_cohort_type), pointer :: currentCohort
    type (ed_patch_type) , pointer :: currentPatch
    real(r8) :: arealayer(nlevcan) ! Amount of canopy in each layer. 
    real(r8) :: inc                   ! Arbitrary daily incremental change in canopy area 
    integer  :: z
    !----------------------------------------------------------------------

    inc = 0.005_r8

    currentPatch => currentSite%oldest_patch

    do while (associated(currentPatch))

       !calculate canopy area in each canopy storey...
       arealayer = 0.0_r8   
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))
          currentCohort%c_area = c_area(currentCohort) 
          if(EDPftvarcon_inst%woody(currentCohort%pft) == 1)then
             arealayer(currentCohort%canopy_layer) = arealayer(currentCohort%canopy_layer) + currentCohort%c_area
          endif
          currentCohort => currentCohort%shorter
       enddo

       !If the canopy area is approaching closure, squash the tree canopies and make them taller and thinner
       do z = 1,nclmax  
         
          if(arealayer(z)/currentPatch%area > 0.9_r8)then
             currentPatch%spread(z) = currentPatch%spread(z) - inc
          else 
             currentPatch%spread(z) = currentPatch%spread(z) + inc 
          endif
          if(currentPatch%spread(z) >= ED_val_maxspread)then 
             currentPatch%spread(z) = ED_val_maxspread
          endif
          if(currentPatch%spread(z) <=  ED_val_minspread)then
             currentPatch%spread(z) = ED_val_minspread
          endif
        enddo !z
        !write(fates_log(),*) 'spread',currentPatch%spread(1:2)
        !currentPatch%spread(:) = ED_val_maxspread
        !FIX(RF,033114) spread is off
        !write(fates_log(),*) 'canopy_spread',currentPatch%area,currentPatch%spread(1:2)
        currentPatch => currentPatch%younger

    enddo !currentPatch

  end subroutine canopy_spread


  ! =====================================================================================

  subroutine canopy_summarization( nsites, sites, bc_in )

     ! ----------------------------------------------------------------------------------
     ! Much of this routine was once ed_clm_link minus all the IO and history stuff
     ! ---------------------------------------------------------------------------------

    use FatesInterfaceMod    , only : bc_in_type
    use EDPatchDynamicsMod   , only : set_patchno
    use EDPatchDYnamicsMod   , only : set_root_fraction
    use EDCohortDynamicsMod  , only : size_and_type_class_index
    use EDGrowthFunctionsMod , only : tree_lai, c_area
    use EDEcophysConType     , only : EDecophyscon
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
          
          call set_root_fraction(currentPatch,bc_in(s)%depth_gl)

          !zero cohort-summed variables. 
          currentPatch%total_canopy_area = 0.0_r8
          currentPatch%total_tree_area = 0.0_r8
          currentPatch%lai = 0.0_r8
          canopy_leaf_area = 0.0_r8
          
          !update cohort quantitie s                                  
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))
             
             ft = currentCohort%pft

             
             ! Update the cohort's index within the size bin classes
             ! Update the cohort's index within the SCPF classification system
             call size_and_type_class_index(currentCohort%dbh,currentCohort%pft, &
                                            currentCohort%size_class,currentCohort%size_by_pft_class)

             
             currentCohort%b = currentCohort%balive+currentCohort%bdead+currentCohort%bstore
             currentCohort%treelai = tree_lai(currentCohort)
             
             currentCohort%c_area = c_area(currentCohort)
             canopy_leaf_area = canopy_leaf_area + currentCohort%treelai *currentCohort%c_area
                  
             if(currentCohort%canopy_layer==1)then
                currentPatch%total_canopy_area = currentPatch%total_canopy_area + currentCohort%c_area
                if(EDPftvarcon_inst%woody(ft)==1)then
                   currentPatch%total_tree_area = currentPatch%total_tree_area + currentCohort%c_area
                endif
             endif
             
             ! Check for erroneous zero values. 
             if(currentCohort%dbh <= 0._r8 .or. currentCohort%n == 0._r8)then
                write(fates_log(),*) 'ED: dbh or n is zero in canopy_summarization', currentCohort%dbh,currentCohort%n
             endif
             if(currentCohort%pft == 0.or.currentCohort%canopy_trim <= 0._r8)then
                write(fates_log(),*) 'ED: PFT or trim is zero in canopy_summarization',currentCohort%pft,currentCohort%canopy_trim
             endif
             if(currentCohort%balive <= 0._r8)then
                write(fates_log(),*) 'ED: balive is zero in canopy_summarization',currentCohort%balive
             endif

             currentCohort => currentCohort%taller
             
          enddo ! ends 'do while(associated(currentCohort))
          
          if ( currentPatch%total_canopy_area-currentPatch%area > 0.000001_r8 ) then
             write(fates_log(),*) 'ED: canopy area bigger than area',currentPatch%total_canopy_area ,currentPatch%area
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
    !
    ! !DESCRIPTION:
    !
    ! !USES:

    use EDGrowthFunctionsMod , only : tree_lai, tree_sai, c_area 
    use EDtypesMod           , only : area, dinc_ed, hitemax, n_hite_bins
    use EDEcophysConType     , only : EDecophyscon
  
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
    integer  :: L                        ! Canopy layer index
    integer  :: p                        ! clm patch index  
    real(r8) :: fraction_exposed         ! how much of this layer is not covered by snow?
    real(r8) :: layer_top_hite           ! notional top height of this canopy layer (m)
    real(r8) :: layer_bottom_hite        ! notional bottom height of this canopy layer (m)
    integer  :: smooth_leaf_distribution ! is the leaf distribution this option (1) or not (0)
    real(r8) :: frac_canopy(N_HITE_BINS) ! amount of canopy in each height class
    real(r8) :: minh(N_HITE_BINS)        ! minimum height in height class (m)
    real(r8) :: maxh(N_HITE_BINS)        ! maximum height in height class (m)
    real(r8) :: dh                       ! vertical detph of height class (m)
    real(r8) :: min_chite                ! bottom of cohort canopy  (m)
    real(r8) :: max_chite                ! top of cohort canopy      (m)
    real(r8) :: lai                      ! summed lai for checking m2 m-2
    real(r8) :: snow_depth_avg           ! avg snow over whole site
    integer  :: NC                       ! number of cohorts, for bug fixing. 
    
    !----------------------------------------------------------------------

    smooth_leaf_distribution = 0

    ! Here we are trying to generate a profile of leaf area, indexed by 'z' and by pft
    ! We assume that each point in the canopy recieved the light attenuated by the average
    ! leaf area index above it, irrespective of PFT identity... 
    ! Each leaf is defined by how deep in the canopy it is, in terms of LAI units.  (FIX(RF,032414), GB)
    
    currentPatch => currentSite%oldest_patch   ! ed patch
    do while(associated(currentPatch))
       
       !Calculate tree and canopy areas. 
       currentPatch%canopy_area = 0._r8
       currentPatch%canopy_layer_lai(:) = 0._r8
       NC = 0
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort))       
          currentCohort%c_area = c_area(currentCohort)
          currentPatch%canopy_area = currentPatch%canopy_area + currentCohort%c_area
          NC = NC+1
          currentCohort => currentCohort%taller    
       enddo
         ! if plants take up all the tile, then so does the canopy.  
       currentPatch%canopy_area = min(currentPatch%canopy_area,currentPatch%area) 
       
       !calculate tree lai and sai.
       currentPatch%ncan(:,:) = 0 
       currentPatch%nrad(:,:) = 0 
       currentPatch%lai = 0._r8
       currentCohort => currentPatch%shortest
       do while(associated(currentCohort)) 
          currentCohort%treelai = tree_lai(currentCohort)    
          currentCohort%treesai = tree_sai(currentCohort)
          currentCohort%lai =  currentCohort%treelai *currentCohort%c_area/currentPatch%canopy_area 
          currentCohort%sai =  currentCohort%treesai *currentCohort%c_area/currentPatch%canopy_area  
          !Calculate the LAI plus SAI in each canopy storey. 
          currentCohort%NV =  ceiling((currentCohort%treelai+currentCohort%treesai)/dinc_ed)  
          
          currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft) = &
                max(currentPatch%ncan(currentCohort%canopy_layer,currentCohort%pft),currentCohort%NV)
          currentPatch%lai = currentPatch%lai +currentCohort%lai
          
          do L = 1,nclmax-1
             if(currentCohort%canopy_layer == L)then
                currentPatch%canopy_layer_lai(L) = currentPatch%canopy_layer_lai(L) + currentCohort%lai + &
                      currentCohort%sai
             endif
          enddo
          
          currentCohort => currentCohort%taller 
          
       enddo !currentCohort
       currentPatch%nrad = currentPatch%ncan

       if(smooth_leaf_distribution == 1)then
          ! we are going to ignore the concept of canopy layers, and put all of the leaf area into height banded bins. 
          ! using the same domains as we had before, except that CL always = 1
          currentPatch%tlai_profile = 0._r8
          currentPatch%tsai_profile = 0._r8  
          currentPatch%elai_profile = 0._r8
          currentPatch%esai_profile = 0._r8  
          
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
             min_chite = currentCohort%hite - currentCohort%hite * EDecophyscon%crown(ft)
             max_chite = currentCohort%hite  
             do iv = 1,N_HITE_BINS  
                frac_canopy(iv) = 0.0_r8
                ! this layer is in the middle of the canopy
                if(max_chite > maxh(iv).and.min_chite < minh(iv))then 
                   frac_canopy(iv)= min(1.0_r8,dh / (currentCohort%hite*EDecophyscon%crown(ft)))
                   ! this is the layer with the bottom of the canopy in it. 
                elseif(min_chite < maxh(iv).and.min_chite > minh(iv).and.max_chite > maxh(iv))then 
                   frac_canopy(iv) = (maxh(iv) -min_chite ) / (currentCohort%hite*EDecophyscon%crown(ft))
                   ! this is the layer with the top of the canopy in it. 
                elseif(max_chite > minh(iv).and.max_chite < maxh(iv).and.min_chite < minh(iv))then 
                   frac_canopy(iv) = (max_chite - minh(iv)) / (currentCohort%hite*EDecophyscon%crown(ft))
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
          
          !check
          currentPatch%lai = 0._r8
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort)) 
             currentPatch%lai = currentPatch%lai +currentCohort%lai
             currentCohort => currentCohort%taller   
          enddo !currentCohort
          lai = 0.0_r8
          do ft = 1,numpft_ed
             lai = lai+ sum(currentPatch%tlai_profile(1,ft,:))
          enddo
          
          if(lai > currentPatch%lai)then
             write(fates_log(), *) 'ED: problem with lai assignments'
          endif
          
          
       else ! smooth leaf distribution  
          !Go through all cohorts and add their leaf area and canopy area to the accumulators. 
          currentPatch%tlai_profile = 0._r8
          currentPatch%tsai_profile = 0._r8  
          currentPatch%elai_profile = 0._r8
          currentPatch%esai_profile = 0._r8 
          currentPatch%layer_height_profile = 0._r8
          currentPatch%canopy_area_profile(:,:,:) = 0._r8       
          currentPatch%ncan(:,:) = 0 
          currentPatch%nrad(:,:) = 0 
          currentCohort => currentPatch%shortest
          
          do while(associated(currentCohort))   
             L = currentCohort%canopy_layer
             ft = currentCohort%pft 
             !Calculate the number of layers of thickness dlai, including the last one. 
             currentCohort%NV =  ceiling((currentCohort%treelai+currentCohort%treesai)/dinc_ed)
             !how much of each tree is stem area index? Assuming that there is 
             if(currentCohort%treelai+currentCohort%treesai > 0._r8)then    
                fleaf = currentCohort%lai / (currentCohort%lai + currentCohort%sai) 
             else
                fleaf = 0._r8
                write(fates_log(), *) 'ED: no stem or leaf area' ,currentCohort%pft,currentCohort%bl, &
                      currentCohort%balive,currentCohort%treelai,currentCohort%treesai,currentCohort%dbh, &
                      currentCohort%n,currentCohort%status_coh
             endif
             currentPatch%ncan(L,ft) = max(currentPatch%ncan(L,ft),currentCohort%NV)  
             currentPatch%nrad(L,ft) = currentPatch%ncan(L,ft)  !fudge - this needs to be altered for snow burial
             if(currentCohort%NV > currentPatch%nrad(L,ft))then
                write(fates_log(), *) 'ED: issue with NV',currentCohort%NV,currentCohort%pft,currentCohort%canopy_layer
             endif
             
             !Whole layers.  Make a weighted average of the leaf area in each layer before dividing it by the total area. 
             !fill up layer for whole layers.  FIX(RF,032414)- for debugging jan 2012
             do iv = 1,currentCohort%NV-1 
                
                ! what is the height of this layer? (for snow burial purposes...)  
                ! EDPftvarcon_inst%vertical_canopy_frac(ft))! fudge - this should be pft specific but i cant get it to compile. 
                layer_top_hite = currentCohort%hite-((iv/currentCohort%NV) * currentCohort%hite * &
                      EDecophyscon%crown(currentCohort%pft) )
                layer_bottom_hite = currentCohort%hite-(((iv+1)/currentCohort%NV) * currentCohort%hite * &
                      EDecophyscon%crown(currentCohort%pft)) ! EDPftvarcon_inst%vertical_canopy_frac(ft))
                
                fraction_exposed =1.0_r8
                
                currentPatch%tlai_profile(L,ft,iv) = currentPatch%tlai_profile(L,ft,iv)+ dinc_ed * fleaf * &
                      currentCohort%c_area/currentPatch%total_canopy_area
                currentPatch%elai_profile(L,ft,iv) = currentPatch%elai_profile(L,ft,iv)+ dinc_ed * fleaf * &
                      currentCohort%c_area/currentPatch%total_canopy_area * fraction_exposed
                
                currentPatch%tsai_profile(L,ft,iv) = currentPatch%tsai_profile(L,ft,iv)+ dinc_ed * (1._r8 - fleaf) * &
                      currentCohort%c_area/currentPatch%total_canopy_area
                currentPatch%esai_profile(L,ft,iv) = currentPatch%esai_profile(L,ft,iv)+ dinc_ed * (1._r8 - fleaf) * &
                      currentCohort%c_area/currentPatch%total_canopy_area * fraction_exposed
                
                currentPatch%canopy_area_profile(L,ft,iv) =  min(1.0_r8,currentPatch%canopy_area_profile(L,ft,iv) + &
                      currentCohort%c_area/currentPatch%total_canopy_area)
                currentPatch%layer_height_profile(L,ft,iv) = currentPatch%layer_height_profile(L,ft,iv) + (dinc_ed * fleaf * &
                      currentCohort%c_area/currentPatch%total_canopy_area *(layer_top_hite+layer_bottom_hite)/2.0_r8) !average height of layer. 
                
                if ( DEBUG ) then
                   write(fates_log(), *) 'calc snow 2', snow_depth_si , frac_sno_eff_si
                   write(fates_log(), *) 'LHP', currentPatch%layer_height_profile(L,ft,iv)
                   write(fates_log(), *) 'EDCLMLink 1246 ', currentPatch%elai_profile(1,ft,iv)
                end if

             end do
                  
             !Bottom layer
             iv = currentCohort%NV
             ! EDPftvarcon_inst%vertical_canopy_frac(ft))! fudge - this should be pft specific but i cant get it to compile.
             layer_top_hite = currentCohort%hite-((iv/currentCohort%NV) * currentCohort%hite * &
                   EDecophyscon%crown(currentCohort%pft) )
             ! EDPftvarcon_inst%vertical_canopy_frac(ft))
             layer_bottom_hite = currentCohort%hite-(((iv+1)/currentCohort%NV) * currentCohort%hite * &
                   EDecophyscon%crown(currentCohort%pft))
             
             fraction_exposed = 1.0_r8 !default. 
             snow_depth_avg = snow_depth_si * frac_sno_eff_si
             if(snow_depth_avg  > layer_top_hite)then
                fraction_exposed = 0._r8
             endif
             if(snow_depth_avg < layer_bottom_hite)then
                fraction_exposed = 1._r8
                
             endif
             if(snow_depth_avg>= layer_bottom_hite.and.snow_depth_avg <= layer_top_hite)then !only partly hidden...                                   
                fraction_exposed =  max(0._r8,(min(1.0_r8,(snow_depth_avg-layer_bottom_hite)/ &
                      (layer_top_hite-layer_bottom_hite ))))
             endif
             fraction_exposed= 1.0_r8

             
             remainder = (currentCohort%treelai + currentCohort%treesai) - (dinc_ed*(currentCohort%NV-1))
             if(remainder > 1.0_r8)then
                write(fates_log(), *)'ED: issue with remainder',currentCohort%treelai,currentCohort%treesai,dinc_ed, & 
                      currentCohort%NV
             endif
             !assumes that fleaf is unchanging FIX(RF,032414)
             
             currentPatch%tlai_profile(L,ft,iv) =  currentPatch%tlai_profile(L,ft,iv)+ remainder * fleaf * &
                   currentCohort%c_area/currentPatch%total_canopy_area
             currentPatch%elai_profile(L,ft,iv) = currentPatch%elai_profile(L,ft,iv) + remainder * fleaf * &
                   currentCohort%c_area/currentPatch%total_canopy_area * fraction_exposed
             !assumes that fleaf is unchanging FIX(RF,032414)
             
             currentPatch%tsai_profile(L,ft,iv) =  currentPatch%tsai_profile(L,ft,iv)+  remainder * &
                   (1.0_r8-fleaf) * currentCohort%c_area/currentPatch%total_canopy_area
             currentPatch%esai_profile(L,ft,iv) = currentPatch%esai_profile(L,ft,iv)+  remainder * &
                   (1.0_r8-fleaf) * currentCohort%c_area/currentPatch%total_canopy_area * fraction_exposed
             
             currentPatch%canopy_area_profile(L,ft,iv) = min(1.0_r8,currentPatch%canopy_area_profile(L,ft,iv) + &
                   currentCohort%c_area/currentPatch%total_canopy_area)
             currentPatch%layer_height_profile(L,ft,iv) = currentPatch%layer_height_profile(L,ft,iv) + (remainder * fleaf * &
                   currentCohort%c_area/currentPatch%total_canopy_area*(layer_top_hite+layer_bottom_hite)/2.0_r8)
             if ( DEBUG ) write(fates_log(), *) 'LHP', currentPatch%layer_height_profile(L,ft,iv)
             if(currentCohort%dbh <= 0._r8.or.currentCohort%n == 0._r8)then
                write(fates_log(), *) 'ED: dbh or n is zero in clmedlink', currentCohort%dbh,currentCohort%n
             endif
             if(currentCohort%pft == 0.or.currentCohort%canopy_trim <= 0._r8)then
                write(fates_log(), *) 'ED: PFT or trim is zero in clmedlink',currentCohort%pft,currentCohort%canopy_trim
             endif
             if(currentCohort%balive <= 0._r8.or.currentCohort%bl < 0._r8)then
                write(fates_log(), *) 'ED: balive is zero in clmedlink',currentCohort%balive,currentCohort%bl
             endif
             
             currentCohort => currentCohort%taller
             
          enddo !cohort
          
          do L = 1,currentPatch%NCL_p
             do ft = 1,numpft_ed
                do iv = 1,currentPatch%nrad(L,ft)
                   !account for total canopy area
                   currentPatch%tlai_profile(L,ft,iv) = currentPatch%tlai_profile(L,ft,iv) / &
                         currentPatch%canopy_area_profile(L,ft,iv)
                   currentPatch%tsai_profile(L,ft,iv) = currentPatch%tsai_profile(L,ft,iv) / &
                         currentPatch%canopy_area_profile(L,ft,iv)
                   
                   if ( DEBUG ) write(fates_log(), *) 'EDCLMLink 1293 ', currentPatch%elai_profile(L,ft,iv)
                   
                   currentPatch%elai_profile(L,ft,iv) = currentPatch%elai_profile(L,ft,iv) / &
                         currentPatch%canopy_area_profile(L,ft,iv)
                   currentPatch%esai_profile(L,ft,iv) = currentPatch%esai_profile(L,ft,iv) / &
                         currentPatch%canopy_area_profile(L,ft,iv)
                   currentPatch%layer_height_profile(L,ft,iv) = currentPatch%layer_height_profile(L,ft,iv) &
                         /currentPatch%tlai_profile(L,ft,iv)
                enddo
                
                currentPatch%tlai_profile(L,ft,currentPatch%nrad(L,ft)+1: nlevcan) = 0._r8
                currentPatch%tsai_profile(L,ft,currentPatch%nrad(L,ft)+1: nlevcan) = 0._r8
                currentPatch%elai_profile(L,ft,currentPatch%nrad(L,ft)+1: nlevcan) = 0._r8 
                currentPatch%esai_profile(L,ft,currentPatch%nrad(L,ft)+1: nlevcan) = 0._r8
                
             enddo
          enddo
          
          currentPatch%nrad = currentPatch%ncan
          do L = 1,currentPatch%NCL_p
             do ft = 1,numpft_ed
                if(currentPatch%nrad(L,ft) > 30)then
                   write(fates_log(), *) 'ED: issue w/ nrad'
                endif
                currentPatch%present(L,ft) = 0
                do  iv = 1, currentPatch%nrad(L,ft);
                   if(currentPatch%canopy_area_profile(L,ft,iv) > 0._r8)then
                      currentPatch%present(L,ft) = 1     
                   endif
                end do !iv
             enddo !ft
             
             if ( L == 1 .and. abs(sum(currentPatch%canopy_area_profile(1,1:numpft_ed,1))) < 0.99999  &
                   .and. currentPatch%NCL_p > 1 ) then
                write(fates_log(), *) 'ED: canopy area too small',sum(currentPatch%canopy_area_profile(1,1:numpft_ed,1))
                write(fates_log(), *) 'ED: cohort areas', currentPatch%canopy_area_profile(1,1:numpft_ed,:)
             endif
             
             if (L == 1 .and. currentPatch%NCL_p > 1 .and.  &
                   abs(sum(currentPatch%canopy_area_profile(1,1:numpft_ed,1))) < 0.99999) then
                write(fates_log(), *) 'ED: not enough area in the top canopy', &
                      sum(currentPatch%canopy_area_profile(L,1:numpft_ed,1)), &
                      currentPatch%canopy_area_profile(L,1:numpft_ed,1)
             endif
             
             if(abs(sum(currentPatch%canopy_area_profile(L,1:numpft_ed,1))) > 1.00001)then
                write(fates_log(), *) 'ED: canopy-area-profile wrong', &
                      sum(currentPatch%canopy_area_profile(L,1:numpft_ed,1)), &
                      currentPatch%patchno, L
                write(fates_log(), *) 'ED: areas',currentPatch%canopy_area_profile(L,1:2,1),currentPatch%patchno
                
                currentCohort => currentPatch%shortest
                
                do while(associated(currentCohort))
                   
                   if(currentCohort%canopy_layer==1)then
                      write(fates_log(), *) 'ED: cohorts',currentCohort%dbh,currentCohort%c_area, &
                            currentPatch%total_canopy_area,currentPatch%area,currentPatch%canopy_area
                      write(fates_log(), *) 'ED: fracarea', currentCohort%pft, &
                            currentCohort%c_area/currentPatch%total_canopy_area
                   endif
                   
                   currentCohort => currentCohort%taller  
                   
                enddo !currentCohort
             endif
          enddo ! loop over L
          
          do L = 1,currentPatch%NCL_p
             do ft = 1,numpft_ed
                if(currentPatch%present(L,FT) > 1)then
                   write(fates_log(), *) 'ED: present issue',L,ft,currentPatch%present(L,FT)
                   currentPatch%present(L,ft) = 1
                endif
             enddo
          enddo
          
       endif !leaf distribution
       
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

     !
     ! !ARGUMENTS    
     integer,            intent(in)            :: nsites
     type(ed_site_type), intent(inout), target :: sites(nsites)
     integer,            intent(in)            :: fcolumn(nsites)
     type(bc_out_type),  intent(inout)         :: bc_out(nsites)

     ! Locals
     integer :: s, ifp, c
     type (ed_patch_type)  , pointer :: currentPatch
     real(r8) :: bare_frac_area
     real(r8) :: total_patch_area

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
           do ft = 1,numpft_ed
              ai = ai + sum(cpatch%canopy_area_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                    cpatch%elai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
           enddo
        enddo
     elseif (trim(ai_type) == 'tlai') then
        do cl = 1,cpatch%NCL_p
           do ft = 1,numpft_ed
              ai = ai + sum(cpatch%canopy_area_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                    cpatch%tlai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
           enddo
        enddo
     elseif (trim(ai_type) == 'esai') then
         do cl = 1,cpatch%NCL_p
           do ft = 1,numpft_ed
              ai = ai + sum(cpatch%canopy_area_profile(cl,ft,1:cpatch%nrad(cl,ft)) * &
                    cpatch%esai_profile(cl,ft,1:cpatch%nrad(cl,ft)))
           enddo
        enddo
     elseif (trim(ai_type) == 'tsai') then
        do cl = 1,cpatch%NCL_p
           do ft = 1,numpft_ed
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

end module EDCanopyStructureMod
