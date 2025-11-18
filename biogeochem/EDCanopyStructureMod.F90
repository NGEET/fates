module EDCanopyStructureMod

  ! =====================================================================================
  ! Code to determine whether the canopy is closed, and which plants are either in the
  ! understorey or overstorey. This is obviosuly far too complicated for it's own good
  ! =====================================================================================

  use FatesConstantsMod     , only : r8 => fates_r8
  use FatesConstantsMod     , only : itrue, ifalse
  use FatesConstantsMod     , only : tinyr8
  use FatesConstantsMod     , only : nearzero, area_error_1
  use FatesConstantsMod     , only : rsnbl_math_prec
  use FatesConstantsMod     , only : nocomp_bareground
  use FatesConstantsMod,      only : i_term_mort_type_canlev
  use FatesGlobals          , only : fates_log
  use EDPftvarcon           , only : EDPftvarcon_inst
  use PRTParametersMod      , only : prt_params
  use FatesAllometryMod     , only : carea_allom
  use EDCohortDynamicsMod   , only : terminate_cohorts, terminate_cohort, fuse_cohorts
  use EDCohortDynamicsMod   , only : InitPRTObject
  use FatesAllometryMod     , only : tree_lai_sai
  use EDTypesMod            , only : ed_site_type
  use EDTypesMod            , only : set_patchno
  use FatesAllometryMod     , only : VegAreaLayer
  use FatesAllometryMod     , only : CrownDepth
  use FatesPatchMod,          only : fates_patch_type
  use FatesCohortMod,         only : fates_cohort_type
  use EDParamsMod           , only : nclmax
  use EDParamsMod           , only : nlevleaf
  use EDParamsMod           , only : GetNVegLayers
  use EDParamsMod           , only : comp_excln_exp
  use EDtypesMod            , only : AREA
  use EDLoggingMortalityMod , only : UpdateHarvestC
  use FatesGlobals          , only : endrun => fates_endrun
  use FatesInterfaceTypesMod     , only : hlm_days_per_year
  use FatesInterfaceTypesMod     , only : hlm_use_planthydro
  use FatesInterfaceTypesMod     , only : hlm_use_cohort_age_tracking
  use FatesInterfaceTypesMod     , only : hlm_use_sp
  use FatesInterfaceTypesMod     , only : numpft
  use FatesInterfaceTypesMod, only : bc_in_type
  use FatesPlantHydraulicsMod, only : UpdateH2OVeg,InitHydrCohort, RecruitWaterStorage
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : leaf_organ
  use PRTGenericMod,          only : fnrt_organ
  use PRTGenericMod,          only : sapw_organ
  use PRTGenericMod,          only : store_organ
  use PRTGenericMod,          only : repro_organ
  use PRTGenericMod,          only : struct_organ
  use PRTGenericMod,          only : SetState
  use PRTGenericMod,          only : carbon12_element
  use FatesTwoStreamUtilsMod, only : FatesConstructRadElements
  use FatesRadiationMemMod  , only : twostr_solver
  use FatesRadiationMemMod  , only : num_rad_stream_types
  
  ! CIME Globals
  use shr_log_mod           , only : errMsg => shr_log_errMsg

  implicit none
  private

  public :: canopy_structure
  public :: canopy_spread
  public :: calc_areaindex
  public :: canopy_summarization
  public :: update_hlm_dynamics
  public :: UpdateFatesAvgSnowDepth
  public :: UpdatePatchLAI
  public :: UpdateCohortLAI

  logical, parameter :: debug=.false.

  character(len=*), parameter, private :: sourcefile = &
       __FILE__

  integer :: istat           ! return status code
  character(len=255) :: smsg ! Message string for deallocation errors

  ! Precision targets for demotion and promotion
  ! We have two:
  ! "pa_area_target_precision" is the required precision at the patch level,
  !    we keep shuffling and splitting cohorts until each layer is within this precision
  ! "co_area_target_precision" is the required precision at the cohort level,
  !    essentially it is the minimum amount of change required to not ignore
  !    a partial promotion or demotion
  
  real(r8), parameter :: pa_area_target_precision = 1.0E-11_r8
  real(r8), parameter :: co_area_target_precision = 1.0E-12_r8 

  integer, parameter :: demotion_phase  = 1
  integer, parameter :: promotion_phase = 2
  
  ! will attempt to reduce errors
  ! below this level

  real(r8), parameter :: area_check_precision  = 1.0E-7_r8     ! Area conservation checks must
  ! be within this absolute tolerance
  real(r8), parameter :: area_check_rel_precision = 1.0E-4_r8  ! Area conservation checks must
  ! be within this relative tolerance

  real(r8), parameter :: similar_height_tol = 1.0E-3_r8    ! I think trees that differ by 1mm
  ! can be roughly considered the same right?

  logical, parameter :: preserve_b4b = .true.


  ! If we want to allow some degree of imperfection
  ! in canopy closure we would add it here
  real(r8), parameter :: imperfect_fraction = 0._r8

  
  ! 10/30/09: Created by Rosie Fisher
  ! 2017/2018/2025: Modifications and updates by Ryan Knox
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
    ! parameter (comp_excln_exp).

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

    use EDTypesMod , only : min_patch_area

    !
    ! !ARGUMENTS
    type(ed_site_type) , intent(inout), target   :: currentSite
    type(bc_in_type), intent(in)                 :: bc_in

    !
    ! !LOCAL VARIABLES:
    type(fates_patch_type) , pointer :: currentPatch
    type(fates_cohort_type), pointer :: currentCohort
    integer  :: i_lyr                  ! current layer index
    integer  :: z                      ! Current number of canopy layers. (1= canopy, 2 = understorey)
    integer  :: ipft
    real(r8) :: arealayer(nclmax+5)    ! Amount of plant area currently in each canopy layer
    integer  :: patch_area_counter     ! count iterations used to solve canopy areas
    logical  :: area_not_balanced      ! logical controlling if the patch layer areas
    real(r8) :: target_area            ! Canopy area that is either in excess/defiency
                                       ! that is slated for demotion/promotion from/into layer
    
    ! have successfully been redistributed
    integer  :: return_code            ! math checks on variables will return>0 if problems exist
    ! We only iterate because of possible imprecisions generated by the cohort
    ! termination process.  These should be super small, so at the most
    ! try to re-balance 3 times.  If that doesn't give layer areas
    ! within tolerance of canopy area, there is something wrong

    integer, parameter  :: max_patch_iterations = 10


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

       ! Make sure we are sorted
       if(debug) call currentPatch%SortCohorts(check_order=.true.)

       ! Terminate cohorts before organizing canopy. That
       ! step will be interested in preserving area, so termination
       ! during that step will be counter productive
       call terminate_cohorts(currentSite, currentPatch, -1,13,bc_in)
       call terminate_cohorts(currentSite, currentPatch, 1,13,bc_in)
       call terminate_cohorts(currentSite, currentPatch, 2,13,bc_in)       
       
       ! ------------------------------------------------------------------------------
       ! Perform numerical checks on some cohort and patch structures
       ! ------------------------------------------------------------------------------

       ! canopy layer has a special bounds check
       currentCohort => currentPatch%tallest
       do while (associated(currentCohort))
          currentCohort%canopy_layer_yesterday = currentCohort%canopy_layer
          if( currentCohort%canopy_layer < 1 ) then
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
          z = NumCanopyLayers(currentPatch)

          do i_lyr = 1,z ! Loop around the currently occupied canopy layers.
             call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer(i_lyr))
             target_area = max(0._r8,arealayer(i_lyr) - (1._r8-imperfect_fraction)*currentPatch%area)
             call PromoteOrDemote(currentSite, currentPatch, i_lyr, demotion_phase, target_area)
          end do

          ! Terminate only for type 1 (near zero number density)
          call terminate_cohorts(currentSite, currentPatch,1,23,bc_in)
          call fuse_cohorts(currentSite, currentPatch, bc_in)

          ! ---------------------------------------------------------------------------------------
          ! Promotion Phase: Identify if any upper-layers are underful and layers below them
          ! have cohorts that can be split and promoted to the layer above.
          ! ---------------------------------------------------------------------------------------

          ! Re-calculate Number of layers
          z = NumCanopyLayers(currentPatch)

          ! We only promote if we have at least two layers
          if (z>1) then
             do i_lyr=2,z
                call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr-1,arealayer(i_lyr-1))
                target_area = max(0._r8,(1._r8-imperfect_fraction)*currentPatch%area - arealayer(i_lyr-1))
                call PromoteOrDemote(currentSite, currentPatch, i_lyr, promotion_phase, target_area)
             end do

             ! Terminate only for type 1 (near zero number density)
             call terminate_cohorts(currentSite, currentPatch,1,24,bc_in)
             call fuse_cohorts(currentSite, currentPatch, bc_in)

          end if

          ! ---------------------------------------------------------------------------------------
          ! Check on Layer Area (if the layer differences are not small
          ! Continue trying to demote/promote. Its possible on the first pass through,
          ! that cohort fusion has nudged the areas a little bit.
          ! On all but the bottom layer, we expect the areas to match the area of the
          ! patch with small precision, since we assume a PPA. On the lowest layer,
          ! we only expect the area to be below the patch area.
          ! ---------------------------------------------------------------------------------------

          z = NumCanopyLayers(currentPatch)
          area_not_balanced = .false.
          do i_lyr = 1,z
             call CanopyLayerArea(currentPatch,currentSite%spread,i_lyr,arealayer(i_lyr))
             if(i_lyr < z)then
                if (abs(arealayer(i_lyr)-(1._r8-imperfect_fraction)*currentPatch%area) > area_check_precision) then
                   area_not_balanced = .true.
                end if
             else
                if ((arealayer(i_lyr)-(1._r8-imperfect_fraction)*currentPatch%area) > area_check_precision) then
                   area_not_balanced = .true.
                end if
             end if
          enddo
          
          ! ---------------------------------------------------------------------------------------
          ! Gracefully exit if too many iterations have gone by
          ! ---------------------------------------------------------------------------------------

          patch_area_counter = patch_area_counter + 1
          if(patch_area_counter > max_patch_iterations .and. area_not_balanced) then
             write(fates_log(),*) 'PATCH AREA CHECK NOT CLOSING'
             write(fates_log(),*) 'patch area:',currentpatch%area
             write(fates_log(),*) 'fraction that is imperfect (unclosed):',imperfect_fraction
             write(fates_log(),*) 'lat:',currentSite%lat
             write(fates_log(),*) 'lon:',currentSite%lon
             write(fates_log(),*) 'spread:',currentSite%spread
             do i_lyr = 1,z
                write(fates_log(),*) '-----------------------------------------'
                write(fates_log(),*) 'layer: ',i_lyr,' area: ',arealayer(i_lyr)
                write(fates_log(),*) 'bias [m2] (layer-patch): ',(arealayer(i_lyr)- &
                     (1._r8-imperfect_fraction)*currentPatch%area)
                currentCohort => currentPatch%tallest
                do while (associated(currentCohort))
                   if(currentCohort%canopy_layer == i_lyr)then
                      write(fates_log(),*) '-----------'
                      write(fates_log(),*) ' co area:',currentCohort%c_area
                      write(fates_log(),*) ' co dbh: ',currentCohort%dbh
                      write(fates_log(),*) ' co pft: ',currentCohort%pft
                      write(fates_log(),*) ' co n: ',currentCohort%n
                   end if
                   currentCohort => currentCohort%shorter
                end do
             enddo

             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

       enddo ! do while(area_not_balanced)


       ! Terminate any cohorts that are still outside the maximum number of
       ! canopy layers.  These terminations only occur in level 3
       call terminate_cohorts(currentSite, currentPatch, 3,17,bc_in)

       z = NumCanopyLayers(currentPatch)
       
       ! Save number of canopy layers to the patch structure
       if(z > nclmax) then
          write(fates_log(),*) 'Termination should have ensured number'
          write(fates_log(),*) 'of canopy layers was not larger than nclmax'
          write(fates_log(),*) 'Predicted: ',z
          write(fates_log(),*) 'nclmax: ',nclmax
          write(fates_log(),*) 'Consider increasing nclmax if this value is to low'
          write(fates_log(),*) 'and you think this number of canopy layers is reasonable.'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       else
          currentPatch%NCL_p = z
       end if
       
       ! -------------------------------------------------------------------------------------------
       ! if we are using "strict PPA", then calculate a z_star value as
       ! the height of the smallest tree in the canopy
       ! loop from top to bottom and locate the shortest cohort in level 1 whose shorter
       ! neighbor is in level 2 set zstar as the ehight of that shortest level 1 cohort
       ! -------------------------------------------------------------------------------------------

       if ( comp_excln_exp .lt. 0.0_r8) then
          currentPatch%zstar = 0._r8
          currentCohort => currentPatch%tallest
          do while (associated(currentCohort))
             if(currentCohort%canopy_layer .eq. 2)then
                if (associated(currentCohort%taller)) then
                   if (currentCohort%taller%canopy_layer .eq. 1 ) then
                      currentPatch%zstar = currentCohort%taller%height
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

  subroutine PromoteOrDemote(site,patch,target_layer,phase,target_area)

    ! --------------------------------------------------------------
    ! This routine will:
    ! 1) Identify the list of cohorts that are in the appropriate
    !    layer for promotion or demotion into the adjacent
    !    layer
    ! 2) Calculate the combined crown area of those cohorts
    !    that will be transferred to the adjacent layer
    ! 3) Perform the transfer either by re-assignment (if whole)
    !    of by splitting the cohort
    ! 4) Track the abundance and mass flows when promoting/demoting
    ! --------------------------------------------------------------

    ! Arguments
    type(ed_site_type)      :: site
    type(fates_patch_type)  :: patch
    integer,intent(in)      :: target_layer ! Canopy layer we draw from
    integer,intent(in)      :: phase        ! promotion or demotion?
    real(r8),intent(in)     :: target_area  ! Area we want to move [m2/ha]

    ! Locals
    type(fates_cohort_type), pointer :: cohort
    type(fates_cohort_type), pointer :: copyc

    real(r8) :: promdem_area    ! Actual area promoted or demoted (minimum of target
                                ! and existing canopy area)
    real(r8) :: sumpd_area      ! Sum crown area of all cohorts in layer [m2/ha]
    real(r8) :: group_area      ! Sum area of cohorts with the same height [m2/ha]
    real(r8) :: remainder_area  ! The area that has not been accounted
    real(r8) :: excess_area     ! The area that could not be accounted
    real(r8) :: attempt_area    ! Amount of area attempted to donate probabilistically
    real(r8) :: max_donate_area ! This is the total area of the layer
                                ! for as seeks to fill out the target_area
    real(r8) :: leaf_c, store_c
    real(r8) :: fnrt_c, sapw_c
    real(r8) :: struct_c
    integer  :: ilyr_change     ! layer offset from current for the destination (+/- 1)
    integer  :: ic,ic_n,ic_nn   ! Cohort indices
    integer  :: n_layer         ! The number of cohorts in the layer


    if (target_area<nearzero) return

    ! Use the patch's scratch vector of cohorts
    ! to help track which cohorts are in the target layer
    associate(layer_co => patch%co_scr)

      ! Step 1: Determine which cohorts are in the layer
      !         and point to them in the scratch vector
      !         Make sure their areas are updated too.
      !         We point to them in the scratch
      !         vector in order of promotion/demotion,
      !         note that this is inconsequential for probabalistic

      ic = 0
      group_area = 0._r8
      if(phase==demotion_phase) then
         cohort => patch%shortest
         ilyr_change = 1
      else
         cohort => patch%tallest
         ilyr_change = -1
      end if
      do while (associated(cohort))
         if(cohort%canopy_layer == target_layer)then
            ic = ic + 1
            call carea_allom(cohort%dbh,cohort%n,site%spread, &
                 cohort%pft,cohort%crowndamage,cohort%c_area)
            group_area = group_area + cohort%c_area
            layer_co(ic)%p => cohort
         end if
         if(phase==demotion_phase) then
            cohort => cohort%taller
         else
            cohort => cohort%shorter
         end if
      end do

      ! We update the target area to be no more than the
      ! area of the layer (can't take more than there is..)
      promdem_area = min(target_area,group_area)

      ! Store the number of cohorts in the layer
      ! and zero out the array of area transfers
      n_layer = ic

      do ic = 1,n_layer
         layer_co(ic)%pd_area = 0._r8
      end do

      ! Step 2: Calculate the promotion or demotion areas      
      comp_excl_type: if (comp_excln_exp .ge. 0.0_r8 ) then

         ! ------------------------------------------------------------------
         ! Stochastic case
         ! ------------------------------------------------------------------

         if_not_trivial: if(target_area >= group_area)then

            ! If promotion/demotion is so large that it
            ! is larger than available area in the layer,
            ! the trivial solution is to just set the
            ! promotion/demotion areas to the cohort areas
            ! in the layer and move on. (We don't need to
            ! do this for rank-ordered, that algorithm
            ! handles this just fine and would just make
            ! it less readable)"
            
            do ic = 1,n_layer
               layer_co(ic)%pd_area = layer_co(ic)%p%c_area
            end do

         else

            sumpd_area = 0._r8
            do ic = 1,n_layer
               cohort => layer_co(ic)%p
               if(phase==demotion_phase) then
                  layer_co(ic)%pd_area = cohort%c_area/(cohort%height**comp_excln_exp)
               elseif(phase==promotion_phase) then
                  layer_co(ic)%pd_area = cohort%c_area*cohort%height**comp_excln_exp
               end if
               sumpd_area = sumpd_area + layer_co(ic)%pd_area
            end do

            ! Distribute areas in a first pass
            ! For those cohorts where more area was to be donated
            ! than it has, accumulate the excess. For those
            ! cohorts that are not filled and still have area to
            ! donate, accumulate remainder. We will use these in
            ! the next step to portion out area.

            excess_area = 0._r8
            remainder_area = 0._r8
            do ic = 1,n_layer
               cohort => layer_co(ic)%p
               attempt_area = promdem_area*layer_co(ic)%pd_area/sumpd_area
               if(attempt_area>cohort%c_area)then
                  excess_area  = excess_area + (attempt_area - cohort%c_area)
               else
                  remainder_area = remainder_area + (cohort%c_area - attempt_area)
               end if
               layer_co(ic)%pd_area = min(cohort%c_area,attempt_area)
            end do

            if(excess_area>nearzero)then

               ! The "if_not_trivial" condition above should prevent
               ! a situation at this point in the code where all
               ! promotion/demotions are larger than the layer area
               ! and thus remainder area is zero.
               if(remainder_area<=0._r8)then
                  write(fates_log(),*) 'prob. prom/dem has encountered a situation'
                  write(fates_log(),*) 'where it thinks all weightings are greater'
                  write(fates_log(),*) 'than the layer areas, but somehow'
                  write(fates_log(),*) 'was not captured at the beginning of the'
                  write(fates_log(),*) 'routine in the if_not_trivial clause.'
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               end if
               
               do ic = 1,n_layer
                  cohort => layer_co(ic)%p
                  ! look at just the cohorts that still have space to give
                  ! remove from them the same fraction of their remaining space
                  if (abs(layer_co(ic)%pd_area-cohort%c_area) > nearzero) then
                     layer_co(ic)%pd_area = layer_co(ic)%pd_area + &
                          (excess_area/remainder_area) * &
                          (cohort%c_area - layer_co(ic)%pd_area)
                  end if
               end do
            end if
         end if if_not_trivial
         
      else !comp_excl_exp < 0

         ! ------------------------------------------------------------------
         ! Rank Ordered Case
         ! ------------------------------------------------------------------

         sumpd_area = 0._r8
         ic  = 1
         do while( ic<=n_layer .and. (promdem_area-sumpd_area)>co_area_target_precision) 

            cohort => layer_co(ic)%p

            ! Determine if the next cohorts in
            ! order have the same height

            group_area = cohort%c_area
            ic_n       = ic
            check_next:do while(ic_n<n_layer)
               if( abs(cohort%height-layer_co(ic_n+1)%p%height) > similar_height_tol ) then
                  exit check_next
               else
                  ic_n = ic_n + 1
                  group_area = group_area+layer_co(ic_n)%p%c_area
               end if
            end do check_next

            remainder_area = min(promdem_area-sumpd_area,group_area)
            do ic_nn = ic,ic_n
               layer_co(ic_nn)%pd_area = remainder_area*layer_co(ic_nn)%p%c_area/group_area
               sumpd_area = sumpd_area + layer_co(ic_nn)%pd_area
            end do

            ic = ic_n + 1

         end do
      end if comp_excl_type

      ! Part 3:
      ! Apply the area changes by splitting the cohort and re-assigning
      ! either all or part of it to a new layer
      ! Check to make sure the changes are within bounds

      ic_loop0: do ic = 1,n_layer

         cohort => layer_co(ic)%p

         ! If the dem/prom area is the same area as the
         !    cohort itself, move the whole thing
         ! If the dem/prom area is less than the cohort area
         !    and not trivialy small (larger than precision
         !    check), then split it and move part of it
         ! If the dem/prom area is less than zero or larger than
         !    the cohort area within precision checks then fail
         
         
         whole_or_part: if( ((layer_co(ic)%pd_area - cohort%c_area) > co_area_target_precision ) .or. &
              (layer_co(ic)%pd_area < 0._r8) ) then
            write(fates_log(),*) 'negative,or more area than the cohort has is being promoted/demoted'
            write(fates_log(),*) 'change: ',layer_co(ic)%pd_area
            write(fates_log(),*) 'existing area:',cohort%c_area
            write(fates_log(),*) 'excess: ',layer_co(ic)%pd_area - cohort%c_area
            call endrun(msg=errMsg(sourcefile, __LINE__))

         
         elseif ( abs(layer_co(ic)%pd_area - cohort%c_area) < co_area_target_precision ) then

            ! Whole cohort promotion/demotion
            cohort%canopy_layer = cohort%canopy_layer + ilyr_change
            
         elseif( layer_co(ic)%pd_area > 0._r8 ) then

            ! Partial cohort promotion/demotion
            ! Make a copy of the current cohort.  The copy and the original
            ! conserve total number density.  The copy
            ! remains in the upper-story.  The original is the one
            ! demoted to the understory

            allocate(copyc)

            ! (keep as an example)
            ! Initialize running means
            !allocate(copyc%tveg_lpa)
            !!allocate(copyc%l2fr_ema)
            !  Note, no need to give a starter value here,
            !  that will be taken care of in copy()
            !!call copyc%l2fr_ema%InitRMean(ema_60day)

            ! Initialize the PARTEH object and point to the
            ! correct boundary condition fields
            copyc%prt => null()
            call InitPRTObject(copyc%prt)

            if( hlm_use_planthydro.eq.itrue ) then
               call InitHydrCohort(site,copyc)
            endif

            call cohort%Copy(copyc)
            call copyc%InitPRTBoundaryConditions()

            remainder_area = cohort%c_area - layer_co(ic)%pd_area
            copyc%n = cohort%n*min(1._r8,max(0._r8,remainder_area/cohort%c_area))
            cohort%n = cohort%n - copyc%n

            ! The copied cohort is the part that remains in-layer
            copyc%canopy_layer = cohort%canopy_layer

            ! The original cohort changes layers
            cohort%canopy_layer = cohort%canopy_layer + ilyr_change

            call carea_allom(copyc%dbh,copyc%n,site%spread,copyc%pft, &
                 copyc%crowndamage, copyc%c_area)
            call carea_allom(cohort%dbh,cohort%n,site%spread, &
                 cohort%pft, cohort%crowndamage, cohort%c_area)

            !----------- Insert copy into linked list ------------------------!
            ! Since we are not changing the heights, no sorting necessary
            !-----------------------------------------------------------------!
            copyc%shorter => cohort
            if(associated(cohort%taller))then
               copyc%taller => cohort%taller
               cohort%taller%shorter => copyc
            else
               patch%tallest => copyc
               copyc%taller => null()
            endif
            cohort%taller => copyc

         end if whole_or_part

         ! Part 4:
         ! keep track of number and biomass promoted/demoted

         if( layer_co(ic)%pd_area > 0._r8 ) then
            leaf_c   = cohort%prt%GetState(leaf_organ,carbon12_element)
            store_c  = cohort%prt%GetState(store_organ,carbon12_element)
            fnrt_c   = cohort%prt%GetState(fnrt_organ,carbon12_element)
            sapw_c   = cohort%prt%GetState(sapw_organ,carbon12_element)
            struct_c = cohort%prt%GetState(struct_organ,carbon12_element)
            
            if(phase==demotion_phase) then
               site%demotion_rate(cohort%size_class) = &
                    site%demotion_rate(cohort%size_class) + cohort%n
               site%demotion_carbonflux = site%demotion_carbonflux + &
                    (leaf_c + store_c + fnrt_c + sapw_c + struct_c) * cohort%n
            else
               site%promotion_rate(cohort%size_class) = &
                    site%promotion_rate(cohort%size_class) + cohort%n
               site%promotion_carbonflux = site%promotion_carbonflux + &
                    (leaf_c + store_c + fnrt_c + sapw_c + struct_c) * cohort%n
            end if
         end if
            
      end do ic_loop0

    end associate

  end subroutine PromoteOrDemote

  ! ==============================================================================================
  
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
    type (fates_cohort_type), pointer :: currentCohort
    type (fates_patch_type) , pointer :: currentPatch
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
          call carea_allom(currentCohort%dbh,currentCohort%n, &
               currentSite%spread,currentCohort%pft,currentCohort%crowndamage,currentCohort%c_area)
          if( ( prt_params%woody(currentCohort%pft) .eq. itrue ) .and. &
               (currentCohort%canopy_layer .eq. 1 ) ) then
             sitelevel_canopyarea = sitelevel_canopyarea + currentCohort%c_area
          endif
          currentCohort => currentCohort%shorter
       enddo

       currentPatch => currentPatch%younger

    enddo !currentPatch

    ! If the canopy area is approaching closure,
    ! squash the tree canopies and make them taller and thinner
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

    use FatesInterfaceTypesMod    , only : hlm_use_cohort_age_tracking
    use FatesInterfaceTypesMod    , only : hlm_radiation_model
    use FatesSizeAgeTypeIndicesMod, only : sizetype_class_index
    use FatesSizeAgeTypeIndicesMod, only : coagetype_class_index
    use EDtypesMod                , only : area
    use FatesConstantsMod         , only : itrue

    ! !ARGUMENTS
    integer                 , intent(in)            :: nsites
    type(ed_site_type)      , intent(inout), target :: sites(nsites)
    type(bc_in_type)        , intent(in)            :: bc_in(nsites)
    !
    ! !LOCAL VARIABLES:
    type (fates_patch_type)  , pointer :: currentPatch
    type (fates_cohort_type) , pointer :: currentCohort
    integer  :: s
    integer  :: ft               ! plant functional type
    integer  :: ifp              ! the number of the vegetated patch (1,2,3). In SP mode bareground patch is 0
    integer  :: patchn           ! identification number for each patch.
    real(r8) :: leaf_c           ! leaf carbon [kg]
    real(r8) :: fnrt_c           ! fineroot carbon [kg]
    real(r8) :: sapw_c           ! sapwood carbon [kg]
    real(r8) :: store_c          ! storage carbon [kg]
    real(r8) :: struct_c         ! structure carbon [kg]

    !----------------------------------------------------------------------

    if ( debug ) then
       write(fates_log(),*) 'in canopy_summarization'
    endif

    do s = 1,nsites

       ! --------------------------------------------------------------------------------
       ! Set the patch indices (this is usefull mostly for communicating with a host or
       ! driving model.  Loops through all patches and sets cpatch%patchno to the integer
       ! order of oldest to youngest where the oldest is 1.
       ! --------------------------------------------------------------------------------
       call set_patchno( sites(s) , .false., 0)

       currentPatch => sites(s)%oldest_patch

       do while(associated(currentPatch))

          !zero cohort-summed variables.
          currentPatch%total_canopy_area = 0.0_r8
          currentPatch%total_tree_area = 0.0_r8

          !update cohort quantitie s
          currentCohort => currentPatch%shortest
          do while(associated(currentCohort))

             ft = currentCohort%pft
             leaf_c   = currentCohort%prt%GetState(leaf_organ, carbon12_element)
             sapw_c   = currentCohort%prt%GetState(sapw_organ, carbon12_element)
             struct_c = currentCohort%prt%GetState(struct_organ, carbon12_element)
             fnrt_c   = currentCohort%prt%GetState(fnrt_organ, carbon12_element)
             store_c  = currentCohort%prt%GetState(store_organ, carbon12_element)

             ! Update the cohort's index within the size bin classes
             ! Update the cohort's index within the SCPF classification system
             call sizetype_class_index(currentCohort%dbh,currentCohort%pft, &
                  currentCohort%size_class,currentCohort%size_by_pft_class)

             if (hlm_use_cohort_age_tracking .eq. itrue) then
                call coagetype_class_index(currentCohort%coage,currentCohort%pft, &
                     currentCohort%coage_class,currentCohort%coage_by_pft_class)
             end if

             if(hlm_use_sp.eq.ifalse)then
                call carea_allom(currentCohort%dbh,currentCohort%n,sites(s)%spread,&
                     currentCohort%pft,currentCohort%crowndamage, currentCohort%c_area)
             endif

             if(currentCohort%canopy_layer==1)then
                currentPatch%total_canopy_area = currentPatch%total_canopy_area + currentCohort%c_area
                if( prt_params%woody(ft) == itrue)then
                   currentPatch%total_tree_area = currentPatch%total_tree_area + currentCohort%c_area
                endif
             endif

             ! adding checks for SP and NOCOMP modes.
             if(currentPatch%nocomp_pft_label.eq.nocomp_bareground)then
                write(fates_log(),*) 'cohorts in barepatch',currentPatch%total_canopy_area,currentPatch%nocomp_pft_label
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if

             if(hlm_use_sp.eq.itrue)then

                if(associated(currentPatch%tallest%shorter))then
                   write(fates_log(),*) 'more than one cohort in SP mode',s,currentPatch%nocomp_pft_label
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if

                if (currentPatch%total_canopy_area - (1._r8-imperfect_fraction)*currentPatch%area > area_error_1) then
                   write(fates_log(),*) 'too much canopy in summary', s, &
                        currentPatch%nocomp_pft_label, currentPatch%total_canopy_area - (1._r8-imperfect_fraction)*currentPatch%area
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
             end if  !sp mode

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
             if( (sapw_c + leaf_c + fnrt_c) <= 0._r8)then
                write(fates_log(),*) 'FATES: alive biomass is zero in canopy_summarization', &
                     sapw_c + leaf_c + fnrt_c
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif

             currentCohort => currentCohort%taller

          enddo ! ends 'do while(associated(currentCohort))

          if ( currentPatch%total_canopy_area>currentPatch%area ) then
             if ( currentPatch%total_canopy_area-currentPatch%area > 0.001_r8 ) then
                write(fates_log(),*) 'FATES: canopy area bigger than area', &
                     currentPatch%total_canopy_area ,currentPatch%area, &
                     currentPatch%total_canopy_area -currentPatch%area,&
                     currentPatch%nocomp_pft_label
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
             currentPatch%total_canopy_area = currentPatch%area
          endif

          currentPatch => currentPatch%younger
       end do !patch loop

       call leaf_area_profile(sites(s))
       
       if(hlm_radiation_model.eq.twostr_solver) then
          call FatesConstructRadElements(sites(s))
       end if
       
    end do ! site loop

    return
  end subroutine canopy_summarization

  ! ====================================================================================

  subroutine UpdateFatesAvgSnowDepth(sites,bc_in)

    ! This routine updates the snow depth used in FATES to occlude vegetation
    ! Currently this average takes into account the depth of snow and the
    ! areal coverage fraction

    type(ed_site_type)      , intent(inout), target :: sites(:)
    type(bc_in_type)        , intent(in)            :: bc_in(:)

    integer  :: s

    do s = 1, size(sites,dim=1)
       sites(s)%snow_depth = bc_in(s)%snow_depth_si * bc_in(s)%frac_sno_eff_si
    end do

    return
  end subroutine UpdateFatesAvgSnowDepth


  ! =====================================================================================

  subroutine leaf_area_profile( currentSite )

    ! -----------------------------------------------------------------------------------
    ! This subroutine calculates how leaf and stem areas are distributed
    ! in vertical and horizontal space.
    !
    ! The following cohort level diagnostics are updated here:
    !
    ! currentCohort%treelai    ! LAI per unit crown area  (m2/m2)
    ! currentCohort%treesai    ! SAI per unit crown area  (m2/m2)
    ! currentCohort%NV         ! The number of discrete vegetation
    !                          ! layers needed to describe this crown
    !
    ! The following patch level diagnostics are updated here:
    !
    ! currentPatch%canopy_layer_tlai(cl)   ! total leaf area index of canopy layer
    ! currentPatch%nleaf(cl,ft)             ! number of vegetation layers needed
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
    ! -----------------------------------------------------------------------------------

    ! !USES:

    use EDtypesMod           , only : area, heightmax, n_height_bins
    use EDParamsMod,           only : dlower_vai,dinc_vai

    !
    ! !ARGUMENTS
    type(ed_site_type)     , intent(inout) :: currentSite


    !
    ! !LOCAL VARIABLES:
    type (fates_patch_type)  , pointer :: cpatch
    type (fates_cohort_type) , pointer :: currentCohort
    real(r8) :: remainder                !Thickness of layer at bottom of canopy.
    real(r8) :: fleaf                    ! fraction of cohort incepting area that is leaves.
    integer  :: ft                       ! Plant functional type index.
    integer  :: iv                       ! Vertical leaf layer index
    integer  :: cl                       ! Canopy layer index
 
    real(r8) :: fraction_exposed         ! how much of this layer is not covered by snow?
    real(r8) :: frac_canopy(N_HEIGHT_BINS) ! amount of canopy in each height class
    real(r8) :: minh(N_HEIGHT_BINS)        ! minimum height in height class (m)
    real(r8) :: maxh(N_HEIGHT_BINS)        ! maximum height in height class (m)
    real(r8) :: dh                       ! vertical detph of height class (m)
    real(r8) :: min_cheight                ! bottom of cohort canopy  (m)
    real(r8) :: max_cheight                ! top of cohort canopy      (m)
    real(r8) :: elai_layer,tlai_layer    ! leaf area per canopy area
    real(r8) :: esai_layer,tsai_layer    ! stem area per canopy area
    real(r8) :: vai_top,vai_bot          ! integrated top down veg area index at boundary of layer
    real(r8) :: crown_depth              ! Current cohort's crown depth
    real(r8) :: layer_bottom_height,layer_top_height,lai,sai  ! Can be removed later
    !----------------------------------------------------------------------


    ! Here we are trying to generate a profile of leaf area, indexed by 'z' and by pft
    ! We assume that each point in the canopy recieved the light attenuated by the average
    ! leaf area index above it, irrespective of PFT identity...
    ! Each leaf is defined by how deep in the canopy it is, in terms of LAI units.  (FIX(RF,032414), GB)

    cpatch => currentSite%oldest_patch
    do while(associated(cpatch))
       
       cpatch%nleaf(:,:) = 0
       cpatch%canopy_layer_tlai(:)        = 0._r8
       ! This routine updates the %nleaf array
       call UpdatePatchLAI(cpatch)
          
          
       ! This call assesses if the large, dynamically allocated
       ! patch arrays need to be allocated for the first time,
       ! or resized
       call cpatch%ReAllocateDynamics()

       ! These calls NaN and zero the above mentioned arrays
       call cpatch%NanDynamics()
       call cpatch%ZeroDynamics()
       
       cpatch%canopy_mask(:,:) = 0

       ! TO-DO: NRAD HYPOTHETICALLY WOULDNT INCLUDE LAYERS
       ! UNDER THE SNOW, BUT WE DONT REALLY USE IT TO FILTER
       ! THEM OUT. CHECK THE CODE AND CONSIDER REMOVING NRAD
       ! ALTOGETHER (RGK 05-2024)
       cpatch%nrad(:,:) = cpatch%nleaf(:,:)
       
       ! ------------------------------------------------------------------------------
       ! It is remotely possible that in deserts we will not have any canopy
       ! area, ie not plants at all...
       ! ------------------------------------------------------------------------------

       if_any_canopy_area: if (cpatch%total_canopy_area > nearzero ) then

          ! -----------------------------------------------------------------------------
          ! Standard canopy layering model.
          ! Go through all cohorts and add their leaf area
          ! and canopy area to the accumulators.
          ! -----------------------------------------------------------------------------

          currentCohort => cpatch%shortest
          do while(associated(currentCohort))
             ft = currentCohort%pft
             cl = currentCohort%canopy_layer

             ! ----------------------------------------------------------------
             ! How much of each tree is stem area index? Assuming that there is
             ! This may indeed be zero if there is a sensecent grass
             ! ----------------------------------------------------------------

             ! preserve_b4b will be removed soon. This is kept here to prevent
             ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23)
             if_preserve_b4b: if(preserve_b4b) then

                lai = currentCohort%treelai * currentCohort%c_area/cpatch%total_canopy_area
                sai = currentCohort%treesai * currentCohort%c_area/cpatch%total_canopy_area
                if( (currentCohort%treelai+currentCohort%treesai) > nearzero)then

                   ! See issue: https://github.com/NGEET/fates/issues/899
                   ! fleaf = currentCohort%treelai / (currentCohort%treelai + currentCohort%treesai)
                   fleaf = lai / (lai+sai)
                else
                   fleaf = 0._r8
                endif

                cpatch%nrad(cl,ft) = cpatch%nleaf(cl,ft)

                if (cpatch%nrad(cl,ft) > nlevleaf ) then
                   write(fates_log(), *) 'Number of radiative leaf layers is larger'
                   write(fates_log(), *) ' than the maximum allowed.'
                   write(fates_log(), *) ' cl: ',cl
                   write(fates_log(), *) ' ft: ',ft
                   write(fates_log(), *) ' nlevleaf: ',nlevleaf
                   write(fates_log(), *) ' cpatch%nrad(cl,ft): ', cpatch%nrad(cl,ft)
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if

                !---~---
                !   Find current crown depth using the allometric function.
                !---~---
                call CrownDepth(currentCohort%height,currentCohort%pft,crown_depth)
                !---~---


                ! --------------------------------------------------------------------------
                ! Whole layers.  Make a weighted average of the leaf area in each layer
                ! before dividing it by the total area. Fill up layer for whole layers.
                ! --------------------------------------------------------------------------

                do iv = 1,currentCohort%NV

                   ! This loop builds the arrays that define the effective (not snow covered)
                   ! and total (includes snow covered) area indices for leaves and stems
                   ! We calculate the absolute elevation of each layer to help determine if the layer
                   ! is obscured by snow.

                   layer_top_height = currentCohort%height - &
                        ( real(iv-1,r8)/currentCohort%NV * crown_depth )

                   layer_bottom_height = currentCohort%height - &
                        ( real(iv,r8)/currentCohort%NV * crown_depth )

                   fraction_exposed = 1.0_r8
                   if(currentSite%snow_depth  > layer_top_height)then
                      fraction_exposed = 0._r8
                   endif
                   if(currentSite%snow_depth < layer_bottom_height)then
                      fraction_exposed = 1._r8
                   endif
                   if(currentSite%snow_depth >= layer_bottom_height .and. &
                        currentSite%snow_depth <= layer_top_height) then !only partly hidden...
                      fraction_exposed =  1._r8 - max(0._r8,(min(1.0_r8,(currentSite%snow_depth -layer_bottom_height)/ &
                           (layer_top_height-layer_bottom_height ))))
                   endif

                   if(iv==currentCohort%NV) then
                      remainder = (currentCohort%treelai + currentCohort%treesai) - dlower_vai(iv)
                   else
                      remainder = dinc_vai(iv)
                   end if

                   cpatch%tlai_profile(cl,ft,iv) = cpatch%tlai_profile(cl,ft,iv) + &
                        remainder * fleaf * currentCohort%c_area/cpatch%total_canopy_area

                   cpatch%elai_profile(cl,ft,iv) = cpatch%elai_profile(cl,ft,iv) + &
                        remainder * fleaf * currentCohort%c_area/cpatch%total_canopy_area * &
                        fraction_exposed

                   cpatch%tsai_profile(cl,ft,iv) = cpatch%tsai_profile(cl,ft,iv) + &
                        remainder * (1._r8 - fleaf) * currentCohort%c_area/cpatch%total_canopy_area

                   cpatch%esai_profile(cl,ft,iv) = cpatch%esai_profile(cl,ft,iv) + &
                        remainder * (1._r8 - fleaf) * currentCohort%c_area/cpatch%total_canopy_area * &
                        fraction_exposed

                   cpatch%canopy_area_profile(cl,ft,iv) = cpatch%canopy_area_profile(cl,ft,iv) + &
                        currentCohort%c_area/cpatch%total_canopy_area


                end do

             else !if_preserve_b4b

                do iv = 1,currentCohort%NV

                   call VegAreaLayer(currentCohort%treelai,     &
                        currentCohort%treesai,                  &
                        currentCohort%height,                   &
                        iv,currentCohort%nv,currentCohort%pft,  &
                        currentSite%snow_depth,                    &
                        vai_top,vai_bot,                          &
                        elai_layer,esai_layer,tlai_layer,tsai_layer)


                   cpatch%tlai_profile(cl,ft,iv) = cpatch%tlai_profile(cl,ft,iv) + &
                        tlai_layer * currentCohort%c_area/cpatch%total_canopy_area

                   cpatch%elai_profile(cl,ft,iv) = cpatch%elai_profile(cl,ft,iv) + &
                        elai_layer * currentCohort%c_area/cpatch%total_canopy_area

                   cpatch%tsai_profile(cl,ft,iv) = cpatch%tsai_profile(cl,ft,iv) + &
                        tsai_layer * currentCohort%c_area/cpatch%total_canopy_area

                   cpatch%esai_profile(cl,ft,iv) = cpatch%esai_profile(cl,ft,iv) + &
                        esai_layer * currentCohort%c_area/cpatch%total_canopy_area

                   cpatch%canopy_area_profile(cl,ft,iv) = cpatch%canopy_area_profile(cl,ft,iv) + &
                        currentCohort%c_area/cpatch%total_canopy_area

                end do

             end if if_preserve_b4b

             currentCohort => currentCohort%taller

          enddo !cohort

          ! --------------------------------------------------------------------------

          ! If there is an upper-story, the top canopy layer
          ! should have a value of exactly 1.0 in its top leaf layer
          ! --------------------------------------------------------------------------

          if ( (cpatch%NCL_p > 1) .and. &
               (sum(cpatch%canopy_area_profile(1,:,1)) < 0.9999 )) then
             write(fates_log(), *) 'FATES: canopy_area_profile was less than 1 at the canopy top'
             write(fates_log(), *) 'cl: ',1
             write(fates_log(), *) 'iv: ',1
             write(fates_log(), *) 'sum(cpatch%canopy_area_profile(1,:,1)): ', &
                  sum(cpatch%canopy_area_profile(1,:,1))
             currentCohort => cpatch%shortest
             do while(associated(currentCohort))
                if(currentCohort%canopy_layer==1)then
                   write(fates_log(), *) 'FATES: cohorts',currentCohort%dbh,currentCohort%c_area, &
                        cpatch%total_canopy_area,cpatch%area
                   write(fates_log(), *) 'ED: fracarea', currentCohort%pft, &
                        currentCohort%c_area/cpatch%total_canopy_area
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

          do cl = 1,cpatch%NCL_p
             do iv = 1,cpatch%nleaf(cl,ft)

                if( debug .and. sum(cpatch%canopy_area_profile(cl,:,iv)) > 1.0001_r8 ) then

                   write(fates_log(), *) 'FATES: A canopy_area_profile exceeded 1.0'
                   write(fates_log(), *) 'cl: ',cl
                   write(fates_log(), *) 'iv: ',iv
                   write(fates_log(), *) 'sum(cpatch%canopy_area_profile(cl,:,iv)): ', &
                        sum(cpatch%canopy_area_profile(cl,:,iv))
                   currentCohort => cpatch%shortest
                   do while(associated(currentCohort))
                      if(currentCohort%canopy_layer==cl)then
                         write(fates_log(), *) 'FATES: cohorts in layer cl = ',cl, &
                              currentCohort%dbh,currentCohort%c_area, &
                              cpatch%total_canopy_area,cpatch%area
                         write(fates_log(), *) 'ED: fracarea', currentCohort%pft, &
                              currentCohort%c_area/cpatch%total_canopy_area
                      endif
                      currentCohort => currentCohort%taller
                   enddo !currentCohort
                   call endrun(msg=errMsg(sourcefile, __LINE__))
                end if
             end do

             do ft = 1,numpft
                do iv = 1,cpatch%nleaf(cl,ft)

                   if( cpatch%canopy_area_profile(cl,ft,iv) > nearzero ) then

                      cpatch%tlai_profile(cl,ft,iv) = cpatch%tlai_profile(cl,ft,iv) / &
                           cpatch%canopy_area_profile(cl,ft,iv)

                      cpatch%tsai_profile(cl,ft,iv) = cpatch%tsai_profile(cl,ft,iv) / &
                           cpatch%canopy_area_profile(cl,ft,iv)

                      cpatch%elai_profile(cl,ft,iv) = cpatch%elai_profile(cl,ft,iv) / &
                           cpatch%canopy_area_profile(cl,ft,iv)

                      cpatch%esai_profile(cl,ft,iv) = cpatch%esai_profile(cl,ft,iv) / &
                           cpatch%canopy_area_profile(cl,ft,iv)
                   end if

                enddo

             enddo
          enddo

          ! --------------------------------------------------------------------------
          ! Set the mask that identifies which PFT x can-layer combinations have
          ! scattering elements in them for radiation.
          ! RGK: I'm not sure we need nrad ... I can't see a scenario where
          !      canopy_area_profile for these layers is not >0 for layers in ncan ...
          !      Leaving this for the time being.
          ! --------------------------------------------------------------------------
          
          cpatch%canopy_mask(:,:) = 0
          ! preserve_b4b will be removed soon. This is kept here to prevent
          ! round off errors in the baseline tests for the two-stream code (RGK 12-27-23)
          if(preserve_b4b) then
             do cl = 1,cpatch%NCL_p
                do ft = 1,numpft
                   do  iv = 1, cpatch%nrad(cl,ft)
                      if(cpatch%canopy_area_profile(cl,ft,iv) > 0._r8)then
                         cpatch%canopy_mask(cl,ft) = 1
                      endif
                   end do !iv
                end do
             end do
          else
             do cl = 1,cpatch%NCL_p
                do ft = 1,numpft
                   if(cpatch%canopy_area_profile(cl,ft,1) > 0._r8 ) cpatch%canopy_mask(cl,ft) = 1
                end do
             end do
          end if

             
       end if if_any_canopy_area

       cpatch => cpatch%younger
    enddo !patch

    
    
    return
  end subroutine leaf_area_profile

  ! ======================================================================================

  subroutine update_hlm_dynamics(nsites,sites,fcolumn,bc_out)

    ! ----------------------------------------------------------------------------------
    ! The purpose of this routine is to package output boundary conditions related
    ! to vegetation coverage to the host land model.
    ! ----------------------------------------------------------------------------------

    use EDTypesMod        , only : ed_site_type, AREA
    use FatesPatchMod,      only : fates_patch_type
    use FatesInterfaceTypesMod , only : bc_out_type

    !
    ! !ARGUMENTS
    integer,            intent(in)            :: nsites
    type(ed_site_type), intent(inout), target :: sites(nsites)
    integer,            intent(in)            :: fcolumn(nsites)
    type(bc_out_type),  intent(inout)         :: bc_out(nsites)

    ! Locals
    type (fates_cohort_type) , pointer :: currentCohort
    integer :: s, ifp, c, p
    type (fates_patch_type)  , pointer :: currentPatch
    real(r8) :: bare_frac_area
    real(r8) :: total_patch_area
    real(r8) :: total_canopy_area
    real(r8) :: total_patch_leaf_stem_area
    real(r8) :: weight  ! Weighting for cohort variables in patch
    
    do s = 1,nsites

       total_patch_area = 0._r8
       total_canopy_area = 0._r8
       bc_out(s)%canopy_fraction_pa(:) = 0._r8
       bc_out(s)%dleaf_pa(:) = 0._r8
       bc_out(s)%z0m_pa(:) = 0._r8
       bc_out(s)%displa_pa(:) = 0._r8
       
       currentPatch => sites(s)%oldest_patch
       c = fcolumn(s)
       do while(associated(currentPatch))

          ifp = currentPatch%patchno
          if_bare: if(currentPatch%nocomp_pft_label.ne.nocomp_bareground)then  ! ignore the bare-ground-PFT patch entirely for these BC outs

             if ( currentPatch%total_canopy_area-currentPatch%area > 0.000001_r8 ) then
                if(debug)then
                   write(fates_log(),*) 'ED: canopy area bigger than area', &
                        currentPatch%total_canopy_area ,currentPatch%area
                end if
                currentPatch%total_canopy_area = currentPatch%area
             endif

             if (associated(currentPatch%tallest)) then
                bc_out(s)%htop_pa(ifp) = currentPatch%tallest%height
             else
                ! FIX(RF,040113) - should this be a parameter for the minimum possible vegetation height?
                bc_out(s)%htop_pa(ifp) = 0.1_r8
             endif

             bc_out(s)%hbot_pa(ifp) = max(0._r8, min(0.2_r8, bc_out(s)%htop_pa(ifp)- 1.0_r8))

             ! Use canopy-only crown area weighting for all cohorts in the patch to define the characteristic
             ! Roughness length and displacement height used by the HLM
             ! use total LAI + SAI to weight the leaft characteristic dimension
             ! Avoid this if running in satellite phenology mode
             ! ----------------------------------------------------------------------------

             if (currentPatch%total_canopy_area > nearzero) then
                currentCohort => currentPatch%shortest
                do while(associated(currentCohort))
                   if (currentCohort%canopy_layer .eq. 1) then
                      weight = min(1.0_r8,currentCohort%c_area/currentPatch%total_canopy_area)
                      bc_out(s)%z0m_pa(ifp) = bc_out(s)%z0m_pa(ifp) + &
                           EDPftvarcon_inst%z0mr(currentCohort%pft) * currentCohort%height * weight
                      bc_out(s)%displa_pa(ifp) = bc_out(s)%displa_pa(ifp) + &
                           EDPftvarcon_inst%displar(currentCohort%pft) * currentCohort%height * weight
                   endif
                   currentCohort => currentCohort%taller
                end do

                ! for lai, scale to total LAI + SAI in patch.  first add up all the LAI and SAI in the patch
                total_patch_leaf_stem_area = 0._r8
                currentCohort => currentPatch%shortest
                do while(associated(currentCohort))
                   total_patch_leaf_stem_area = total_patch_leaf_stem_area + &
                        (currentCohort%treelai + currentCohort%treesai) * currentCohort%c_area
                   currentCohort => currentCohort%taller
                end do

                ! make sure there is some leaf and stem area
                if (total_patch_leaf_stem_area > nearzero) then
                   currentCohort => currentPatch%shortest
                   do while(associated(currentCohort))
                      ! weight dleaf by the relative totals of leaf and stem area
                      weight = (currentCohort%treelai + currentCohort%treesai) * currentCohort%c_area / total_patch_leaf_stem_area
                      bc_out(s)%dleaf_pa(ifp) = bc_out(s)%dleaf_pa(ifp) + &
                           EDPftvarcon_inst%dleaf(currentCohort%pft) * weight
                      currentCohort => currentCohort%taller
                   end do
                else
                   ! dummy case
                   bc_out(s)%dleaf_pa(ifp)  = EDPftvarcon_inst%dleaf(1)
                endif
             else
                ! if no canopy, then use dummy values (first PFT) of aerodynamic properties
                bc_out(s)%z0m_pa(ifp)    = EDPftvarcon_inst%z0mr(1) * bc_out(s)%htop_pa(ifp)
                bc_out(s)%displa_pa(ifp) = EDPftvarcon_inst%displar(1) * bc_out(s)%htop_pa(ifp)
                bc_out(s)%dleaf_pa(ifp)  = EDPftvarcon_inst%dleaf(1)
             endif
             ! -----------------------------------------------------------------------------

             ! We are assuming here that grass is all located underneath tree canopies.
             ! The alternative is to assume it is all spatial distinct from tree canopies.
             ! In which case, the bare area would have to be reduced by the grass area...
             ! currentPatch%total_canopy_area/currentPatch%area is fraction of this patch cover by plants
             ! currentPatch%area/AREA is the fraction of the soil covered by this patch.

             if(currentPatch%area.gt.0.0_r8)then
                bc_out(s)%canopy_fraction_pa(ifp) = &
                     min(1.0_r8,currentPatch%total_canopy_area/currentPatch%area)*(currentPatch%area/AREA)
             else
                bc_out(s)%canopy_fraction_pa(ifp) = 0.0_r8
             endif

             bare_frac_area = (1.0_r8 - min(1.0_r8,currentPatch%total_canopy_area/currentPatch%area)) * &
                  (currentPatch%area/AREA)

             total_patch_area = total_patch_area + bc_out(s)%canopy_fraction_pa(ifp) + bare_frac_area

             total_canopy_area = total_canopy_area + bc_out(s)%canopy_fraction_pa(ifp)

             bc_out(s)%nocomp_pft_label_pa(ifp) = currentPatch%nocomp_pft_label

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

          else  ! nocomp or SP, and currentPatch%nocomp_pft_label .eq. 0

             total_patch_area = total_patch_area + currentPatch%area/AREA

          end if if_bare
          currentPatch => currentPatch%younger
       end do

       ! Apply patch and canopy area corrections
       ! If the difference is above reasonable math precision, apply a fix
       ! If the difference is way above reasonable math precision, gracefully exit

       if(abs(total_patch_area-1.0_r8) > rsnbl_math_prec ) then

          if(abs(total_patch_area-1.0_r8) > 1.0e-8_r8 )then
             write(fates_log(),*) 'total area is wrong in update_hlm_dynamics',total_patch_area
             call endrun(msg=errMsg(sourcefile, __LINE__))
          end if

          if(debug) then
             write(fates_log(),*) 'imprecise patch areas in update_hlm_dynamics',total_patch_area
          end if

          currentPatch => sites(s)%oldest_patch
          do while(associated(currentPatch))
             ifp = currentPatch%patchno
             if(currentPatch%nocomp_pft_label.ne.nocomp_bareground)then ! for vegetated patches only
                bc_out(s)%canopy_fraction_pa(ifp) = bc_out(s)%canopy_fraction_pa(ifp)/total_patch_area
             endif ! veg patch
             currentPatch => currentPatch%younger
          end do

       endif

       ! If running hydro, perform a final check to make sure that we
       ! have conserved water. Since this is the very end of the dynamics
       ! cycle. No water should had been added or lost to the site during dynamics.
       ! With growth and death, we may have shuffled it around.
       ! For recruitment, we initialized their water, but flagged them
       ! to not be included in the site level balance yet, for they
       ! will demand the water for their initialization on the first hydraulics time-step

       if (hlm_use_planthydro.eq.itrue) then
          call UpdateH2OVeg(sites(s),bc_out(s),bc_out(s)%plant_stored_h2o_si,1)
       end if

       ! Pass FATES Harvested C to bc_out.
       call UpdateHarvestC(sites(s),bc_out(s))

    end do

    ! This call to RecruitWaterStorage() makes an accounting of
    ! how much water is used to intialize newly recruited plants.
    ! However, it does not actually move water from the soil or create
    ! a flux, it is just accounting for diagnostics purposes.  The water
    ! will not actually be moved until the beginning of the first hydraulics
    ! call during the fast timestep sequence

    if (hlm_use_planthydro.eq.itrue) then
       call RecruitWaterStorage(nsites,sites)
    end if

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
    type(fates_patch_type),intent(in), target :: cpatch
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
    type(fates_patch_type),intent(inout), target   :: currentPatch
    real(r8),intent(in)                         :: site_spread
    integer,intent(in)                          :: layer_index
    real(r8),intent(inout)                      :: layer_area

    type(fates_cohort_type), pointer :: currentCohort
    
    
    layer_area = 0.0_r8
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort))
       call carea_allom(currentCohort%dbh,currentCohort%n,site_spread, &
            currentCohort%pft,currentCohort%crowndamage, currentCohort%c_area)
       if (currentCohort%canopy_layer .eq. layer_index) then
          layer_area = layer_area + currentCohort%c_area
       end if
       currentCohort => currentCohort%shorter
    enddo
    return
   end subroutine CanopyLayerArea
  
  ! ===============================================================================================

  subroutine UpdatePatchLAI(currentPatch)

   ! --------------------------------------------------------------------------------------------
   ! This subroutine works through the current patch cohorts and updates the canopy_layer_tlai
   ! and related variables
   ! ---------------------------------------------------------------------------------------------

   ! Arguments
   type(fates_patch_type),intent(inout), target   :: currentPatch

   ! Local Variables
   type(fates_cohort_type), pointer :: currentCohort
   integer  :: cl                                  ! Canopy layer index
   integer  :: ft                                  ! Plant functional type index
   
   ! Calculate LAI of layers above.  Because it is possible for some understory cohorts
   ! to be taller than cohorts in the top canopy layer, we must iterate through the 
   ! patch by canopy layer first.  Given that canopy_layer_tlai is a patch level variable
   ! we could iterate through each cohort in any direction as long as we go down through
   ! the canopy layers.

   canopyloop: do cl = 1,nclmax
      currentCohort => currentPatch%tallest
      cohortloop: do while(associated(currentCohort))

         ! Only update the current cohort tree lai if lai of the above layers have been calculated
         if (currentCohort%canopy_layer .eq. cl) then
            
            ft     = currentCohort%pft
            ! Update the cohort level lai and related variables
            call UpdateCohortLAI(currentCohort,currentPatch%canopy_layer_tlai,  &
                 currentPatch%total_canopy_area)
            
            ! Update the number of number of vegetation layers
            currentPatch%nleaf(cl,ft) = max(currentPatch%nleaf(cl,ft),currentCohort%NV)

            ! Update the patch canopy layer tlai (LAI per canopy area)
            currentPatch%canopy_layer_tlai(cl) = currentPatch%canopy_layer_tlai(cl) +  &
                 currentCohort%treelai *currentCohort%c_area/currentPatch%total_canopy_area
            
         end if
         currentCohort => currentCohort%shorter

      end do cohortloop
   end do canopyloop

  end subroutine UpdatePatchLAI
  ! ===============================================================================================
  
  subroutine UpdateCohortLAI(currentCohort, canopy_layer_tlai, total_canopy_area)
   
   ! Update LAI and related variables for a given cohort
   
   ! Arguments
   type(fates_cohort_type),intent(inout), target   :: currentCohort
   real(r8), intent(in) :: canopy_layer_tlai(nclmax)  ! total leaf area index of each canopy layer
   real(r8), intent(in) :: total_canopy_area                  ! either patch%total_canopy_area or patch%area
   
   ! Local variables
   real(r8) :: leaf_c                              ! leaf carbon [kg]
   real(r8) :: treesai                             ! stem area index within crown m2/m2
   
   ! Obtain the leaf carbon
   leaf_c = currentCohort%prt%GetState(leaf_organ,carbon12_element)

   ! Note that tree_lai has an internal check on the canopy location
   call  tree_lai_sai(leaf_c, currentCohort%pft, currentCohort%c_area, currentCohort%n,           &
          currentCohort%canopy_layer, canopy_layer_tlai, currentCohort%vcmax25top, currentCohort%dbh, currentCohort%crowndamage,          &
          currentCohort%canopy_trim, currentCohort%efstem_coh, 4, currentCohort%treelai, treesai )

   ! Do not update stem area index of SP vegetation
   if (hlm_use_sp .eq. ifalse) then
      currentCohort%treesai = treesai
   end if
   
   ! Number of actual vegetation layers in this cohort's crown
   currentCohort%nv = GetNVegLayers(currentCohort%treelai+currentCohort%treesai)
   
  end subroutine UpdateCohortLAI
  
  ! ===============================================================================================

  function NumCanopyLayers(currentPatch) result(z)

    ! --------------------------------------------------------------------------------------------
    ! Calculate the number of canopy layers in this patch.
    ! This simple call only determines total layering by querying the cohorts
    ! which layer they are in, it doesn't do any size evaluation.
    ! --------------------------------------------------------------------------------------------

    type(fates_patch_type)          :: currentPatch
    type(fates_cohort_type),pointer :: currentCohort

    integer :: z
    real(r8) :: c_area
    real(r8) :: arealayer

    z = 1
    currentCohort => currentPatch%tallest
    do while (associated(currentCohort))
       z = max(z,currentCohort%canopy_layer)
       currentCohort => currentCohort%shorter
    enddo

  end function NumCanopyLayers

end module EDCanopyStructureMod
