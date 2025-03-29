module FatesLandUseChangeMod

  ! Controls the transfer and initialization of patch structure to land use types

  use FatesGlobals              , only : fates_log
  use FatesConstantsMod         , only : primaryland, secondaryland, pastureland, rangeland, cropland
  use FatesConstantsMod         , only : n_landuse_cats
  use FatesConstantsMod         , only : nearzero
  use FatesGlobals              , only : endrun => fates_endrun
  use FatesConstantsMod         , only : r8 => fates_r8
  use FatesConstantsMod         , only : itrue, ifalse
  use FatesConstantsMod         , only : fates_unset_int
  use FatesConstantsMod         , only : years_per_day
  use FatesInterfaceTypesMod    , only : bc_in_type
  use FatesInterfaceTypesMod    , only : hlm_use_luh
  use FatesInterfaceTypesMod    , only : hlm_num_luh2_states
  use FatesInterfaceTypesMod    , only : hlm_num_luh2_transitions
  use FatesInterfaceTypesMod    , only : hlm_use_potentialveg
  use FatesUtilsMod             , only : FindIndex
  use EDTypesMod                , only : area_site => area

  ! CIME globals
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  
  !
  implicit none
  private

  character(len=*), parameter :: sourcefile = __FILE__

  public :: GetLanduseTransitionRates
  public :: GetLanduseChangeRules
  public :: GetLUHStatedata
  public :: GetInitLanduseTransitionRates
  public :: GetInitLanduseHarvestRate
  public :: FatesGrazing

  ! module data
  integer, parameter :: max_luh2_types_per_fates_lu_type = 5

  ! Define the mapping from the luh2 state names to the aggregated fates land use categories
  type :: luh2_fates_lutype_map

     character(len=5), dimension(12) :: state_names = &
                   [character(len=5) ::  'primf','primn','secdf','secdn', &
                                         'pastr','range', 'urban', &
                                         'c3ann','c4ann','c3per','c4per','c3nfx']
     integer, dimension(12) :: landuse_categories = &
                                [primaryland, primaryland, secondaryland, secondaryland, &
                                 pastureland, rangeland, fates_unset_int, &
                                 cropland, cropland, cropland, cropland, cropland]

     contains

       procedure :: GetIndex => GetLUCategoryFromStateName

  end type luh2_fates_lutype_map


  ! 03/10/2023 Created By Charlie Koven
  ! ============================================================================

contains

  ! ============================================================================
  subroutine GetLanduseTransitionRates(bc_in, min_allowed_landuse_fraction, landuse_transition_matrix, &
       landuse_vector_gt_min)


    ! The purpose of this routine is to ingest the land use transition rate information that the host
    ! model has read in from a dataset,aggregate land use types to those being used in the simulation,
    ! and output a transition matrix that can be used to drive patch disturbance rates.

    ! !ARGUMENTS:
    type(bc_in_type) , intent(in) :: bc_in
    real(r8), intent(in)          :: min_allowed_landuse_fraction
    real(r8), intent(inout)       :: landuse_transition_matrix(n_landuse_cats, n_landuse_cats)  ! [m2/m2/day]
    logical,  intent(inout)       :: landuse_vector_gt_min(n_landuse_cats)

    ! !LOCAL VARIABLES:
    type(luh2_fates_lutype_map) :: lumap
    integer :: i_donor, i_receiver, i_luh2_transitions, i_luh2_states, i_urban
    character(5) :: donor_name, receiver_name
    character(14) :: transition_name
    real(r8) :: urban_fraction
    real(r8) :: temp_vector(hlm_num_luh2_transitions)
    logical  :: modified_flag
    real(r8) :: state_vector(n_landuse_cats)  ! [m2/m2]
    integer  :: i_lu

    ! zero the transition matrix and the urban fraction
    landuse_transition_matrix(:,:) = 0._r8
    urban_fraction = 0._r8

    ! if we are using potential veg only, then keep all transitions equal to zero.
    if (hlm_use_potentialveg .eq. ifalse) then

       ! Check the LUH data incoming to see if any of the transitions are NaN
       temp_vector = bc_in%hlm_luh_transitions
       call CheckLUHData(temp_vector,modified_flag)
       if (.not. modified_flag) then
          ! identify urban fraction so that it can be factored into the land use state output
          urban_fraction = bc_in%hlm_luh_states(FindIndex(bc_in%hlm_luh_state_names,'urban'))
       end if

       !! TODO: may need some logic here to ask whether or not ot perform land use change on this timestep.
       !! current code occurs every day. If not doing transition every day, need to update units.

       transitions_loop: do i_luh2_transitions = 1, hlm_num_luh2_transitions

          ! transition names are written in form xxxxx_to_yyyyy where x and y are donor and receiver state names
          transition_name = bc_in%hlm_luh_transition_names(i_luh2_transitions)
          donor_name = transition_name(1:5)
          receiver_name = transition_name(10:14)

          ! Get the fates land use type index associated with the luh2 state types
          i_donor= lumap%GetIndex(donor_name)
          i_receiver = lumap%GetIndex(receiver_name)

          ! Avoid transitions with 'urban' as those are handled seperately
          ! Also ignore diagonal elements of transition matrix.
          if (.not.(i_donor .eq. fates_unset_int .or. i_receiver .eq. fates_unset_int .or. &
               i_donor .eq. i_receiver)) then
             landuse_transition_matrix(i_donor,i_receiver) = &
                  landuse_transition_matrix(i_donor,i_receiver) +  temp_vector(i_luh2_transitions) &
                  * years_per_day / (1._r8 - urban_fraction)

          end if
       end do transitions_loop

       ! zero all transitions where the receiving land use type state vector is less than the minimum allowed,
       ! and otherwise if this is the first timestep where the minimum was exceeded,
       ! then apply all transitions from primary to this type and reset the flag
       ! note that the flag resetting should not happen for secondary lands, as this is handled in the
       ! logging logic
       call GetLUHStatedata(bc_in, state_vector)
       do i_lu = secondaryland, n_landuse_cats
          if ( state_vector(i_lu) .le. min_allowed_landuse_fraction ) then
             landuse_transition_matrix(:,i_lu) = 0._r8
          else if ((.not. landuse_vector_gt_min(i_lu)) .and. (i_lu .ne. secondaryland)) then
             landuse_transition_matrix(:,i_lu) = 0._r8
             landuse_transition_matrix(primaryland,i_lu) = state_vector(i_lu)
             landuse_vector_gt_min(i_lu) = .true.
          end if
       end do
    end if

  end subroutine GetLanduseTransitionRates

  !----------------------------------------------------------------------------------------------------

  function GetLUCategoryFromStateName(this, state_name) result(landuse_category)

    class(luh2_fates_lutype_map) :: this
    character(len=5), intent(in) :: state_name
    integer :: landuse_category
    integer :: index_statename

    index_statename = FindIndex(this%state_names,state_name)

    ! Check that the result from the landuse_categories is not zero, which indicates that no
    ! match was found.
    if (index_statename .eq. 0) then
       write(fates_log(),*) 'The input state name from the HLM does not match the FATES landuse state name options'
       write(fates_log(),*) 'input state name: ', state_name
       write(fates_log(),*) 'state name options: ', this%state_names
       call endrun(msg=errMsg(sourcefile, __LINE__))
    else
       landuse_category = this%landuse_categories(index_statename)
    end if

  end function GetLUCategoryFromStateName

  !----------------------------------------------------------------------------------------------------

  subroutine GetLanduseChangeRules(clearing_matrix)

    ! the purpose of this is to define a ruleset for when to clear the vegetation in transitioning
    ! from one land use type to another

    logical, intent(out) :: clearing_matrix(n_landuse_cats,n_landuse_cats)
    
    ! default value of ruleset 4 above means that plants are not cleared during land use change
    ! transitions to rangeland, whereas plants are cleared in transitions to pasturelands and croplands.
    integer, parameter    :: ruleset = 4   ! ruleset to apply from table 1 of Ma et al (2020)
    ! https://doi.org/10.5194/gmd-13-3203-2020

    ! clearing matrix applies from the donor to the receiver land use type of the newly-transferred
    ! patch area values of clearing matrix: false => do not clear; true => clear

    clearing_matrix(:,:) = .false.

    select case(ruleset)

    case(1)

       ! note that this ruleset isnt exactly what is in Ma et al. rulesets 1 and 2, because FATES
       ! does not make the distinction between forested and non-forested lands from a land use/land
       ! cover perspective.

       clearing_matrix(:,cropland) = .true.
       clearing_matrix(:,pastureland) = .true.
       clearing_matrix(primaryland,rangeland) = .true.
       clearing_matrix(secondaryland,rangeland) = .true.

    case(2)

       ! see comment on number 1 above
       clearing_matrix(:,cropland) = .true.
       clearing_matrix(primaryland,pastureland) = .true.
       clearing_matrix(secondaryland,pastureland) = .true.
       clearing_matrix(primaryland,rangeland) = .true.
       clearing_matrix(secondaryland,rangeland) = .true.

    case(3)

       clearing_matrix(:,cropland) = .true.
       clearing_matrix(:,pastureland) = .true.
       clearing_matrix(:,rangeland) = .true.

    case(4)

       clearing_matrix(:,cropland) = .true.
       clearing_matrix(:,pastureland) = .true.
       clearing_matrix(:,rangeland) = .false.

    case(5)

       clearing_matrix(:,cropland) = .true.
       clearing_matrix(:,pastureland) = .false.
       clearing_matrix(:,rangeland) = .true.

    case(6)

       clearing_matrix(:,cropland) = .true.
       clearing_matrix(:,pastureland) = .false.
       clearing_matrix(:,rangeland) = .false.

    case(7)

       clearing_matrix(:,cropland) = .false.
       clearing_matrix(:,pastureland) = .true.
       clearing_matrix(:,rangeland) = .true.

    case(8)

       clearing_matrix(:,cropland) = .false.
       clearing_matrix(:,pastureland) = .true.
       clearing_matrix(:,rangeland) = .false.

    case(9)

       clearing_matrix(:,cropland) = .false.
       clearing_matrix(:,pastureland) = .false.
       clearing_matrix(:,rangeland) = .true.

    case default

       write(fates_log(),*) 'unknown clearing ruleset?'
       write(fates_log(),*) 'ruleset: ', ruleset
       call endrun(msg=errMsg(sourcefile, __LINE__))

    end select

  end subroutine GetLanduseChangeRules

  !----------------------------------------------------------------------------------------------------

  subroutine GetLUHStatedata(bc_in, state_vector)

    type(bc_in_type) , intent(in) :: bc_in
    real(r8), intent(out) :: state_vector(n_landuse_cats)  ! [m2/m2]

    ! LOCALS
    type(luh2_fates_lutype_map) :: lumap
    real(r8) :: temp_vector(hlm_num_luh2_states)  ! [m2/m2]
    real(r8) :: urban_fraction
    integer  :: i_luh2_states
    integer  :: ii
    character(5) :: state_name
    logical :: modified_flag

    ! zero state vector and urban fraction
    state_vector(:) = 0._r8
    urban_fraction = 0._r8

    if (hlm_use_potentialveg .eq. itrue) then
       state_vector(primaryland) = 1._r8
    else
       ! Check to see if the incoming state vector is NaN.
       temp_vector = bc_in%hlm_luh_states
       call CheckLUHData(temp_vector,modified_flag)
       if (.not. modified_flag) then
          ! identify urban fraction so that it can be factored into the land use state output
          urban_fraction = bc_in%hlm_luh_states(FindIndex(bc_in%hlm_luh_state_names,'urban'))
       end if

       ! loop over all states and add up the ones that correspond to a given fates land use type
       do i_luh2_states = 1, hlm_num_luh2_states

          ! Get the luh2 state name and determine fates aggregated land use
          ! type index from the state to lutype map
          state_name = bc_in%hlm_luh_state_names(i_luh2_states)
          ii = lumap%GetIndex(state_name)

          ! Avoid 'urban' states whose indices have been given unset values
          if (ii .ne. fates_unset_int) then
             state_vector(ii) = state_vector(ii) + &
                  temp_vector(i_luh2_states) / (1._r8 - urban_fraction)
          end if
       end do

       ! if all zeros, make all primary lands
       if ( sum(state_vector(:)) .gt. nearzero ) then

          ! check to ensure total area == 1, and correct if not
          if ( abs(sum(state_vector(:)) - 1._r8) .gt. nearzero ) then
             state_vector(:) = state_vector(:) / sum(state_vector(:))
          end if
       else
          state_vector(primaryland) = 1._r8
       endif
    end if

  end subroutine GetLUHStatedata

  !----------------------------------------------------------------------------------------------------

  subroutine CheckLUHData(luh_vector,modified_flag)

    use shr_infnan_mod   , only : isnan => shr_infnan_isnan

    real(r8), intent(inout) :: luh_vector(:)  ! [m2/m2]
    logical, intent(out)    :: modified_flag

    ! Check to see if the incoming luh2 vector is NaN.
    ! This suggests that there is a discepency where the HLM and LUH2 states
    ! there is vegetated ground. E.g. LUH2 data is missing for glacier-margin
    ! regions such as Antarctica. In this case, states should be Nan.  If so,
    ! set the current state to be all primary forest, and all transitions to be zero.
    ! If only a portion of the vector is NaN, there is something  amiss with
    ! the data, so end the run.

    modified_flag = .false.
    if (all(isnan(luh_vector))) then
       luh_vector(:) = 0._r8
       ! Check if this is a state vector, otherwise leave transitions as zero
       if (size(luh_vector) .eq. hlm_num_luh2_states) then
          luh_vector(primaryland) = 1._r8
       end if
       modified_flag = .true.
       !write(fates_log(),*) 'WARNING: land use state is all NaN;
       !setting state as all primary forest.' ! GL DIAG
    else if (any(isnan(luh_vector))) then
       if (any(.not. isnan(luh_vector))) then
          write(fates_log(),*) 'ERROR: land use vector has NaN'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine CheckLUHData


  subroutine GetInitLanduseHarvestRate(bc_in, min_allowed_landuse_fraction, harvest_rate, &
       landuse_vector_gt_min)

    ! the purpose of this subroutine is, only under the case where we are transitioning from a spinup
    ! run that did not have land use to a run that does, to apply the land-use changes needed to get
    ! to the state vector in a single daily instance. this is for the hrvest rate from primary lands,
    ! i.e. the transition from primary to secondary lands. thus instead of using the harvest dataset
    ! itself, it only uses the state vector for what land use compositoin we want to achieve, and log
    ! the forests accordingly.

    ! !ARGUMENTS:
    type(bc_in_type) , intent(in) :: bc_in
    real(r8), intent(in)          :: min_allowed_landuse_fraction
    real(r8), intent(out)         :: harvest_rate  ! [m2/ m2 / day]
    logical,  intent(inout)       :: landuse_vector_gt_min(n_landuse_cats)

    ! LOCALS
    real(r8) ::  state_vector(n_landuse_cats)  ! [m2/m2]
    
    call GetLUHStatedata(bc_in, state_vector)

    ! only do this if the state vector exceeds the minimum viable patch size, and if so, note that in the
    ! landuse_vector_gt_min flag (which will be coming in as .false. because of the use_potentialveg logic).
    if ( state_vector(secondaryland) .gt. min_allowed_landuse_fraction) then
       harvest_rate = state_vector(secondaryland)
       landuse_vector_gt_min(secondaryland) = .true.
    else
       harvest_rate = 0._r8
    endif

  end subroutine GetInitLanduseHarvestRate

  subroutine GetInitLanduseTransitionRates(bc_in, min_allowed_landuse_fraction, &
       landuse_transition_matrix, landuse_vector_gt_min)
    
    ! The purpose of this subroutine is, only under the case where we are transitioning from a spinup
    ! run that did not have land use to a run that does, to apply the land-use changes needed to get
    ! to the state vector in a single daily instance. This is for the transitions other than harvest,
    ! i.e. from primary lands to all other categories aside from secondary lands. 

    ! !ARGUMENTS:
    type(bc_in_type) , intent(in) :: bc_in
    real(r8), intent(in)          :: min_allowed_landuse_fraction
    real(r8), intent(inout)       :: landuse_transition_matrix(n_landuse_cats, n_landuse_cats)  ! [m2/m2/day]
    logical,  intent(inout)       :: landuse_vector_gt_min(n_landuse_cats)

    ! LOCALS
    real(r8) ::  state_vector(n_landuse_cats)  ! [m2/m2]
    integer  :: i

    landuse_transition_matrix(:,:) = 0._r8
    
    call GetLUHStatedata(bc_in, state_vector)

    ! only do this if the state vector exceeds the minimum viable patch size, and if so, note that
    ! in the landuse_vector_gt_min flag (which will be coming in as .false. because of the
    ! use_potentialveg logic).
    do i = secondaryland+1,n_landuse_cats
       if ( state_vector(i) .gt. min_allowed_landuse_fraction) then
          landuse_transition_matrix(primaryland,i) = state_vector(i)
          landuse_vector_gt_min(i) = .true.
       end if
    end do
    
  end subroutine GetInitLanduseTransitionRates
    
  !----------------------------------------------------------------------------------------------------

  subroutine FatesGrazing(prt, ft, land_use_label, height)

    use PRTGenericMod,    only : leaf_organ
    use PRTGenericMod,    only : prt_vartypes
    use PRTLossFluxesMod, only : PRTHerbivoryLosses
    use EDParamsMod     , only : landuse_grazing_rate
    use EDParamsMod     , only : landuse_grazing_maxheight
    use EDPftvarcon     , only : EDPftvarcon_inst
    use PRTParametersMod, only : prt_params
    use FatesAllometryMod,only : CrownDepth

    ! apply grazing and browsing to plants as a function of PFT, height (for woody plants), and land use label.

    class(prt_vartypes), intent(inout), pointer :: prt
    integer, intent(in)  :: ft
    integer, intent(in)  :: land_use_label
    real(r8), intent(in) :: height

    real(r8) :: grazing_rate    ! rate of grazing (or browsing) of leaf tissue [day -1]
    real(r8) :: crown_depth
    
    grazing_rate = landuse_grazing_rate(land_use_label) * EDPftvarcon_inst%landuse_grazing_palatability(ft)
    
    if ( grazing_rate .gt. 0._r8) then
       if (prt_params%woody(ft) == itrue) then

          call CrownDepth(height,ft,crown_depth)
          
          grazing_rate = grazing_rate * &
               max(0._r8, min(1._r8, &
               (landuse_grazing_maxheight - (height - crown_depth )) / crown_depth))
       endif

       call PRTHerbivoryLosses(prt, leaf_organ, grazing_rate)
    end if

  end subroutine FatesGrazing

end module FatesLandUseChangeMod
