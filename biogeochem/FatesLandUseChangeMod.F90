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
  use EDTypesMod           , only : area_site => area

  ! CIME globals
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  
  !
  implicit none
  private

  character(len=*), parameter :: sourcefile = __FILE__

  public :: get_landuse_transition_rates
  public :: get_landusechange_rules
  public :: get_luh_statedata


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
  subroutine get_landuse_transition_rates(bc_in, landuse_transition_matrix)


    ! The purpose of this routine is to ingest the land use transition rate information that the host model has read in from a dataset,
    ! aggregate land use types to those being used in the simulation, and output a transition matrix that can be used to drive patch
    ! disturbance rates.

    ! !ARGUMENTS:
    type(bc_in_type) , intent(in) :: bc_in
    real(r8), intent(inout) :: landuse_transition_matrix(n_landuse_cats, n_landuse_cats)  ! [m2/m2/day]

    ! !LOCAL VARIABLES:
    type(luh2_fates_lutype_map) :: lumap
    integer :: i_donor, i_receiver, i_luh2_transitions, i_luh2_states
    character(5) :: donor_name, receiver_name
    character(14) :: transition_name
    real(r8) :: urban_fraction
    real(r8) :: temp_vector(hlm_num_luh2_transitions)
    logical  :: modified_flag

    ! zero the transition matrix and the urban fraction
    landuse_transition_matrix(:,:) = 0._r8
    urban_fraction = 0._r8

    use_luh_if: if ( hlm_use_luh .eq. itrue ) then

       ! Check the LUH data incoming to see if any of the transitions are NaN
       temp_vector = bc_in%hlm_luh_transitions
       call CheckLUHData(temp_vector,modified_flag)
       if (.not. modified_flag) then
          ! identify urban fraction so that it can be factored into the land use state output
          urban_fraction = bc_in%hlm_luh_states(findloc(bc_in%hlm_luh_state_names,'urban',dim=1))
       end if

       !!TODO: may need some logic here to ask whether or not ot perform land use change on this timestep. current code occurs every day.
       !!If not doing transition every day, need to update units.

       transitions_loop: do i_luh2_transitions = 1, hlm_num_luh2_transitions

          ! transition names are written in form xxxxx_to_yyyyy where x and y are donor and receiver state names
          transition_name = bc_in%hlm_luh_transition_names(i_luh2_transitions)
          donor_name = transition_name(1:5)
          receiver_name = transition_name(10:14)

          ! Get the fates land use type index associated with the luh2 state types
          i_donor= lumap%GetIndex(donor_name)
          i_receiver = lumap%GetIndex(receiver_name)

          ! Avoid transitions with 'urban' as those are handled seperately
          if (.not.(i_donor .eq. fates_unset_int .or. i_receiver .eq. fates_unset_int)) then
             landuse_transition_matrix(i_donor,i_receiver) = &
                  landuse_transition_matrix(i_donor,i_receiver) +  temp_vector(i_luh2_transitions) * years_per_day / (1._r8 - urban_fraction)

          end if
       end do transitions_loop
    end if use_luh_if
  end subroutine get_landuse_transition_rates

  !----------------------------------------------------------------------------------------------------

  function GetLUCategoryFromStateName(this, state_name) result(landuse_category)

    class(luh2_fates_lutype_map) :: this
    character(len=5), intent(in) :: state_name
    integer :: landuse_category

    landuse_category = this%landuse_categories(findloc(this%state_names,state_name,dim=1))

  end function GetLUCategoryFromStateName

  !----------------------------------------------------------------------------------------------------

  subroutine get_landusechange_rules(clearing_matrix)

    ! the purpose of this is to define a ruleset for when to clear the vegetation in transitioning from one land use type to another

    logical, intent(out) :: clearing_matrix(n_landuse_cats,n_landuse_cats)
    integer, parameter    :: ruleset = 1   ! ruleset to apply from table 1 of Ma et al (2020) https://doi.org/10.5194/gmd-13-3203-2020

    ! clearing matrix applies from the donor to the receiver land use type of the newly-transferred patch area
    ! values of clearing matrix: false => do not clear; true => clear

    clearing_matrix(:,:) = .false.

    select case(ruleset)

    case(1)

       clearing_matrix(:,cropland) = .true.
       clearing_matrix(:,pastureland) = .true.
       clearing_matrix(pastureland,rangeland) = .true.
       clearing_matrix(cropland,rangeland) = .true.

    case(2)

       clearing_matrix(:,cropland) = .true.
       clearing_matrix(rangeland,pastureland) = .true.
       clearing_matrix(cropland,pastureland) = .true.
       clearing_matrix(pastureland,rangeland) = .true.
       clearing_matrix(cropland,rangeland) = .true.

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

  end subroutine get_landusechange_rules

  !----------------------------------------------------------------------------------------------------

  subroutine get_luh_statedata(bc_in, state_vector)

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

    ! Check to see if the incoming state vector is NaN.
    temp_vector = bc_in%hlm_luh_states
    call CheckLUHData(temp_vector,modified_flag)
    if (.not. modified_flag) then
       ! identify urban fraction so that it can be factored into the land use state output
       urban_fraction = bc_in%hlm_luh_states(findloc(bc_in%hlm_luh_state_names,'urban',dim=1))
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

    ! check to ensure total area == 1, and correct if not
    if ( abs(sum(state_vector(:)) - 1._r8) .gt. nearzero ) then
       write(fates_log(),*) 'warning: sum(state_vector) = ', sum(state_vector(:))
       state_vector = state_vector / sum(state_vector)
    end if

  end subroutine get_luh_statedata

  !----------------------------------------------------------------------------------------------------

  subroutine CheckLUHData(luh_vector,modified_flag)

    use shr_infnan_mod   , only : isnan => shr_infnan_isnan

    real(r8), intent(inout) :: luh_vector(:)  ! [m2/m2]
    logical, intent(out)    :: modified_flag

    ! Check to see if the incoming luh2 vector is NaN.
    ! This suggests that there is a discepency where the HLM and LUH2 states
    ! there is vegetated ground. E.g. LUH2 data is missing for glacier-margin regions such as Antarctica.
    ! In this case, states should be Nan.  If so,
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
       write(fates_log(),*) 'WARNING: land use state is all NaN; setting state as all primary forest.'
    else if (any(isnan(luh_vector))) then
       if (any(.not. isnan(luh_vector))) then
          write(fates_log(),*) 'ERROR: land use vector has NaN'
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end if

  end subroutine CheckLUHData

end module FatesLandUseChangeMod
