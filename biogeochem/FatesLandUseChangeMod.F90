module FatesLandUseChangeMod

  ! Controls the transfer and initialization of patch structure to land use types

  use FatesGlobals         , only : fates_log
  use FatesConstantsMod    , only : primarylands, secondarylands, pasture, rangelands, crops
  use FatesConstantsMod    , only : n_landuse_cats
  use FatesGlobals         , only : endrun => fates_endrun
  use FatesConstantsMod    , only : r8 => fates_r8
  use FatesConstantsMod    , only : itrue, ifalse
  use FatesInterfaceTypesMod    , only : bc_in_type
  use FatesInterfaceTypesMod    , only : hlm_num_luh2_transitions
  use EDTypesMod           , only : area_site => area

  ! CIME globals
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  
  !
  implicit none
  private
  !
  public :: get_landuse_transition_rates
  public :: init_luh2_fates_mapping

  ! module data
  integer :: max_luh2_types_per_fates_lu_type = 5
  CHARACTER(len = 5), protected, DIMENSION(n_landuse_cats,max_luh2_types_per_fates_lu_type) :: luh2_fates_luype_map

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
    real(r8), intent(inout) :: landuse_transition_matrix(n_landuse_cats, n_landuse_cats)

    ! !LOCAL VARIABLES:
    integer :: i_donor, i_receiver, i_luh2_transitions
    character(5) :: donor_name, receiver_name
    character(14) :: transition_name

    ! zero the transition matrix
    landuse_transition_matrix(:,:) = 0._r8

    ! loop over FATES donor and receiver land use types
    do i_donor = 1,n_landuse_cats
       do i_receiver = 1,n_landuse_cats

          ! ignore diagonals of transition matrix
          if ( i_donor .ne. i_receiver ) then

             ! ignore special case of primary -> secondary, which is handled by harvest mechanism
             if ( .not. ((i_donor .eq. primarylands) .and. (i_receiver .eq. secondarylands)) ) then

                do i_luh2_transitions = 1, hlm_num_luh2_transitions

                   ! transition names are written in form xxxxx_to_yyyyy where x and y are donor and receiver state names
                   transition_name = bc_in%hlm_luh_transition_names(i_luh2_transitions)
                   donor_name = transition_name(1:5)
                   receiver_name = transition_name(10:14)

                   if (any(luh2_fates_luype_map(:,i_donor) == donor_name) .and. &
                        any(luh2_fates_luype_map(:,i_receiver) == receiver_name)) then

                      landuse_transition_matrix(i_donor,i_receiver) = &
                           landuse_transition_matrix(i_donor,i_receiver) +  bc_in%hlm_luh_transitions(i_luh2_transitions)

                   end if
                end do
             end if
          end if
       end do
    end do

  end subroutine get_landuse_transition_rates

  subroutine init_luh2_fates_mapping

    ! initialize the character mapping of the LUH2 : FATES correspondance
    luh2_fates_luype_map(:,:) = ''
    
    luh2_fates_luype_map(1,primarylands) = 'primf'
    luh2_fates_luype_map(2,primarylands) = 'primn'

    luh2_fates_luype_map(1,secondarylands) = 'secdf'
    luh2_fates_luype_map(2,secondarylands) = 'secdn'

    luh2_fates_luype_map(1,crops) = 'c3ann'
    luh2_fates_luype_map(2,crops) = 'c4ann'
    luh2_fates_luype_map(3,crops) = 'c3per'
    luh2_fates_luype_map(4,crops) = 'c4per'
    luh2_fates_luype_map(5,crops) = 'c3nfx'

    luh2_fates_luype_map(1,pasture) = 'pastr'

    luh2_fates_luype_map(1,rangelands) = 'range'
    
  end subroutine init_luh2_fates_mapping

end module FatesLandUseChangeMod
