module FatesLandUseChangeMod

  ! Controls the transfer and initialization of patch structure to land use types

  use FatesGlobals         , only : fates_log
  use FatesConstantsMod    , only : primarylands, secondarylands, pasture_rangelands, crops
  use FatesConstantsMod    , only : n_landuse_cats
  use FatesGlobals         , only : endrun => fates_endrun
  use FatesConstantsMod    , only : r8 => fates_r8
  use FatesConstantsMod    , only : itrue, ifalse
  use FatesInterfaceTypesMod    , only : bc_in_type
  use EDTypesMod           , only : area_site => area

  ! CIME globals
  use shr_infnan_mod       , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod          , only : errMsg => shr_log_errMsg
  
  !
  implicit none
  private
  !
  public :: get_landuse_transition_rates

  ! 03/10/2023 Created By Charlie Koven
  ! ============================================================================

contains

  ! ============================================================================
  subroutine get_landuse_transition_rates(site_in, bc_in, landuse_transition_matrix)


    ! The purpose of this routine is to ingest the land use transition rate information that the host model has read in from a dataset,
    ! aggregate land use types to those being used in the simulation, and output a transition matrix that can be used to drive patch
    ! disturbance rates.

    ! !ARGUMENTS:
    type(ed_site_type) , intent(in) :: site_in
    type(bc_in_type) , intent(in) :: bc_in


  end subroutine get_landuse_transition_rates


end module FatesLandUseChangeMod
