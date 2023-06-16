module FatesMassBalTypeMod

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none 
  private 

  type, public :: site_massbal_type

  ! --------------------------------------------------------------------------------------
  ! This type is used for accounting purposes to ensure that we are not
  ! loosing or creating mass. This type is supposed to be allocated for each element 
  ! we simulate (e.g. carbon12_element, etc)
  ! Note that the unit of "site", is nominally equivalent to 1 hectare
  !
  ! This set of mass checks are for INCREMENTAL checks during the dynamics step.
  ! --------------------------------------------------------------------------------------

  real(r8) :: old_stock ! remember biomass stock from last time  [Kg/site]
  real(r8) :: err_fates ! total mass balance error for FATES processes [kg/site]

  ! --------------------------------------------------------------------------------------
  ! Components of the total site level mass fluxes
  ! --------------------------------------------------------------------------------------

  real(r8) :: gpp_acc          ! accumulated gross primary productivity [kg/site/day]
  real(r8) :: aresp_acc        ! accumulated autotrophic respiration [kg/site/day]
  real(r8) :: net_root_uptake  ! net uptake of carbon or nutrients through the roots [kg/site/day]
                               !   could include exudation, and for N this also includes symbiotic
                               !   fixation
  real(r8) :: seed_in          ! total mass of external seed rain into fates site [kg/site/day]
                               !   this is from external grid-cells or from user parameterization
                               !   (user param seed rain, or dispersal model)
  real(r8) :: seed_out         ! total mass of seeds exported outside of fates site [kg/site/day]
                               !   (this is not used currently, placeholder, rgk feb-2019)
  real(r8) :: frag_out         ! litter and coarse woody debris fragmentation flux [kg/site/day]
  real(r8) :: wood_product     ! total mass exported as wood product [kg/site/day]
  real(r8) :: burn_flux_to_atm ! total mass burned and exported to the atmosphere [kg/site/day]

  real(r8) :: flux_generic_in  ! used for prescribed or artificial input fluxes
                               !  and initialization [kg/site/day]
  real(r8) :: flux_generic_out ! used for prescribed or artificial output fluxes
                               !  for instance when prescribed physiology is on
  real(r8) :: patch_resize_err ! amount of mass gained (or loss when negative)
                               !  due to re-sizing patches when area math starts to lose
                               !  precision
  contains

  procedure :: ZeroMassBalState
  procedure :: ZeroMassBalFlux

  end type site_massbal_type

  contains
       
    subroutine ZeroMassBalState(this)
      !
      !  DESCRIPTION:
      !  Zero state values
      !

      ! ARGUMENTS:
      class(site_massbal_type) :: this ! site mass balance object
      
      this%old_stock = 0._r8
      this%err_fates = 0._r8

    end subroutine ZeroMassBalState

    ! ====================================================================================
    
    subroutine ZeroMassBalFlux(this)
      !
      !  DESCRIPTION:
      !  Zero flux values
      !

      ! ARGUMENTS:
      class(site_massbal_type) :: this ! site mass balance object

      this%gpp_acc          = 0._r8
      this%aresp_acc        = 0._r8
      this%net_root_uptake  = 0._r8
      this%seed_in          = 0._r8
      this%seed_out         = 0._r8
      this%frag_out         = 0._r8
      this%wood_product     = 0._r8
      this%burn_flux_to_atm = 0._r8
      this%flux_generic_in  = 0._r8
      this%flux_generic_out = 0._r8
      this%patch_resize_err = 0._r8

  end subroutine ZeroMassBalFlux

  ! ======================================================================================

end module FatesMassBalTypeMod