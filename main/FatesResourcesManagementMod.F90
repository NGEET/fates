module FatesResourcesManagementMod

  use FatesConstantsMod, only : r8 => fates_r8

  implicit none 
  private

  type, public :: fates_resources_management_type
    
  real(r8) :: trunk_product_site ! actual trunk product at site level [kgC/site]
  real(r8) :: harvest_debt       ! amount of kgC per site that was not successfully harvested [kgC/site]
  real(r8) :: harvest_debt_sec   ! amount of kgC per site from secondary patches that was
                                 !  not successfully harvested [kgC/site]

  ! DEBUG VARIABLES
  real(r8) :: delta_litter_stock  ! kgC/site = kgC/ha
  real(r8) :: delta_biomass_stock ! kgC/site
  real(r8) :: delta_individual
  
  contains 
    
    procedure :: ZeroVals

end type fates_resources_management_type

contains 

  subroutine ZeroVals(this)
    !
    ! DESCRIPTION:
    !     sets specific variables in the object to zero
    !

    ! ARGUMENTS:
    class(fates_resources_management_type), intent(inout) :: this ! resources management object
    
    this%harvest_debt           = 0.0_r8
    this%harvest_debt_sec       = 0.0_r8
    this%trunk_product_site     = 0.0_r8

  end subroutine ZeroVals

end module FatesResourcesManagementMod