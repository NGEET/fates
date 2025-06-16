module FatesRxFireMod


  ! ============================================================================
  ! Methods to help with prescribed fire
  ! ============================================================================

  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nearzero

  implicit none
  private
  
  public :: is_prescribed_burn
  public :: is_wildfire
  
  logical function is_prescribed_burn(wildfire_FI, wildfire_ignitions, rx_min_FI,         &
    rx_max_FI, wildfire_FI_thresh)
    !
    !  DESCRIPTION:
    !  Determines if a prescribed burn is happening
    !
    
    ! ARGUMENTS:
    real(r8), intent(in) :: wildfire_FI        ! wildfire fire intensity [kW/m]
    real(r8), intent(in) :: wildfire_ignitions ! wildfire ignitions [count/km2/day]
    real(r8), intent(in) :: rx_min_FI          ! minimum fire energy of prescribed fire [kW/m]
    real(r8), intent(in) :: rx_max_FI          ! maximum fire energy of prescribed fire [kW/m]
    real(r8), intent(in) :: wildfire_FI_thresh ! threshold for fires that spread or go out [kW/m]
    
    ! LOCALS:
    logical :: rx_man             ! prescribed fire using human ignitions
    logical :: rx_hyb             ! prescribed fire due to both lightning strike and human ignitions
    logical :: within_rx_FI_range ! fire intensity is within prescribed burn limits

    ! check if fire intensity falls within prescribed burn range
    within_rx_FI_range = wildfire_FI > rx_min_FI .and. wildfire_FI < rx_max_FI
  
    ! condition for prescribed burn solely due to human ignitions
    rx_man = within_rx_FI_range .and. wildfire_ignitions < nearzero

    ! condition for hybrid prescribed burn (low-intensity fire + human ignitions)
    rx_hyb = within_rx_FI_range .and. wildfire_FI < wildfire_FI_thresh .and. &
      wildfire_ignitions > nearzero

    is_prescribed_burn = rx_man .or. rx_hyb
  
  end logical function is_prescribed_burn
  
  !---------------------------------------------------------------------------------------
  
  logical function fire_has_ignitions_and_intensity(wildfire_FI, wildfire_ignitions, rxfire_maxFI,         &
    wildfire_intensity_thresh)
    !
    !  DESCRIPTION:
    !  Determines if a wildfire is happening
    !
    
    ! ARGUMENTS:
    real(r8), intent(in) :: wildfire_FI        ! wildfire fire intensity [kW/m]
    real(r8), intent(in) :: wildfire_ignitions ! wildfire ignitions [count/km2/day]
    real(r8), intent(in) :: rx_max_FI          ! maximum fire energy of prescribed fire [kW/m]
    real(r8), intent(in) :: wildfire_FI_thresh ! threshold for fires that spread or go out [kW/m]
    
    ! LOCALS:
    logical :: managed_wildfire      ! is it a wildfire with FI lower than the max rxfire intensity? [can either be Rx fire or wildfire]
    logical :: true_wildfire         ! is it a wildfire that cannot be managed?
    logical :: has_ignitions         ! any natural ignitions at the site?
    logical :: above_wildfire_thresh ! above the wildfire energy threshold
    
    has_ignitions = wildfire_ignitions > nearzero
    above_wildfire_thresh = wildfire_FI > wildfire_FI_thresh
  
    managed_wildfire = has_ignitions .and. above_wildfire_thresh .and.  &
      wildfire_FI < rx_max_FI

    true_wildfire = has_ignitions .and. above_wildfire_thresh .and.  &
      wildfire_FI > rx_max_FI
      
    is_wildfire = managed_wildfire .or. true_wildfire
  
  end logical function fire_has_ignitions_and_intensity
  
end module FatesRxFireMod