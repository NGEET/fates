module FatesTestEnvironmentMod
  !
  ! DESCRIPTION:
  ! Prescribed atmospheric and soil boundary conditions for functional tests drivers.
  ! 
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : t_water_freeze_k_1atm
  use LeafBiophysicsMod, only : QSat, GetConstrainedVPress
  
  implicit none 
  private
  
  ! CONSTANTS:
  real(r8), parameter :: sea_level_press      = 101325.0_r8 ! [Pa]
  real(r8), parameter :: gb_well_ventilated   = 2.0e6_r8    ! [umol/m2/s] (~20 s/m equivalent)
  real(r8), parameter :: default_vpd          = 1000.0_r8   ! leaf-to-air VPD [Pa]
  real(r8), parameter :: default_veg_tempk    = 25.0_r8 + t_water_freeze_k_1atm ! [K]
  real(r8), parameter :: default_co2_molfrac  = 380.0e-6_r8 ! [mol/mol] (380 umol/mol)
  real(r8), parameter :: default_o2_molfrac   = 210.0e-3_r8 ! [mol/mol] (210 mmol/mol)
  real(r8), parameter :: default_dayl_factor  = 1.0_r8      ! [0-1] (no seasonal daylength change assumed)
  real(r8), parameter :: default_btran        = 1.0_r8      ! [0-1] (non-limiting water assumed)
  real(r8), parameter :: default_par          = 1500.0_r8   ! [umol/m2/s]
  
  ! shared defaults
  real(r8), public, parameter :: default_nscaler      = 1.0_r8 ! [0-1]
  real(r8), public, parameter :: default_rdark_scaler = 1.0_r8 ! [0-1] (canopy top, no leaf area above)
  
  public :: CanopyVaporPressure
  public :: BtranFromSMP
  public :: SoilMatricPotential
  
  type, public :: environment_type
  
    real(r8) :: tempk          ! vegetation/leaf temperature [K]
    real(r8) :: can_press      ! air pressure at the leaf surface [Pa]
    real(r8) :: can_co2_ppress ! CO2 partial pressure at the leaf surface [Pa]
    real(r8) :: can_o2_ppress  ! O2 partial pressure at the leaf surface [Pa]
    real(r8) :: can_vpress     ! vapor pressure of the canopy air [Pa]
    real(r8) :: gb             ! leaf boundary layer conductance [umol/m2/s]
    real(r8) :: btran          ! soil moisture stress factor [0-1]
    real(r8) :: dayl_factor    ! day-length photosynthetic-capacity acclimation factor [0-1]
    real(r8) :: veg_esat       ! saturation vapor pressure
    real(r8) :: par            ! PAR [umol/m2/s]
  
  contains 
  
    procedure, public :: Init
  
  end type environment_type
  
  contains 
  
  
    ! ==========================================================================
  
    subroutine Init(this)
    !
    ! DESCRIPTION:
    ! Set the prescribed atmospheric and soil boundary conditions

    ! ARGUMENTS:
    class(environment_type), intent(out) :: this ! environment object
    
    ! LOCALS:
    real(r8) :: qs_dummy ! saturation specific humidity output from QSat (unused here)

    ! initialize defaults
    this%tempk = default_veg_tempk ! [K]
    this%can_press = sea_level_press ! [Pa] 
    this%can_co2_ppress = default_co2_molfrac * this%can_press ! [Pa]
    this%can_o2_ppress = default_o2_molfrac * this%can_press ! [Pa]
    this%par = default_par
    this%gb = gb_well_ventilated ! [umol/m2/s]
    this%dayl_factor = default_dayl_factor ! [0-1]
    this%btran = default_btran ! [0-1]
    
    ! the standard reference condition's vapor-pressure state
    call QSat(this%tempk, this%can_press, qs_dummy, this%veg_esat)
    this%can_vpress = CanopyVaporPressure(this%veg_esat)

  end subroutine Init
  
  ! ==========================================================================
  
  function CanopyVaporPressure(veg_esat, vpress_unconstrained, vpd) result(can_vpress)
    !
    ! DESCRIPTION:
    ! Prescribed canopy-air vapor pressure
    !
    ! This value is constrained through GetConstrainedVPress
    !

    ! ARGUMENTS:
    real(r8), intent(in)            :: veg_esat             ! saturation vapor pressure at the current tempk [Pa]
    real(r8), intent(out), optional :: vpress_unconstrained ! uncontrained vapor pressure
    real(r8), intent(in),  optional :: vpd                  ! leaf-to-air VPD to impose [Pa]; default_vpd if absent

    ! RESULT:
    real(r8) :: can_vpress ! canopy air vapor pressure [Pa]

    ! LOCALS:
    real(r8) :: vpd_local ! this call's imposed VPD [Pa]

    vpd_local = default_vpd
    if (present(vpd)) vpd_local = vpd

    can_vpress = veg_esat - vpd_local
    if (present(vpress_unconstrained)) vpress_unconstrained = can_vpress
    can_vpress = GetConstrainedVPress(can_vpress, veg_esat)

  end function CanopyVaporPressure
  
  ! ==========================================================================
  
  pure function BtranFromSMP(smp, smpsc, smpso) result(btran)
    !
    ! DESCRIPTION:
    ! Soil moisture stress factor (btran) from soil matric potential, for a
    ! driver with no soil column.
    !
    ! Production subroutine EDBtranMod.F90::btran_ed builds btran as a rootfrac-weighted
    ! sum of a per-layer term over the whole soil column:
    !
    !   smp_node = max(smpsc, smp_sl(j))
    !   rresis   = min( (eff_porosity_sl(j)/watsat_sl(j)) *
    !                   (smp_node - smpsc)/(smpso - smpsc), 1 )
    !   btran    = SUM_j rootfrac(j)*rresis(j)
    !
    ! evaluated only over layers holding liquid water.
    !
    ! With a single unfrozen layer at root fraction 1, eff_porosity equals watsat 
    ! (effective porosity is saturated water content minus the ice volume, and we assume
    ! no ice here), so that factor = 1.0, so we can remove it here
    !

    ! ARGUMENTS:
    real(r8), intent(in) :: smp   ! soil matric potential [mm], negative
    real(r8), intent(in) :: smpsc ! soil matric potential at full stomatal closure [mm], negative
    real(r8), intent(in) :: smpso ! soil matric potential at full stomatal opening [mm], negative

    ! RESULT:
    real(r8) :: btran ! soil moisture stress factor [0-1]

    ! LOCALS:
    real(r8) :: smp_node ! smp clamped at the full-closure threshold [mm]

    smp_node = max(smpsc, smp)
    btran = min((smp_node - smpsc)/(smpso - smpsc), 1.0_r8)

  end function BtranFromSMP

  ! ==========================================================================

  pure function SoilMatricPotential(soilfrac, smpsc) result(smp)
    !
    ! DESCRIPTION:
    ! Soil matric potential at a given soil water content, expressed as a
    ! fraction of saturation.
    !
    ! This function interpolates linearly between smpsc at zero water
    ! content and that at saturation

    ! ARGUMENTS:
    real(r8), intent(in) :: soilfrac ! soil water content, fraction of saturation [0-1]
    real(r8), intent(in) :: smpsc    ! soil matric potential at full stomatal closure [mm], negative

    ! RESULT:
    real(r8) :: smp ! soil matric potential [mm], negative

    smp = smpsc*(1.0_r8 - soilfrac)

  end function SoilMatricPotential

  ! ==========================================================================
  

end module FatesTestEnvironmentMod