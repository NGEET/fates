module FatesMossMod

! ==============================================================================================
! This module contains the relevant code for moss.
!
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
!
!  MOSS IS AN EXPERIMENTAL OPTION THAT IS STILL UNDERGOING TESTING.
!
! WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
!
! ==============================================================================================

  use FatesConstantsMod , only: r8 => fates_r8
  use FatesGlobals, only      : endrun => fates_endrun
  use FatesGlobals, only      : fates_log


  implicit none

  ! Constants from Bonan & Korzukhin (1989)
  real(r8), parameter :: Q_KG_PER_KGMOSS = 0.12       ! Respiration? parameter
  real(r8), parameter :: B_KG_PER_KGMOSS = 0.136      ! Mortality? parameter
  real(r8), parameter :: EXT = 0.5      ! Light extinction coefficient
  real(r8), parameter :: SLA_M2_PER_KG = 1.0       ! Specific leaf area (m2/kg)
  real(r8), parameter :: SPORES_KG_PER_KGMOSS = 0.001 ! Proportion of biomass spent on reproduction
  real(r8), parameter :: PMAX_KG_PER_M2_PER_YR = 0.35    ! Maximum moss production (kg/m2/yr)
  real(r8), parameter :: LRMIN = 0.01   ! Light compensation point
  real(r8), parameter :: LRMAX = 0.05   ! Light compensation point

  real(r8), parameter :: BULK_MOSS_KG_PER_M3 = 18       ! Moss bulk density (kg/m3)

  ! TODO: Replace these with FATES versions (or delete)
  real(r8), parameter :: HEC_TO_M2 = 10000.0  ! Convert from hectares to m2
  real(r8), parameter :: KG_TO_T = 0.001      ! Convert from kg to metric tons
  real(r8), parameter :: plotsize_m2 = 500.0     ! Area of plots (m2)


  ! PUBLIC MEMBER FUNCTIONS:
  public :: moss

contains

!------------------------------------------------------------------------------
subroutine moss(alff, cla_m2, decLit_t_per_ha, moss_biom_kg, moss_litter_flux_t_per_ha, &
                livemoss_depth_m)
  !
  !  Calculates annual moss growth and mortality
  !  Adapted from Bonan and Korzukhin 1989 Vegetatio 84:31-44
  !  Further adapted from Foster et al. (2019, Ecol. Mod., doi: 10.1016/j.ecolmodel.2019.108765)
  !

  ! Arguments
  real(r8), intent(in)    :: alff    ! Available light on the forest floor (0-1)
  real(r8), intent(in)    :: cla_m2     ! Cumulative leaf area on forest floor (m2)
  real(r8), intent(in)    :: decLit_t_per_ha  ! Fresh deciduous leaf litter (t/ha)
  real(r8), intent(inout) :: moss_biom_kg  ! Moss biomass (kg, not kg/m2)
  real(r8), intent(out)   :: moss_litter_flux_t_per_ha  ! Moss biomass flux to litter (t/ha)
  real(r8), intent(out)   :: livemoss_depth_m  ! Depth (m) of live moss layer

  ! Local variables
  real(r8) :: moss_biom_kg_per_m2     ! Moss biomass (kg/m2)
  real(r8) :: al        ! Available light
  real(r8) :: algf      ! Available light growth factor
  real(r8) :: fcgf      ! Forest cover growth factor
  real(r8) :: dlgf      ! Deciduous leaf litter growth factor
  real(r8) :: ddgf      ! Moisture growth factor
  real(r8) :: assim_kg_per_m2     ! Moss assimilation rate (kg/m2)
  real(r8) :: assim_eff_kg_per_kgmoss ! Effective assimilation (kg/kg)
  real(r8) :: repro_eff_kg_per_kgmoss ! Effective reproduction (kg/kg)
  real(r8) :: prod_kg_per_m2      ! Moss production (kg/m2)
  real(r8) :: litter_kg    ! Moss litter (kg)

  ! Convert moss biomass in kg to kg/m2
  moss_biom_kg_per_m2 = moss_biom_kg/plotsize_m2

  ! Light growth multiplier
  al = exp(-1.0*EXT*(cla_m2 / plotsize_m2 + moss_biom_kg_per_m2 * SLA_M2_PER_KG))
  algf = (al - LRMIN)/(LRMAX - LRMIN)
  algf = max(0.0, algf)
  algf = min(1.0, algf)

  ! Forest cover growth multiplier (alff > 0.75)
  if (alff > 0.75) then
     fcgf = 1.5625 - alff**2
  else
      fcgf = 1.0
  end if
  if (fcgf > 1.0) fcgf = 1.0
  if (fcgf < 0.0) fcgf = 0.0

  ! Deciduous leaf litter growth multiplier
  if (decLit_t_per_ha > 0.0) then
      dlgf = exp(-0.2932*decLit_t_per_ha)
      if (dlgf <= 0.0) dlgf = 0.0
      if (dlgf >= 1.0) dlgf = 1.0
  else
      dlgf = 1.0
  end if

  ! Moisture growth factor
  ! TODO: Implement this
  ddgf = 1.0

  ! Moss assimilation
  assim_kg_per_m2 = PMAX_KG_PER_M2_PER_YR*algf*fcgf*dlgf*ddgf

  ! Effective reproduction (fraction of live biomass)
  repro_eff_kg_per_kgmoss = SPORES_KG_PER_KGMOSS*dlgf*ddgf

  ! Effective assimilation
  assim_eff_kg_per_kgmoss = SLA_M2_PER_KG*assim_kg_per_m2*(1.0 - repro_eff_kg_per_kgmoss)
  prod_kg_per_m2 = assim_eff_kg_per_kgmoss*moss_biom_kg_per_m2 - moss_biom_kg_per_m2*Q_KG_PER_KGMOSS - moss_biom_kg_per_m2*B_KG_PER_KGMOSS + repro_eff_kg_per_kgmoss

  if (moss_biom_kg_per_m2 + prod_kg_per_m2 < 0.0) then
      ! Not enough moss to account for mortality/respiration
      ! Set moss loss to all of current biomass
      prod_kg_per_m2 = -1.0*moss_biom_kg_per_m2
  end if

  moss_biom_kg = (moss_biom_kg_per_m2 + prod_kg_per_m2)*plotsize_m2

  ! Calculate litter (kg)
  litter_kg = (assim_eff_kg_per_kgmoss*moss_biom_kg_per_m2 + repro_eff_kg_per_kgmoss - prod_kg_per_m2)*plotsize_m2
  litter_kg = max(0.0, litter_kg)

  ! Convert to tonnes/ha from kg/plot
  moss_litter_flux_t_per_ha = litter_kg/plotsize_m2*HEC_TO_M2*KG_TO_T

  ! Thickness of live moss layer (m)
  livemoss_depth_m = moss_biom_kg/plotsize_m2/BULK_MOSS_KG_PER_M3


end subroutine moss

! =====================================================================================


end module FatesMossMod
