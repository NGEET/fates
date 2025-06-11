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
  real(r8), parameter :: SLA_M2LEAF_PER_KGMOSS = 1.0       ! Specific leaf area (m2/kg)
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
subroutine moss(alff, cla_m2_per_plot, decLit_t_per_haplot, moss_biom_kg_per_plot_inout, moss_litter_flux_t_per_haplot, &
                livemoss_depth_m)
  !
  !  Calculates annual moss growth and mortality
  !  Adapted from Bonan and Korzukhin 1989 Vegetatio 84:31-44
  !  Further adapted from Foster et al. (2019, Ecol. Mod., doi: 10.1016/j.ecolmodel.2019.108765)
  !

  ! Arguments
  real(r8), intent(in)    :: alff    ! Available light on the forest floor (0-1)
  real(r8), intent(in)    :: cla_m2_per_plot     ! Cumulative leaf area on forest floor (m2)
  real(r8), intent(in)    :: decLit_t_per_haplot  ! Fresh deciduous leaf litter (t/ha)
  real(r8), intent(inout) :: moss_biom_kg_per_plot_inout  ! Moss biomass (kg, not kg/m2)
  real(r8), intent(out)   :: moss_litter_flux_t_per_haplot  ! Moss biomass flux to litter (t/ha)
  real(r8), intent(out)   :: livemoss_depth_m  ! Depth (m) of live moss layer

  ! Local variables
  real(r8) :: moss_biom_kg_per_plot_before    ! Moss biomass (kg per plot) before this timestep
  real(r8) :: moss_biom_kg_per_m2plot_before  ! Moss biomass (kg/m2) before this timestep
  real(r8) :: moss_biom_kg_per_m2plot_after   ! Moss biomass (kg/m2) after this timestep
  real(r8) :: al        ! Available light
  real(r8) :: algf      ! Available light growth factor
  real(r8) :: fcgf      ! Forest cover growth factor
  real(r8) :: dlgf      ! Deciduous leaf litter growth factor
  real(r8) :: ddgf      ! Moisture growth factor
  real(r8) :: assim_kg_per_m2leaf     ! Moss assimilation rate (kg/m2)
  real(r8) :: assim_eff_kg_per_kgmoss ! Effective assimilation (kg/kg)
  real(r8) :: assim_eff_kg_per_m2plot ! Assimilation (kg/m2)
  real(r8) :: repro_eff_kg_per_kgmoss ! Effective reproduction (kg/kg)
  real(r8) :: repro_eff_kg_per_m2plot ! Effective reproduction (kg/m2)
  real(r8) :: prod_kg_per_m2plot      ! Moss production (kg/m2)
  real(r8) :: moss_to_litter_flux_kg_per_m2plot  ! Flux from moss to litter (kg/m2)
  real(r8) :: moss_to_litter_flux_kg_per_plot    ! Flux from moss to litter (kg)

  ! Save this for later
  moss_biom_kg_per_plot_before = moss_biom_kg_per_plot_inout

  ! Convert moss biomass in kg to kg/m2
  moss_biom_kg_per_m2plot_before = moss_biom_kg_per_plot_before/plotsize_m2

  ! Light growth multiplier
  al = exp(-1.0*EXT*(cla_m2_per_plot / plotsize_m2 + moss_biom_kg_per_m2plot_before * SLA_M2LEAF_PER_KGMOSS))
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
  if (decLit_t_per_haplot > 0.0) then
      dlgf = exp(-0.2932*decLit_t_per_haplot)
      if (dlgf <= 0.0) dlgf = 0.0
      if (dlgf >= 1.0) dlgf = 1.0
  else
      dlgf = 1.0
  end if

  ! Moisture growth factor
  ! TODO: Implement this
  ddgf = 1.0

  ! Moss assimilation
  assim_kg_per_m2leaf = PMAX_KG_PER_M2_PER_YR*algf*fcgf*dlgf*ddgf

  ! Effective reproduction (fraction of live biomass)
  repro_eff_kg_per_kgmoss = SPORES_KG_PER_KGMOSS*dlgf*ddgf
  repro_eff_kg_per_m2plot = repro_eff_kg_per_kgmoss * moss_biom_kg_per_m2plot_before

  ! Effective assimilation
  assim_eff_kg_per_kgmoss = SLA_M2LEAF_PER_KGMOSS*assim_kg_per_m2leaf*(1.0 - repro_eff_kg_per_kgmoss)
  assim_eff_kg_per_m2plot = assim_eff_kg_per_kgmoss * moss_biom_kg_per_m2plot_before
  prod_kg_per_m2plot = assim_eff_kg_per_m2plot - moss_biom_kg_per_m2plot_before*(Q_KG_PER_KGMOSS + B_KG_PER_KGMOSS) + repro_eff_kg_per_m2plot

  if (moss_biom_kg_per_m2plot_before + prod_kg_per_m2plot < 0.0) then
      ! Not enough moss to account for mortality/respiration
      ! Set moss loss to all of current biomass
      prod_kg_per_m2plot = -1.0*moss_biom_kg_per_m2plot_before
  end if

  ! Calculate litter flux (kg)
  ! TODO: Flux to litter should only come from mortality. Respiration should go to atmosphere.
  moss_to_litter_flux_kg_per_m2plot = assim_eff_kg_per_m2plot + repro_eff_kg_per_m2plot - prod_kg_per_m2plot
  moss_to_litter_flux_kg_per_m2plot = max(0.0, moss_to_litter_flux_kg_per_m2plot)

  ! Update biomass
  moss_biom_kg_per_m2plot_after = moss_biom_kg_per_m2plot_before + prod_kg_per_m2plot

  ! Thickness of live moss layer (m)
  livemoss_depth_m = moss_biom_kg_per_m2plot_after / BULK_MOSS_KG_PER_M3

  ! Convert certain variables to their UVAFME outputs
  ! TODO: Change these to what FATES needs
  moss_to_litter_flux_kg_per_plot = moss_to_litter_flux_kg_per_m2plot * plotsize_m2
  moss_litter_flux_t_per_haplot = moss_to_litter_flux_kg_per_plot/plotsize_m2*HEC_TO_M2*KG_TO_T
  moss_biom_kg_per_plot_inout = moss_biom_kg_per_m2plot_after * plotsize_m2



end subroutine moss

! =====================================================================================


end module FatesMossMod
