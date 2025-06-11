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
  real(r8), parameter :: FRAC_ASSIM_TO_SPORES = 0.001 ! Proportion of assimilation spent on reproduction. UNUSED but keeping for possible future work.
  real(r8), parameter :: PMAX_KG_PER_M2_PER_YR = 0.35    ! Maximum moss production (kg/m2/yr)
  real(r8), parameter :: LRMIN = 0.01   ! Light compensation point
  real(r8), parameter :: LRMAX = 0.05   ! Light compensation point

  ! Constants from UVAFME (maybe also from Bonan & Korzukhin, 1989, but unsure)
  real(r8), parameter :: BULK_MOSS_KG_PER_M3 = 18       ! Moss bulk density (kg/m3)
  real(r8), parameter :: FCGF_ALFF_THRESH = 0.75  ! Above this level of light on forest floor, moss assimilation starts to decrease
  real(r8), parameter :: FCGF_INTERCEPT = 1.5625  ! Intercept parameter in fcgf equation
  real(r8), parameter :: DLGF_DECLIT_THRESH = 0  ! Above this level of deciduous litter (t/ha plot), moss assimilation starts to decrease
  real(r8), parameter :: DLGF_DECAY = 0.2932  ! Decay parameter of moss assimilation with increasing deciduous leaf litter

  ! TODO: Replace these with FATES versions (or delete)
  real(r8), parameter :: HEC_TO_M2 = 10000.0  ! Convert from hectares to m2
  real(r8), parameter :: KG_TO_T = 0.001      ! Convert from kg to metric tons
  real(r8), parameter :: plotsize_m2 = 500.0     ! Area of plots (m2)


  ! PUBLIC MEMBER FUNCTIONS:
  public :: light_growth_multiplier
  public :: forest_cover_growth_multiplier
  public :: litter_growth_multiplier
  public :: moss

contains

!------------------------------------------------------------------------------
function available_light_under_canopy_and_moss(cla_m2_per_plot, moss_biom_kg_per_m2plot_before) result(al)
  real(r8), intent(in)    :: cla_m2_per_plot  ! Cumulative leaf area on forest floor (m2)
  real(r8), intent(in)    :: moss_biom_kg_per_m2plot_before  ! Moss biomass (kg/m2) before this timestep
  real(r8) :: canopy_lai  ! Leaf area index of canopy (i.e., excluding moss) (m2 leaves / m2 plot)
  real(r8) :: moss_lai  ! Leaf area index of moss (m2 leaves / m2 plot)
  real(r8) :: lai       ! Total leaf area index, canopy + moss (m2 leaves / m2 plot)
  real(r8) :: al        ! Available light

  canopy_lai = cla_m2_per_plot / plotsize_m2
  moss_lai = moss_biom_kg_per_m2plot_before * SLA_M2LEAF_PER_KGMOSS
  lai = canopy_lai + moss_lai
  al = exp(-1.0*EXT*lai)
end function available_light_under_canopy_and_moss

!------------------------------------------------------------------------------
function light_growth_multiplier(cla_m2_per_plot, moss_biom_kg_per_m2plot_before) result(algf)
  real(r8), intent(in)    :: cla_m2_per_plot  ! Cumulative leaf area on forest floor (m2)
  real(r8), intent(in)    :: moss_biom_kg_per_m2plot_before  ! Moss biomass (kg/m2) before this timestep
  real(r8) :: al        ! Available light
  real(r8) :: algf      ! Light growth multiplier

  al = available_light_under_canopy_and_moss(cla_m2_per_plot, moss_biom_kg_per_m2plot_before)
  algf = (al - LRMIN)/(LRMAX - LRMIN)
  algf = max(0.0, algf)
  algf = min(1.0, algf)
end function light_growth_multiplier

!------------------------------------------------------------------------------
function forest_cover_growth_multiplier(alff) result(fcgf)
  real(r8), intent(in)    :: alff    ! Available light on the forest floor (0-1)
  real(r8) :: fcgf      ! Forest cover growth factor
  if (alff > FCGF_ALFF_THRESH) then
      fcgf = FCGF_INTERCEPT - alff**2
  else
      fcgf = 1.0
  end if
  if (fcgf > 1.0) fcgf = 1.0
  if (fcgf < 0.0) fcgf = 0.0
end function forest_cover_growth_multiplier

!------------------------------------------------------------------------------
function litter_growth_multiplier(decLit_t_per_haplot) result(dlgf)
  real(r8), intent(in)    :: decLit_t_per_haplot    ! Fresh deciduous leaf litter (t/ha)
  real(r8) :: dlgf      ! Deciduous leaf litter growth multiplier
  if (decLit_t_per_haplot > DLGF_DECLIT_THRESH) then
      dlgf = exp(-DLGF_DECAY * decLit_t_per_haplot)
      if (dlgf <= 0.0) dlgf = 0.0
      if (dlgf >= 1.0) dlgf = 1.0
  else
      dlgf = 1.0
  end if
end function litter_growth_multiplier

subroutine moss_biomass_change_kg_per_m2(q_kg_per_kg_moss_in, b_kg_per_kg_moss_in, assim_eff, moss_biom_before, moss_resp, moss_mort, moss_biom_after)
  ! ALL UNITS KG / M2 PLOT UNLESS OTHERWISE SPECIFIED
  real(r8), intent(in) :: q_kg_per_kg_moss_in  ! Respiration? parameter
  real(r8), intent(in) :: b_kg_per_kg_moss_in  ! Mortality? parameter
  real(r8), intent(in) :: assim_eff ! Assimilation (kg/m2)
  real(r8), intent(in) :: moss_biom_before   ! Moss biomass (kg/m2) before this timestep
  real(r8), intent(out) :: moss_resp  ! Moss respiration (kg/m2) during this timestep
  real(r8), intent(out) :: moss_mort  ! Moss mortality (kg/m2) during this timestep
  real(r8), intent(out) :: moss_biom_after   ! Moss biomass (kg/m2) after this timestep

  ! Local variables
  real(r8) :: moss_respmort  ! Sum of moss respiration and mortality (kg/m2) during this timestep
  real(r8) :: qb_sum  ! Sum of Q and B
  real(r8) :: moss_biom_change     ! Change in moss biomass (kg/m2) during this timestep

  ! Total losses to respiration and mortality
  moss_resp = q_kg_per_kg_moss_in * moss_biom_before
  moss_mort = b_kg_per_kg_moss_in * moss_biom_before
  moss_respmort = moss_resp + moss_mort

  ! Not enough moss to account for mortality/respiration
  ! Set moss loss to all of current biomass
  if (moss_respmort > moss_biom_before) then
      qb_sum = q_kg_per_kg_moss_in + b_kg_per_kg_moss_in
      moss_resp = moss_biom_before * (q_kg_per_kg_moss_in / qb_sum)
      moss_mort = moss_biom_before * (b_kg_per_kg_moss_in / qb_sum)
      moss_respmort = moss_resp + moss_mort
  end if

  ! Net change in moss biomass
  moss_biom_change = assim_eff - moss_respmort

  ! Update biomass
  moss_biom_after = moss_biom_before + moss_biom_change
end subroutine moss_biomass_change_kg_per_m2

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
  real(r8) :: algf      ! Available light growth factor
  real(r8) :: fcgf      ! Forest cover growth factor
  real(r8) :: dlgf      ! Deciduous leaf litter growth factor
  real(r8) :: ddgf      ! Moisture growth factor
  real(r8) :: assim_kg_per_m2leaf     ! Moss assimilation rate (kg/m2)
  real(r8) :: assim_eff_kg_per_kgmoss ! Effective assimilation (kg/kg)
  real(r8) :: assim_eff_kg_per_m2plot ! Assimilation (kg/m2)
  real(r8) :: moss_biom_change_kg_per_m2plot     ! Change in moss biomass (kg/m2) during this timestep
  real(r8) :: moss_resp  ! Moss respiration (kg/m2) during this timestep
  real(r8) :: moss_mort  ! Moss mortality (kg/m2) during this timestep

  real(r8) :: moss_to_litter_flux_kg_per_m2plot  ! Flux from moss to litter (kg/m2)
  real(r8) :: moss_to_litter_flux_kg_per_plot    ! Flux from moss to litter (kg)
  real(r8) :: moss_respmort_kg_per_kgmoss  ! Moss respiration and mortality (kg / kg moss) during this timestep
  real(r8) :: moss_respmort_kg_per_m2plot  ! Moss respiration and mortality (kg/m2) during this timestep

  ! Save this for later
  moss_biom_kg_per_plot_before = moss_biom_kg_per_plot_inout

  ! Convert moss biomass in kg to kg/m2
  moss_biom_kg_per_m2plot_before = moss_biom_kg_per_plot_before/plotsize_m2

  ! Light growth multiplier
  algf = light_growth_multiplier(cla_m2_per_plot, moss_biom_kg_per_m2plot_before)

  ! Forest cover growth multiplier
  fcgf = forest_cover_growth_multiplier(alff)

  ! Deciduous leaf litter growth multiplier
  dlgf = litter_growth_multiplier(decLit_t_per_haplot)

  ! Moisture growth factor
  ! TODO: Implement this
  ddgf = 1.0

  ! Moss assimilation
  assim_kg_per_m2leaf = PMAX_KG_PER_M2_PER_YR*algf*fcgf*dlgf*ddgf
  assim_eff_kg_per_kgmoss = SLA_M2LEAF_PER_KGMOSS * assim_kg_per_m2leaf
  assim_eff_kg_per_m2plot = (assim_eff_kg_per_kgmoss * moss_biom_kg_per_m2plot_before)

  ! Get fluxes from moss and change in moss biomass
  call moss_biomass_change_kg_per_m2(Q_KG_PER_KGMOSS, B_KG_PER_KGMOSS, assim_eff_kg_per_m2plot, moss_biom_kg_per_m2plot_before, moss_resp, moss_mort, moss_biom_kg_per_m2plot_after)

  ! Get flux from live moss to litter
  moss_to_litter_flux_kg_per_m2plot = moss_resp + moss_mort

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
