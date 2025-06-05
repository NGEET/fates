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

  use FatesGlobals, only      : endrun => fates_endrun
  use FatesGlobals, only      : fates_log


  implicit none


  ! PUBLIC MEMBER FUNCTIONS:
  public :: moss

contains

!------------------------------------------------------------------------------
subroutine moss(alff, cla, decLit, drydays, moss_biom_kg, moss_litter_flux, &
                livemoss_depth)
  !
  !  Calculates annual moss growth and mortality
  !  Adapted from Bonan and Korzukhin 1989 Vegetatio 84:31-44
  !  Further adapted from Foster et al. (2019, Ecol. Mod., doi: 10.1016/j.ecolmodel.2019.108765)
  !

  ! Data dictionary: constants - from Bonan & Korzukhin (1989)
  real, parameter :: Q = 0.12       ! Mortality parameter
  real, parameter :: B = 0.136      ! Mortality parameter
  real, parameter :: EXT = 0.5      ! Light extinction coefficient
  real, parameter :: SLA = 1.0       ! Specific leaf area (m2/kg)
  real, parameter :: SPORES = 0.001 ! Proportion of biomass spent on reproduction
  real, parameter :: PMAX = 0.35    ! Maximum moss production (kg/m2/yr)
  real, parameter :: LRMIN = 0.01   ! Light compensation point
  real, parameter :: LRMAX = 0.05   ! Light compensation point

  ! Data dictionary: calling arguments
  real,             intent(in)    :: alff    ! Available light on the forest floor (0-1)
  real,             intent(in)    :: cla     ! Cumulative leaf area on forest floor (m2)
  real,             intent(in)    :: decLit  ! Fresh deciduous leaf litter (t/ha)
  real,             intent(in)    :: drydays ! Drought index (0-1)
  real,             intent(inout) :: moss_biom_kg  ! Moss biomass (kg, not kg/m2)
  real,             intent(out)   :: moss_litter_flux  ! Moss biomass flux to litter (t/ha)
  real,             intent(out)   :: livemoss_depth  ! Depth (m) of live moss layer

  ! Data dictionary: local variables
  real :: biokg     ! Moss biomass (kg/m2)
  real :: al        ! Available light
  real :: algf      ! Available light growth factor
  real :: fcgf      ! Forest cover growth factor
  real :: dlgf      ! Deciduous leaf litter growth factor
  real :: ddgf      ! Moisture growth factor
  real :: assim     ! Moss assimilation rate (kg/m2)
  real :: assim_eff ! Effective assimilation (kg/kg)
  real :: prod      ! Moss production (kg/m2)
  real :: litter    ! Moss litter (kg)

  ! Convert moss biomass in kg to kg/m2
  biokg = moss_biom_kg/plotsize

  ! Light growth multiplier
  al = exp(-1.0*EXT*(cla/plotsize + biokg*SLA))
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
  if (decLit > 0.0) then
      dlgf = exp(-0.2932*decLit)
      if (dlgf <= 0.0) dlgf = 0.0
      if (dlgf >= 1.0) dlgf = 1.0
  else
      dlgf = 1.0
  end if

  ! Soil moisture growth multiplier
  if (drydays > 0.025) then
      ddgf = 0.0
  else
      ddgf = 1.0
  end if

  ! Moss assimilation
  assim = PMAX*algf*fcgf*dlgf*ddgf

  ! Effective assimilation
  assim_eff = SLA*assim*(1.0 - SPORES*dlgf*ddgf)
  prod = assim_eff*biokg - biokg*Q - biokg*B + SPORES*dlgf*ddgf

  if (biokg + prod < 0.0) then
      ! Not enough moss to account for mortality/respiration
      ! Set moss loss to all of current biomass
      prod = -1.0*biokg
  end if

  moss_biom_kg = (biokg + prod)*plotsize

  ! Calculate litter (kg)
  litter = (assim_eff*biokg + SPORES*dlgf*ddgf - prod)*plotsize
  litter = max(0.0, litter)

  ! Convert to tonnes/ha from kg/plot
  moss_litter_flux = litter/plotsize*HEC_TO_M2*KG_TO_T

  ! Thickness of live moss layer (m)
  livemoss_depth = moss_biom_kg/plotsize/BULK_MOSS


end subroutine moss

! =====================================================================================


end module FatesMossMod
