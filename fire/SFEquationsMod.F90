module SFEquationsMod
   ! ============================================================================
   ! Helper methods for the FATES SPITFIRE model
   ! Most equations come from:
   !  Thonicke et al. 2010, Biogeosciences 7:1991-2011
   !  Rothermel et al. 1972, Research Paper INT 115
   !  Albini et al. 1976, Research Report INT 30
   ! ============================================================================

   use FatesFuelClassesMod, only : num_fuel_classes, fuel_classes
   use FatesConstantsMod,   only : r8 => fates_r8
   use FatesConstantsMod,   only : nearzero
   use FatesGlobals,        only : endrun => fates_endrun
   use shr_log_mod,         only : errMsg => shr_log_errMsg

   implicit none
   private

   public :: MaximumReactionVelocity
   public :: OptimumReactionVelocity
   public :: OptimumPackingRatio
   public :: MoistureCoefficient
   public :: ReactionIntensity
   public :: HeatofPreignition
   public :: EffectiveHeatingNumber
   public :: HeatSink
   public :: WindFactor
   public :: PropagatingFlux
   public :: ForwardRateOfSpread
   public :: BackwardRateOfSpread
   public :: FireDuration
   public :: LengthToBreadth
   public :: FireSize
   public :: AreaBurnt
   public :: FireIntensity
   public :: ScorchHeight
   public :: CrownFractionBurnt
   public :: BarkThickness
   public :: CambialMortality
   public :: TotalFireMortality
   public :: CrownFireMortality
   public :: CriticalResidenceTime
   public :: cambial_mort

   ! for error message writing
   character(len=*), parameter :: sourcefile = __FILE__

contains

   real(r8) function OptimumPackingRatio(SAV)
      !
      !  DESCRIPTION:
      !  Calculates optimum packing ratio [unitless]
      !  Equation A6 in Thonicke et al. 2010
      !  Rothermel 1972 Eq. 37
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: SAV ! surface area to volume ratio of fuel [/cm]

      ! CONSTANTS:
      real(r8), parameter :: a = 0.20395_r8
      real(r8), parameter :: b = -0.8189_r8

      if (SAV < nearzero) then
         OptimumPackingRatio = 0.0_r8
      else
         OptimumPackingRatio = a*(SAV**b)
      end if

   end function OptimumPackingRatio

   !-------------------------------------------------------------------------------------

   real(r8) function MaximumReactionVelocity(SAV)
      !
      !  DESCRIPTION:
      !  Calculates maximum reaction velocity in /min
      !
      !  From Equation 36 in Rothermel 1972; Fig. 12
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: SAV ! fuel surface area to volume ratio [/cm]

      if (SAV < nearzero) then
         MaximumReactionVelocity = 0.0_r8
      else
         MaximumReactionVelocity = 1.0_r8/(0.0591_r8 + 2.942_r8*(SAV**(-1.5_r8)))
      end if

   end function MaximumReactionVelocity

   !-------------------------------------------------------------------------------------

   real(r8) function OptimumReactionVelocity(max_reaction_vel, SAV, beta_ratio)
      !
      !  DESCRIPTION:
      !  Calculates optimum reaction velocity in /min
      !
      !  Reaction velocity (i.e. rate of fuel consumption) that would exist if the
      !  fuel were free of moisture and contained minerals at the same reaction
      !  concentration as alpha cellulose
      !
      !  From Equation 38 in Rothermel 1972; Fig. 11
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: max_reaction_vel ! maximum reaction velocity [/min]
      real(r8), intent(in) :: SAV              ! fuel surface area to volume ratio [/cm]
      real(r8), intent(in) :: beta_ratio       ! ratio of packing ratio to optimum packing ratio [0-1]

      ! LOCALS:
      real (r8) :: a, a_beta ! intermediate variables

      ! Equations in Table A1 Thonicke et al. 2010
      a = 8.9033_r8*(SAV**(-0.7913_r8))
      a_beta = exp(a*(1.0_r8 - beta_ratio))

      OptimumReactionVelocity = max_reaction_vel*(beta_ratio**a)*a_beta

   end function OptimumReactionVelocity

   !-------------------------------------------------------------------------------------

   real(r8) function MoistureCoefficient(moisture, MEF)
      !
      !  DESCRIPTION:
      !  Calculates the moisture dampening coefficient for reaction intensity
      !  Based on Equation in table A1 Thonicke et al. 2010.
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: moisture ! fuel moisture [m3/m3]
      real(r8), intent(in) :: MEF      ! fuel moisture of extinction [m3/m3]

      ! LOCALS:
      real(r8) :: mw_weight ! relative fuel moisture/fuel moisture of extinction

      if (MEF < nearzero) then
         ! this really should never happen - essentially this means fuel can never burn
         ! but we are putting this here to avoid divide by zeros
         MoistureCoefficient = 0.0_r8
         return
      end if

      ! average values for litter pools (dead leaves, twigs, small and large branches), plus grass
      mw_weight = moisture/MEF

      ! MoistureCoefficient is unitless
      MoistureCoefficient = 1.0_r8 - (2.59_r8*mw_weight) + (5.11_r8*(mw_weight**2.0_r8)) - &
         (3.52_r8*(mw_weight**3.0_r8))

      if (MoistureCoefficient < nearzero) MoistureCoefficient = 0.0_r8

   end function MoistureCoefficient

   !-------------------------------------------------------------------------------------

   real(r8) function ReactionIntensity(fuel_loading, SAV,  beta_ratio, moisture, MEF)
      !
      !  DESCRIPTION:
      !  Calculates reaction intensity in kJ/m2/min
      !
      !  Rate of energy release per unit area within the flaming front
      !

      ! USES
      use SFParamsMod, only : SF_val_fuel_energy, SF_val_miner_damp

      ! ARGUMENTS:
      real(r8), intent(in) :: fuel_loading ! net fuel loading [kg/m2]
      real(r8), intent(in) :: SAV          ! fuel surface area to volume ratio [/cm]
      real(r8), intent(in) :: beta_ratio   ! ratio of packing ratio to optimum packing ratio [0-1]
      real(r8), intent(in) :: moisture     ! fuel moisture [m3/m3]
      real(r8), intent(in) :: MEF          ! fuel moisture of extinction [m3/m3]

      ! LOCALS:
      real(r8) :: max_reaction_vel ! maximum reaction velocity
      real(r8) :: opt_reaction_vel ! optimum reaction velocity
      real(r8) :: moist_coeff      ! moisture dampening coefficient [0-1]

      ! calculate maximum reaction velocity [/min]
      max_reaction_vel = MaximumReactionVelocity(SAV)

      ! calculate optimum reacion velocity [/min]
      opt_reaction_vel = OptimumReactionVelocity(max_reaction_vel, SAV, beta_ratio)

      ! calculate moisture dampening coefficient [0-1]
      moist_coeff = MoistureCoefficient(moisture, MEF)

      ReactionIntensity = opt_reaction_vel*fuel_loading*SF_val_fuel_energy*              &
         moist_coeff*SF_val_miner_damp

   end function ReactionIntensity

   !-------------------------------------------------------------------------------------

   real(r8) function HeatofPreignition(fuel_moisture)
      !
      !  DESCRIPTION:
      !  Calculates heat of pre-ignition in kJ/kg
      !
      !  Heat of pre-ignition is the heat required to bring a unit weight of fuel to
      !  ignition
      !
      !  Equation A4 in Thonicke et al. 2010
      !  Rothermel EQ12 = 250 Btu/lb + 1116 Btu/lb * average_moisture
      !  conversion of Rothermel (1972) EQ12 in BTU/lb to current kJ/kg
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: fuel_moisture ! fuel moisture [m3/m3]

      ! CONTANTS:
      real(r8), parameter :: q_dry = 581.0_r8 ! heat of pre-ignition of dry fuels [kJ/kg]

      HeatofPreignition = q_dry + 2594.0_r8*fuel_moisture

   end function HeatofPreignition

   !-------------------------------------------------------------------------------------

   real(r8) function EffectiveHeatingNumber(SAV)
      !
      !  DESCRIPTION:
      !  Calculates effective heating number [unitless]
      !
      !  Proportion of a fuel particle that is heated to ignition temperature at the time
      !  flaming combustion starts
      !
      !  Equation A3 in Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: SAV ! fuel surface area to volume ratio [/cm]

      if (SAV < nearzero) then
         EffectiveHeatingNumber = 0.0_r8
      else
         EffectiveHeatingNumber = exp(-4.528_r8/SAV)
      end if

   end function EffectiveHeatingNumber

   !-------------------------------------------------------------------------------------
   real(r8) function HeatSink(q_ig, eps, wf, mean_FBD, wf_dead, wf_live)
      !
      ! DESCRIPTION
      ! Calculate the total heat sink that are required to pre-heat the fuel bed
      ! and ignite the fuels per unit volume [kJ m-3]
      ! EQ 78 Rothermel 1972
      !
      ! ARGUMENTS
      real(r8),  intent(in)  :: q_ig(num_fuel_classes)         ! heat of pre-ignition per fuel class [kJ kg-1]
      real(r8),  intent(in)  :: eps(num_fuel_classes)          ! effective heating number per fue class [unitless]
      real(r8),  intent(in)  :: wf(num_fuel_classes)           ! weighting factor for each fuel class given SA:V [unitess]
      real(r8),  intent(in)  :: mean_FBD                       ! weighted fuel bed bulk density [kg m-3]
      real(r8),  intent(in)  :: wf_dead                        ! fractional total surface area contributed by dead fuels [unitless]
      real(r8),  intent(in)  :: wf_live                        ! fractional total surface area contributed by live fuels [unitless]
      !
      ! LOCALS
      real(r8)    :: sum_dead     ! sum of the product of wf, q_ig, and effective heating number across dead fuel classes
      real(r8)    :: sum_live     ! sum of the product of wf, q_ig, and effective heating number across live fuel classes
      integer     :: i            ! looping index
      sum_dead = 0.0_r8
      sum_live = 0.0_r8
      do i = 1, num_fuel_classes
         if(i /= fuel_classes%live_grass())then
            sum_dead = sum_dead + wf(i)*eps(i)*q_ig(i)
         else
            sum_live = sum_live + wf(i)*eps(i)*q_ig(i)
         end if
      end do
      HeatSink = mean_FBD * (wf_dead*sum_dead + wf_live*sum_live)


   end function HeatSink


   !-------------------------------------------------------------------------------------

   real(r8) function WindFactor(wind_speed, beta_ratio, SAV)
      !
      !  DESCRIPTION:
      !  Calculates wind factor for the rate of spread equation [unitless]
      !
      !  Accounts for effect of wind speed increasing ROS
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: wind_speed ! wind speed [m/min]
      real(r8), intent(in) :: beta_ratio ! relative packing ratio [unitless]
      real(r8), intent(in) :: SAV        ! fuel surface area to volume ratio [/cm]

      ! LOCALS:
      real(r8) :: b, c, e ! temporary variables

      ! Equation A7 in Thonicke et al. 2010 per eqn 49 from Rothermel 1972
      b = 0.15988_r8*(SAV**0.54_r8)

      ! Equation A8 in Thonicke et al. 2010 per eqn 48 from Rothermel 1972
      c = 7.47_r8*(exp(-0.8711_r8*(SAV**0.55_r8)))

      ! Equation A9 in Thonicke et al. 2010 (appears to have typo, using coefficient Eq. 50 Rothermel 1972)
      e = 0.715_r8*(exp(-0.01094_r8*SAV))

      ! Equation A5 in Thonicke et al. 2010
      ! convert wind_speed (wind at elev relevant to fire) from m/min to ft/min for Rothermel ROS Eq.
      WindFactor = c*((3.281_r8*wind_speed)**b)*(beta_ratio**(-e))

   end function WindFactor

   !-------------------------------------------------------------------------------------

   real(r8) function PropagatingFlux(beta, SAV)
      !
      !  DESCRIPTION:
      !  Calculates propagating flux ratio [unitless]
      !
      !  Proportion of reaction intensity that heats adjacent fuel particles to ignition
      !
      !  Equation A2 in Thonicke et al. 2010 and Eq. 42 Rothermel 1972
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: beta ! packing ratio [unitless]
      real(r8), intent(in) :: SAV  ! fuel surface area to volume ratio [/cm]

      PropagatingFlux = (exp((0.792_r8 + 3.7597_r8*(SAV**0.5_r8))*(beta + 0.1_r8)))/     &
         (192.0_r8 + 7.9095_r8*SAV)

   end function PropagatingFlux

   !-------------------------------------------------------------------------------------

   real(r8) function ForwardRateOfSpread(heat_sink, i_r, xi, phi_wind)
      !
      !  DESCRIPTION:
      !  Calculates forward rate of spread [m/min]
      !
      !  Flaming front of a surface fire
      !
      !  Equation 9. Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: heat_sink ! total heat required to ignite fuel bed [kJ m-3]
      real(r8), intent(in) :: i_r          ! reaction intensity [kJ/m2/min]
      real(r8), intent(in) :: xi           ! propagating flux [unitless]
      real(r8), intent(in) :: phi_wind     ! wind factor [unitless]

      if (heat_sink <= nearzero) then
         ForwardRateOfSpread = 0.0_r8
      else
         ForwardRateOfSpread = (i_r*xi*(1.0_r8 + phi_wind))/heat_sink
      endif

   end function ForwardRateOfSpread

   !-------------------------------------------------------------------------------------

   real(r8) function BackwardRateOfSpread(ros_front, wind_speed)
      !
      !  DESCRIPTION:
      !  Calculates backwards rate of spread [m/min]
      !
      !  Equation 10 in Thonicke et al. 2010
      !  backward ROS from Can FBP System (1992)
      !  backward ROS wind not changed by vegetation
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: ros_front  ! forward rate of spread [m/min]
      real(r8), intent(in) :: wind_speed ! wind speed [m/min]

      BackwardRateOfSpread = ros_front*exp(-0.012_r8*wind_speed)

   end function BackwardRateOfSpread

   !-------------------------------------------------------------------------------------

   real(r8) function FireDuration(FDI)
      !
      !  DESCRIPTION:
      !  Calculates fire duration [min]
      !
      !  Equation 14 in Thonicke et al. 2010
      !
      use SFParamsMod, only : SF_val_max_durat, SF_val_durat_slope

      ! ARGUMENTS:
      real(r8), intent(in) :: FDI  ! fire danger index [0-1]

      FireDuration = (SF_val_max_durat + 1.0_r8)/(1.0_r8 + SF_val_max_durat*             &
         exp(SF_val_durat_slope*FDI))

   end function FireDuration

   !-------------------------------------------------------------------------------------

   real(r8) function LengthToBreadth(effective_windspeed, tree_fraction)
      !
      !  DESCRIPTION:
      !  Calculates length to breadth ratio of fire ellipse [unitless], used for calculating area burnt
      !
      !     Canadian Forest Fire Behavior Prediction System Ont.Int.Rep. ST-X-3, 1992
      !     Information Report GLC-X-10, Wotten et al. 2009
      !
      use FatesConstantsMod, only : m_per_km, min_per_hr

      ! ARGUMENTS:
      real(r8), intent(in) :: effective_windspeed ! effective windspeed [m/min]
      real(r8), intent(in) :: tree_fraction       ! tree fraction [0-1]

      ! LOCALS:
      real(r8) :: windspeed_km_hr ! effective windspeed, converted to correct units [km/hr]

      ! CONSTANTS:
      real(r8), parameter :: lb_threshold = 0.55_r8 ! tree canopy fraction below which to use grassland length-to-breadth eqn

      windspeed_km_hr = effective_windspeed/m_per_km*min_per_hr

      if (windspeed_km_hr < 1.0_r8) then
         LengthToBreadth = 1.0_r8
      else
         if (tree_fraction > lb_threshold) then
            LengthToBreadth = 1.0_r8 + 8.729_r8*((1.0_r8 - exp(-0.03_r8*windspeed_km_hr))**2.155_r8)
         else
            LengthToBreadth = 1.1_r8*(windspeed_km_hr**0.464_r8)
         endif
      endif

   end function LengthToBreadth

   !-------------------------------------------------------------------------------------

   real(r8) function FireSize(length_to_breadth, ros_back, ros_forward, fire_duration)
      !
      !  DESCRIPTION:
      !  Calculates fire size [m2]
      !
      !  Eq 14 Arora and Boer JGR 2005 (area of an ellipse)
      !
      use FatesConstantsMod, only : pi_const

      ! ARGUMENTS:
      real(r8), intent(in) :: length_to_breadth ! length to breadth ratio of fire ellipse [unitless]
      real(r8), intent(in) :: ros_back          ! backwards rate of spread [m/min]
      real(r8), intent(in) :: ros_forward       ! forward rate of spread [m/min]
      real(r8), intent(in) :: fire_duration     ! fire duration [min]

      ! LOCALS:
      real(r8) :: dist_back    ! distance fire has travelled backwards [m]
      real(r8) :: dist_forward ! distance fire has travelled forward [m]
      real(r8) :: fire_size    ! area of fire [m2]

      dist_back = ros_back*fire_duration
      dist_forward = ros_forward*fire_duration

      ! Eq 14 Arora and Boer JGR 2005 (area of an ellipse)
      if (length_to_breadth < nearzero) then
         FireSize = 0.0_r8
      else
         FireSize = (pi_const/(4.0_r8*length_to_breadth))*((dist_forward + dist_back)**2.0_r8)
      end if

   end function FireSize

   !-------------------------------------------------------------------------------------

   real(r8) function AreaBurnt(fire_size, num_ignitions, FDI)
      !
      !  DESCRIPTION:
      !  Calculates area burnt [m2/km2/day]
      !
      ! daily area burnt = size fires in m2 * num ignitions per day per km2 * prob ignition starts fire
      ! Thonicke 2010 Eq. 1
      !
      ! the denominator in the units of currentSite%NF is total gridcell area, but since we assume that ignitions
      ! are equally probable across patches, currentSite%NF is equivalently per area of a given patch
      ! thus AreaBurnt has units of m2 burned area per km2 patch area per day
      !
      ! TO DO: Connect here with the Li & Levis GDP fire suppression algorithm.
      !     Equation 16 in arora and boer model JGR 2005
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: fire_size     ! fire size [m2]
      real(r8), intent(in) :: num_ignitions ! number of ignitions [/km2/day]
      real(r8), intent(in) :: FDI           ! fire danger index [0-1]

      AreaBurnt = fire_size*num_ignitions*FDI

   end function AreaBurnt

   !-------------------------------------------------------------------------------------

   real(r8) function FireIntensity(fuel_consumed, ros)
      !
      !  DESCRIPTION:
      !  Calculates fire intensity [kW/m]
      !  Eq 15 Thonicke et al 2010

      use SFParamsMod, only : SF_val_fuel_energy

      ! ARGUMENTS:
      real(r8), intent(in) :: fuel_consumed ! fuel consumed [kg/m2]
      real(r8), intent(in) :: ros           ! rate of spread [m/s]

      FireIntensity = SF_val_fuel_energy*fuel_consumed*ros

   end function FireIntensity

   !---------------------------------------------------------------------------------------

   real(r8) function ScorchHeight(alpha_SH, FI)
      !
      !  DESCRIPTION:
      !  Calculates scorch height [m]
      !
      !  Equation 16 in Thonicke et al. 2010
      !  Van Wagner 1973 Eq. 8; Byram (1959)
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: alpha_SH ! alpha parameter for scorch height equation
      real(r8), intent(in) :: FI       ! fire intensity [kW/m]

      if (FI < nearzero) then
         ScorchHeight = 0.0_r8
      else
         ScorchHeight = alpha_SH*(FI**0.667_r8)
      end if

   end function ScorchHeight

   !---------------------------------------------------------------------------------------

   real(r8) function CrownFractionBurnt(SH, height, crown_depth)
      !
      !  DESCRIPTION:
      !  Calculates fraction of the crown burnt of woody plants
      !  Equation 17 in Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: SH          ! scorch height [m]
      real(r8), intent(in) :: height      ! plant height [m]
      real(r8), intent(in) :: crown_depth ! crown depth [m]

      if (crown_depth < nearzero) then
         CrownFractionBurnt = 0.0_r8
      else
         CrownFractionBurnt = (SH - height + crown_depth)/crown_depth
         CrownFractionBurnt = min(1.0_r8, max(0.0_r8, CrownFractionBurnt))
      end if

   end function CrownFractionBurnt

   !---------------------------------------------------------------------------------------

   real(r8) function BarkThickness(bark_scalar, dbh)
      !
      !  DESCRIPTION:
      !  Calculates bark thickness [cm]
      !  Equation 21 in Thonicke et al 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: bark_scalar ! bark per dbh [cm/cm]
      real(r8), intent(in) :: dbh         ! diameter at breast height [cm]

      BarkThickness = bark_scalar*dbh

      if (BarkThickness < nearzero) then
         call endrun(msg="bark thickness is negative", &
            additional_msg=errMsg(sourcefile, __LINE__))
      end if

   end function BarkThickness

   !---------------------------------------------------------------------------------------

   real(r8) function CriticalResidenceTime(bark_thickness)
      !
      !  DESCRIPTION:
      !  Calculates critical fire residence time for cambial damage [min]
      !  Equation 19 in Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: bark_thickness ! bark thickness [cm]

      CriticalResidenceTime = 2.9_r8*bark_thickness**2.0_r8

   end function CriticalResidenceTime

   !---------------------------------------------------------------------------------------

   real(r8) function CambialMortality(bark_scalar, dbh, tau_l)
      !
      !  DESCRIPTION:
      !  Calculates rate of cambial damage mortality [0-1]
      !  Equation 19 in Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: bark_scalar ! cm bark per cm dbh [cm/cm]
      real(r8), intent(in) :: dbh         ! diameter at breast height [cm]
      real(r8), intent(in) :: tau_l       ! residence time of fire [min]

      ! LOCALS:
      real(r8) :: bark_thickness ! bark thickness [cm]
      real(r8) :: tau_c          ! critical fire residence time for cambial damage [min]
      real(r8) :: tau_r          ! relative fire residence time (actual / critical)

      ! calculate bark thickness based of bark scalar parameter and DBH
      bark_thickness = BarkThickness(bark_scalar, dbh)

      ! calculate critical residence time for cambial damage [min]
      tau_c = CriticalResidenceTime(bark_thickness)

      ! relative residence time
      tau_r = tau_l/tau_c

      CambialMortality = cambial_mort(tau_r)

   end function CambialMortality

   !---------------------------------------------------------------------------------------

   real(r8) function cambial_mort(tau_r)
      !
      !  DESCRIPTION:
      !  Helper function for CambialMortality
      !  Calculates rate of cambial damage mortality [0-1]
      !  Equation 19 in Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: tau_r ! relative residence time of fire

      if (tau_r >= 2.0_r8) then
         cambial_mort = 1.0_r8
      else if (tau_r < 2.0_r8 .and. tau_r > 0.22_r8) then
         cambial_mort = 0.563_r8*tau_r - 0.125_r8
      else
         cambial_mort = 0.0_r8
      end if

   end function cambial_mort

   !---------------------------------------------------------------------------------------

   real(r8) function CrownFireMortality(crown_kill, fraction_crown_burned)
      !
      !  DESCRIPTION:
      !  Calculates rate of mortality from crown scorching [0-1]
      !  Equation 19 in Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: crown_kill            ! parameter for crown kill cm bark per cm dbh [cm/cm]
      real(r8), intent(in) :: fraction_crown_burned ! fraction of the crown burned [0-1]

      CrownFireMortality = crown_kill*fraction_crown_burned**3.0_r8
      if (CrownFireMortality > 1.0_r8) CrownFireMortality = 1.0_r8
      if (CrownFireMortality < nearzero) CrownFireMortality = 0.0_r8

   end function CrownFireMortality

   !---------------------------------------------------------------------------------------

   real(r8) function TotalFireMortality(crownfire_mort, cambial_damage_mort)
      !
      !  DESCRIPTION:
      !  Calculates rate of mortality from wildfire [0-1]
      !  Equation 18 in Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: crownfire_mort      ! mortality rate from crown scorching [0-1]
      real(r8), intent(in) :: cambial_damage_mort ! mortality rate from cambial damage [0-1]

      if (crownfire_mort > 1.0_r8 .or. cambial_damage_mort > 1.0_r8) then
         TotalFireMortality = 1.0_r8
      else
         TotalFireMortality = crownfire_mort + cambial_damage_mort - (crownfire_mort*cambial_damage_mort)
      end if

      if (TotalFireMortality > 1.0_r8) TotalFireMortality = 1.0_r8
      if (TotalFireMortality < nearzero) TotalFireMortality = 0.0_r8

   end function TotalFireMortality

   !---------------------------------------------------------------------------------------

end module SFEquationsMod
