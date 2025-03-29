module SFEquationsMod
  ! ============================================================================
  ! Helper methods for the FATES SPITFIRE model
  ! Most equations come from:
  !  Thonicke et al. 2010, Biogeosciences 7:1991-2011
  !  Rothermel et al. 1972, Research Paper INT 115
  !  Albini et al. 1976, Research Report INT 30
  ! ============================================================================
  
  use FatesConstantsMod, only : r8 => fates_r8
  use FatesConstantsMod, only : nearzero
  
  implicit none
  private
  
  public :: MaximumReactionVelocity
  public :: OptimumReactionVelocity
  public :: OptimumPackingRatio
  public :: MoistureCoefficient
  public :: ReactionIntensity
  public :: HeatofPreignition
  public :: EffectiveHeatingNumber
  public :: WindFactor
  public :: PropagatingFlux
  public :: ForwardRateOfSpread
  public :: BackwardRateOfSpread
  public :: FireDuration
  public :: LengthToBreadth
  public :: FireSize
  public :: AreaBurnt
  public :: FireIntensity
  public :: PassiveCrownFireIntensity
  public :: HeatReleasePerArea
  public :: TorchingIndex
  public :: CrowningIndex
  public :: CrownFireIntensity
  public :: LiveFuelMoistureContent
  
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
    real(r8), parameter :: a = 0.200395_r8
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
        MaximumReactionVelocity = 1.0_r8/(0.0591_r8 + 2.926_r8*(SAV**(-1.5_r8)))
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
    
    real(r8) function ForwardRateOfSpread(bulk_density, eps, q_ig, i_r, xi, phi_wind)
      !
      !  DESCRIPTION:
      !  Calculates forward rate of spread [m/min]
      !  
      !  Flaming front of a surface fire
      !
      !  Equation 9. Thonicke et al. 2010
      !

      ! ARGUMENTS:
      real(r8), intent(in) :: bulk_density ! fulk bulk density [kg/m3]
      real(r8), intent(in) :: eps          ! effective heating number [unitless]
      real(r8), intent(in) :: q_ig         ! heat of preignition [kJ/kg] 
      real(r8), intent(in) :: i_r          ! reaction intensity [kJ/m2/min]
      real(r8), intent(in) :: xi           ! propagating flux [unitless]     
      real(r8), intent(in) :: phi_wind     ! wind factor [unitless]
    
      if ((bulk_density <= nearzero) .or. (eps <= nearzero) .or. (q_ig <= nearzero)) then
        ForwardRateOfSpread = 0.0_r8
      else
        ForwardRateOfSpread = (i_r*xi*(1.0_r8 + phi_wind))/(bulk_density*eps*q_ig)
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

    FireDuration = (SF_val_max_durat + 1.0_r8)/(1.0_r8 + SF_val_max_durat*               &
      exp(SF_val_durat_slope*FDI))

  end function FireDuration

  !---------------------------------------------------------------------------------------
 
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
  real(r8), parameter :: lb_threshold = 0.55_r8                ! tree canopy fraction below which to use grassland length-to-breadth eqn
  real(r8), parameter :: m_per_min__to__km_per_hour = 0.06_r8  ! convert wind speed from m/min to km/hr

  windspeed_km_hr = effective_windspeed/m_per_km*min_per_hr

  if (windspeed_km_hr < 1.0_r8) then 
    LengthToBreadth = 1.0_r8
 else
    if (tree_fraction > lb_threshold) then 
      LengthToBreadth = (1.0_r8 + (8.729_r8* &
        ((1.0_r8 -(exp(-0.03_r8*m_per_min__to__km_per_hour*effective_windspeed)))**2.155_r8)))
    else  
      LengthToBreadth = (1.1_r8*((m_per_min__to__km_per_hour*effective_windspeed)**0.464_r8))
    endif
  endif

  end function LengthToBreadth

  !---------------------------------------------------------------------------------------
  
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
      return
    end if
    FireSize = (pi_const/(4.0_r8*length_to_breadth))*((dist_forward + dist_back)**2.0_r8)
  
  end function FireSize
  
  !---------------------------------------------------------------------------------------
  
  real(r8) function AreaBurnt(fire_size, num_ignitions, FDI)
    !
    !  DESCRIPTION:
    !  Calculates area burnt [m2/m2/day]
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
  
  !---------------------------------------------------------------------------------------
  
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
  
  real(r8) function PassiveCrownFireIntensity(canopy_base_height, canopy_water_content)
  ! DESCRIPTION:
  ! Calculate the energy threshold for igniting crown fuels [kW/m or kJ/m/s]
  ! EQ. 11 in Scott & Reinhardt 2001

  ! ARGUMENTS:
  real(r8), intent(in) :: canopy_base_height   ! canopy base height at which biomass density > minimum density 0.011 kg/m3
  real(r8), intent(in) :: canopy_water_content ! canopy water content [%]
  

  ! Locals:
  real(r8)            :: crown_ignition_energy  ! surface fire intensity required to ignite crown fuels [kJ/kg]
  !real(r8), parameter :: canopy_water_content = 100.0_r8  ! canopy fuel water content in %, to be replaced by transient value later 

  ! Note: crown_ignition_energy to be calculated based on PFT foliar moisture content from FATES-Hydro
  ! or create foliar moisture % based on BTRAN
  ! Use foliar_moisture(currentCohort%pft) and compute weighted PFT average with Eq 3 Van Wagner 1977
  ! in place of canopy_water_content parameter

  ! Eq 3 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
  crown_ignition_energy = 460.0_r8 + 25.9_r8 * canopy_water_content

  ! Crown fuel ignition potential (kW/m), Eq 4 Van Wagner 1977, Eq 11 Scott & Reinhardt 2001
  ! FI = (Czh)**3/2 where z=canopy base height,h=heat of crown ignite energy, FI=fire intensity
  ! 0.01 = C, empirical constant Van Wagner 1977 Eq 4 for 6m canopy base height, 100% FMC, FI 2500kW/m
  ! passive_crown_FI = min fire intensity to ignite canopy fuel (kW/m or kJ/m/s)
  PassiveCrownFireIntensity = (0.01_r8 * canopy_base_height * crown_ignition_energy)**1.5_r8
end function PassiveCrownFireIntensity

!---------------------------------------------------------------------------------------
real(r8) function HeatReleasePerArea(SAV, i_r)
!
! DESCRIPTION:
! Calculate heat release per unit area of surface fire (HPA, kJ/m2)
! EQ. 2 in Scott & Reinhardt 2001

! ARGUMENTS:
real(r8), intent(in)  :: SAV      ! fuel surface area to volume ratio [/cm]
real(r8), intent(in)  :: i_r      ! reaction intensity [kJ/m2/min]


! Locals:
real(r8)   :: time_r       ! fire residence time [min]
if (SAV < nearzero) then
  time_r = 0.0_r8
else
  time_r = 12.595_r8 / SAV  ! EQ 3 in Scott & Reinhardt 2001
end if

HeatReleasePerArea = i_r * time_r

end function HeatReleasePerArea

!---------------------------------------------------------------------------------------

real(r8) function TorchingIndex(bulk_density, SAV, eps, q_ig, i_r, xi, beta_ratio, &
                                passive_crown_FI, HPA)
!
! DESCRIPTION:
! Calculate open wind speed [in km/hour] at which surface fire intensity equals 
! to the minimum fire intensity for initiating crown fire, AKA Torching Index (TI).
!  EQ. 18 in Scott & Reinhardt 2001
! XLG: currently we are ignoring the slope effect 
! ARGUMENTS:

real(r8), intent(in) :: bulk_density      ! fulk bulk density [kg/m3]
real(r8), intent(in) :: SAV               ! fuel surface area to volume ratio [/cm]
real(r8), intent(in) :: eps               ! effective heating number [unitless]
real(r8), intent(in) :: q_ig              ! heat of preignition [kJ/kg] 
real(r8), intent(in) :: i_r               ! reaction intensity [kJ/m2/min]
real(r8), intent(in) :: xi                ! propagating flux [unitless]      
real(r8), intent(in) :: beta_ratio        ! relative packing ratio [unitless]
real(r8), intent(in) :: passive_crown_FI  ! fire intensity threshold for initiating crown fire [kW/m or kJ/m/s] 
real(r8), intent(in) :: HPA               ! heat release per unit area [kJ/m2]

! Locals:
real(r8)              :: wind_coef                   ! critical wind coefficient for crown fire initiation [unitless]
real(r8)              :: b, c, e                     ! temporary variables
real(r8), parameter   :: wind_reduce_factor = 0.2_r8 ! wind reduction factor. XLG: can be incluede as a user defined param later?

! EQ. 16 in Scott & Reinhardt 2001 ignoring slope factor
wind_coef = (60._r8 * passive_crown_FI * bulk_density * eps * q_ig) / &
(HPA * xi * i_r) - 1.0_r8  

! Equation A7 in Thonicke et al. 2010 per eqn 49 from Rothermel 1972
b = 0.15988_r8*(SAV**0.54_r8)
      
! Equation A8 in Thonicke et al. 2010 per eqn 48 from Rothermel 1972 
c = 7.47_r8*(exp(-0.8711_r8*(SAV**0.55_r8)))

! Equation A9 in Thonicke et al. 2010 (appears to have typo, using coefficient Eq. 50 Rothermel 1972)
e = 0.715_r8*(exp(-0.01094_r8*SAV))

! EQ. 18 in Scott & Reinhardt 2001
TorchingIndex = (wind_coef / c * (beta_ratio**(-e)))**(1.0_r8 / b)

end function TorchingIndex

!---------------------------------------------------------------------------------------
real(r8) function CrowningIndex(eps, q_ig, i_r, canopy_bulk_density)
!
! DESCRIPTION:
! Calculate open wind speed [km/hr] at which a fully active crown fire is sustained
! EQ. 20 in Scott & Reinhardt 2001
! XLG: currently we are ignoring slope effect 
!
! ARGUMENTS:
real(r8), intent(in) :: eps                  ! effective heating number [unitless]
real(r8), intent(in) :: q_ig                 ! heat of preignition [kJ/kg] 
real(r8), intent(in) :: i_r                  ! reaction intensity [kJ/m2/min]
real(r8), intent(in) :: canopy_bulk_density  ! canopy fuel bulk density [kg biomass / m3]
! Locals:
real(r8)   :: CI_temp      ! temporary variables

if(i_r <= nearzero .or. canopy_bulk_density <= nearzero) then
  CI_temp = 0.0_r8
else
  CI_temp = (164.8_r8 * eps * q_ig) / (i_r * canopy_bulk_density) - 1.0_r8
end if
CrowningIndex = 0.0457_r8 * ((CI_temp / 0.001612_r8)**0.7_r8)

end function CrowningIndex


!---------------------------------------------------------------------------------------

real(r8) function CrownFireIntensity(HPA, canopy_fuel, CFB, ROS_final)
!
! Description
! Calculate fire intentisy for crown fire using
! EQ. 22 in Scott & Reinhardt 2001
!
! ARGUMENTS:
real(r8), intent(in) :: HPA             ! heat release per unit area [kJ/m2]
real(r8), intent(in) :: canopy_fuel     ! canopy fuel load [kg biomass] 
real(r8), intent(in) :: CFB             ! crown fraction burnt [fraction]
real(r8), intent(in) :: ROS_final       ! final rate of spread after a crown fire happens [m/min]
! Locals:
real(r8), parameter  :: H_canopy = 18000.0_r8 ! heat yield for canopy fuels [kJ/kg biomass]

CrownFireIntensity = (HPA + (canopy_fuel * H_canopy * CFB)) * &
 ROS_final / 60.0_r8 

end function CrownFireIntensity

!---------------------------------------------------------------------------------------

real(r8) function LiveFuelMoistureContent(lai, swc, max_lfmc, min_lfmc, &
swc_alpha, lai_beta, gamma_int)
!
! DESCRIPTION
! Calculates live fuel moisture content  
! using EQ. 4 in McNorton and Giuseppe 2024 'A global fuel characteristic model
! and dataset for wildfire prediction'
! LFMC = max_lfmc - min_lfmc * exp(-swc_alpha*swc + lai_beta*lai + gamma_int*swc*lai)
! swc is mass-based in the original Eq., we now switch to volumetric soil water content.
! As swc_alpha and gamma_int are fit using observations, we can make them as model 
! parameters thus this change in swc unit and resulted error hopefully can be
! adjusted through parameter tuning
!
! ARGUMENTS:
real(r8),         intent(in)    :: lai                   ! cohort level leaf area index [m2 m-2]
real(r8),         intent(in)    :: swc                   ! volumetric soil water content [m3 m-3]
real(r8),         intent(in)    :: max_lfmc              ! max. live fuel moisure content [%]
real(r8),         intent(in)    :: min_lfmc              ! min. live fuel moisture content [%]
real(r8),         intent(in)    :: swc_alpha             ! model coef. associated with swc [unitless]
real(r8),         intent(in)    :: lai_beta              ! model coef. associated with lai [unitless]
real(r8),         intent(in)    :: gamma_int             ! model coef. for interaction effect of swc and lai on LFMC [unitless]

! Locals:
real(r8)          :: effect_temp        ! the temporary variable for calculating effect of SWC and LAI

effect_temp = swc_alpha*swc + lai_beta*lai + gamma_int*swc*lai
LiveFuelMoistureContent = max_lfmc - min_lfmc*exp(-effect_temp)
      

end function LiveFuelMoistureContent


  
end module SFEquationsMod
